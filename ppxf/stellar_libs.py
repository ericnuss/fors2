import os, glob 
import astropy.io.fits as pyfits
from astropy.io import ascii
from matplotlib import pylab as plt
import numpy as np
import ppxf_util as util
from scipy import ndimage

basedir="/home/cohen/lsst/photoz/sedfitting"

def load_bc03_library(wave=None, velscale=None, FWHM_gal=1, version='N'):
    if version=="S":
        basefile=os.path.join(basedir,'STARLIGHT/STARLIGHTv04/Base.BC03.S')
        nAges = 25
        nMetal = 6
    else:
        basefile=os.path.join(basedir,'STARLIGHT/STARLIGHTv04/Base.BC03.N')
        nAges = 15
        nMetal = 3

    logAge_grid = np.empty((nAges, nMetal))
    metal_grid = np.empty((nAges, nMetal))

    create=True
    data = ascii.read(basefile, names=['specfile','age','Z','code','Mstar','yav','a/Fe'], fast_reader=False, data_start=1)
    for i, d in enumerate(data):
        tmp=os.path.join(os.path.dirname(basefile),'BasesDir', d['specfile'])
        lam, ssp = np.loadtxt(tmp, unpack=True)
        if create:
            full = np.empty((lam.size, nAges, nMetal))
            full_data = np.empty((nAges, nMetal), dtype=object)
            if wave is not None:
                templates = np.empty((wave.size, nAges, nMetal))
                FWHM_tem = 1 # this number is crap as the spectra are synthetic in part...
                FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                sigma = FWHM_dif/2.355/np.diff(lam)[0] # Sigma difference in pixels

            create=False
        full[:,i%nAges,int(i/nAges)] = ssp
        full_data[i%nAges,int(i/nAges)] = d
        logAge_grid[i%nAges,int(i/nAges)] = np.log10(d['age']/1.e9)
        metal_grid[i%nAges,int(i/nAges)] = np.log10(d['Z'])
        if wave is not None:
            templates[:,i%nAges,int(i/nAges)] = \
              build_templates(lam, ssp, wave, velscale, sigma)
    #lam is the same for all spectra
    if wave is None:
        return lam, full, full_data, logAge_grid, metal_grid
    else:
        return lam, full, full_data, logAge_grid, metal_grid, templates
    
def build_templates(lam, ssp, wave, velscale, sigma):
        new_ssp = np.interp(wave, lam, ssp, right=0, left=0)
        lam=wave
        ssp=new_ssp
        ssp = ndimage.gaussian_filter1d(new_ssp,sigma)
        sspNew, logLam2, velscale = util.log_rebin([wave[0],wave[-1]], ssp, velscale=velscale)
 
        return sspNew

    
def setup_bc03_library(wave, velscale=None, FWHM_gal=1, version='N', in_log=True):
    FWHM_tem = 1 # this number is crap as the spectra are synthetic in part...
    lam, full, full_data, logAge_grid, metal_grid = load_bc03_library(version=version)

    templates = np.zeros_like(full)
    
    new_ssp = np.interp(wave,lam, ssp, right=0, left=0)
    lamRange_temp = [wave[0],wave[-1]]
    delta = np.diff(lam)[0]
    if in_log:
        sspNew, logLam2, velscale = util.log_rebin(lamRange_temp, ssp, velscale=0.99999*velscale)
    else:
        sspNew=ssp
        logLam2=np.log(lam)
    #

    templates = np.empty((sspNew.size, nAges, nMetal))

    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
    sigma = FWHM_dif/2.355/delta # Sigma difference in pixels

    for i,d in enumerate(data):
        tmp=os.path.join(os.path.dirname(basefile),'BasesDir',d['specfile'])
        lam,ssp = np.loadtxt(tmp, unpack=True)
        full[:,i%nAges,int(i/nAges)] = ssp
        new_lam = wave
        #new_lam = np.arange(new_lam[0], new_lam[-1],1)
        new_ssp = np.interp(new_lam,lam, ssp, right=0, left=0)
        # new_lam = np.arange(lam[0], lam[-1],1)
        # new_ssp = np.interp(new_lam,lam, ssp, right=0, left=0)
        # mask = (new_lam > wave.min()-1) & (new_lam < wave.max()+2)
        lam=new_lam#[mask]
        ssp=new_ssp#[mask]
        ssp = ndimage.gaussian_filter1d(ssp,sigma)
        if in_log:
            sspNew, logLam2, velscale = util.log_rebin(lamRange_temp, ssp, velscale=velscale)
        else:
            sspNew=ssp
        templates[:,i%nAges,int(i/nAges)] = sspNew
        logAge_grid[i%nAges,int(i/nAges)] = np.log10(d['age']/1.e9)
        metal_grid[i%nAges,int(i/nAges)] = np.log10(d['Z'])
            
    return templates, lamRange_temp, logAge_grid, metal_grid#, lam_full, full

#------------------------------------------------------------------------------
def setup_miles_library(velscale, FWHM_gal):

    # Read the list of filenames from the Single Stellar Population library
    # by Vazdekis et al. (2010, MNRAS, 404, 1639) http://miles.iac.es/.
    #
    # For this example I downloaded from the above website a set of
    # model spectra with default linear sampling of 0.9A/pix and default
    # spectral resolution of FWHM=2.51A. I selected a Salpeter IMF
    # (slope 1.30) and a range of population parameters:
    #
    #     [M/H] = [-1.71, -1.31, -0.71, -0.40, 0.00, 0.22]
    #     Age = np.linspace(np.log10(1), np.log10(17.7828), 26)
    #
    # This leads to a set of 156 model spectra with the file names like
    #
    #     Mun1.30Zm0.40T03.9811.fits
    #
    # IMPORTANT: the selected models form a rectangular grid in [M/H]
    # and Age: for each Age the spectra sample the same set of [M/H].
    #
    # We assume below that the model spectra have been placed in the
    # directory "miles_models" under the current directory.
    #
    vazdekis = glob.glob(os.path.join(basedir, 'ppxf','miles_models/Mun1.30*.fits'))
    vazdekis.sort()
    FWHM_tem = 2.51 # Vazdekis+10 spectra have a resolution FWHM of 2.51A.

    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the SDSS galaxy spectrum, to determine
    # the size needed for the array which will contain the template spectra.
    #
    hdu = pyfits.open(vazdekis[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lamRange_temp = h2['CRVAL1'] + np.array([0., h2['CDELT1']*(h2['NAXIS1']-1)])
    sspNew, logLam2, velscale = util.log_rebin(lamRange_temp, ssp, velscale=velscale)

    # Create a three dimensional array to store the
    # two dimensional grid of model spectra
    #
    nAges = 26
    nMetal = 6
    templates = np.empty((sspNew.size, nAges, nMetal))

    # Convolve the whole Vazdekis library of spectral templates
    # with the quadratic difference between the SDSS and the
    # Vazdekis instrumental resolution. Logarithmically rebin
    # and store each template as a column in the array TEMPLATES.

    # Quadratic sigma difference in pixels Vazdekis --> SDSS
    # The formula below is rigorously valid if the shapes of the
    # instrumental spectral profiles are well approximated by Gaussians.
    #
    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
    sigma = FWHM_dif/2.355/h2['CDELT1'] # Sigma difference in pixels

    # These are the array where we want to store
    # the characteristics of each SSP model
    #
    logAge_grid = np.empty((nAges, nMetal))
    metal_grid = np.empty((nAges, nMetal))

    # These are the characteristics of the adopted rectangular grid of SSP models
    #
    logAge = np.linspace(np.log10(1), np.log10(17.7828), nAges)
    metal = [-1.71, -1.31, -0.71, -0.40, 0.00, 0.22]

    # Here we make sure the spectra are sorted in both [M/H]
    # and Age along the two axes of the rectangular grid of templates.
    # A simple alphabetical ordering of Vazdekis's naming convention
    # does not sort the files by [M/H], so we do it explicitly below
    #
    metal_str = ['m1.71', 'm1.31', 'm0.71', 'm0.40', 'p0.00', 'p0.22']
    for k, mh in enumerate(metal_str):
        files = [s for s in vazdekis if mh in s]
        for j, filename in enumerate(files):
            hdu = pyfits.open(filename)
            ssp = hdu[0].data
            ssp = ndimage.gaussian_filter1d(ssp,sigma)
            sspNew, logLam2, velscale = util.log_rebin(lamRange_temp, ssp, velscale=velscale)
            templates[:, j, k] = sspNew  # Templates are *not* normalized here
            logAge_grid[j, k] = logAge[j]
            metal_grid[j, k] = metal[k]

    return templates, lamRange_temp, logAge_grid, metal_grid

#------------------------------------------------------------------------------

def plot_miles():
    miles_dir=os.path.join(basedir,"ppxf/miles_models/")
    miles_list=glob.glob(os.path.join(miles_dir,'Mun1*'))
    for star in miles_list:
        f=pyfits.open(star)
        h=f[0].header
        lam=h['CRVAL1']+np.arange(h['NAXIS1'])*h['CDELT1']
        plt.plot(lam, f[0].data)


# starlight_dir=os.path.join(basedir,"SPEC/STARLIGHT/STARLIGHTv04/BasesDir/")
# starlight_list=glob.glob(os.path.join(starlight_dir,'bc2003*'))
# for star in starlight_list:
#     x,y=np.loadtxt(star, unpack=True)
#     plt.loglog(x,y)

# plt.show()

