import os
from astropy.io import ascii
import numpy as np
import ppxf_util as util
from scipy import ndimage

def check_bc03file(full_path):
    """
    Check the consistency of the selection file.
    Input :
    - full_path : full path to the selection file, that lists the
    desired seds to use
    Output :
    - data : astropy record of the selection file
    - nAges, nMetal : # of ages and metalicities selected
    - sampling : size of the spectra, 6900 for BC03
    - track : either padova93 or padova00, which determines 
    the metalicity mapping
    """
    basename = os.path.basename(full_path)
    if basename in ['Base.BC03.N',
                        'Base.BC03.S',
                        'Base.bc03.Padova1994.chab.All']:
        #bypass checks if the files are the known 2
        #starlight examples, or the full file :
        #these are subsamples from Padova94 track
        if 'N' in basename:
            nAges=15
            nMetal=3
        elif 'S' in basename:
            nAges = 25
            nMetal = 6
        elif 'All' in basename:
            nAges = 221
            nMetal = 6            
        data = ascii.read(full_path, 
                names=['specfile','age','Z','code','Mstar','yav','a/Fe'], 
                    fast_reader=False, data_start=1)
        sampling = 6900
        track = 'padova94'
    else:
        #need to code more general tests
        pass
    return data, nAges, nMetal, sampling, track


def z2metal(z, track):
    """
    Return the metalicity of a given Z, 
    for track padova93 or padova00
    """
    metal_pd94 ={'0.0001':-2.2490, '0.0004':-1.6464, '0.004':-0.6392, '0.008':-0.3300, '0.02':0.0932, '0.05':0.5595}
    metal_pd00 ={'0.0004':-1.6469, '0.001':-1.2486, '0.004':-0.6392, '0.008':-0.3300, '0.019':0.0660, '0.03':0.2883}
    if track == "padova94":
        return metal_pd94[str(z)]
    else:
        return metal_pd00[str(z)]
    


def load_bc03(bc03_path, selection_file):
    """
    Load a set of templates from the Bruzual&Charlot 2003 work.
    Input :
    - bc03_path : path of the directory containing the BC03 files 
    (typically a directory named BasesDir if directly obtained 
    by untarring one of the BC03 distributions from the STARLIGHT 
    site)
    - selection_file : file listing the subset of files that are wanted
    Output : 
    - bc03_lambda : the array of wavelength in the BC03 templates 
    (identical for all files)
    - bc03_sed : grid of templates
    dimensions #of points x #of ages x #of metalicities
    - logAge_grid : grid of ages
    dimension  #of ages x #of metalicities
    - metal_grid : grid of metalicities
    dimension  #of ages x #of metalicities
    """
    data, nAges, nMetal, sampling, track = check_bc03file(selection_file)
    logAge_grid = np.empty((nAges, nMetal))
    metal_grid = np.empty((nAges, nMetal))
    bc03_seds = np.empty((sampling, nAges, nMetal))
    
    for i, d in enumerate(data):
        tmp=os.path.join(bc03_path, d['specfile'])
        lam, ssp = np.loadtxt(tmp, unpack=True)
        bc03_seds[:,i%nAges,int(i/nAges)] = ssp
        logAge_grid[i%nAges,int(i/nAges)] = np.log10(d['age']/1.e9)
        metal_grid[i%nAges,int(i/nAges)] = (z2metal(d['Z'], track))
        
    bc03_lambda = lam #take the last one, as they are all identical
    return bc03_lambda, bc03_seds, logAge_grid, metal_grid



def smooth_and_rebin(bc03, bc03_lambda, gal_lambda, velscale, FWHM_gal=1, FWHM_tem=0):
    """
    Smooth and rebin the BC03 spectra based on the sampling 
    of the galaxy to be fitted.
    Input :
    - bc03 : the spectral templates
    - bc03_lambda : the template wavelength array
    - gal_lambda : the galaxy wavelength array, which the rebinning 
    is going to match. Note tha the galaxy must have already 
    been log-rebinned. 
    - velscale : velocity scale of the galaxy, for log rebinning
    - FWHM_gal : FWHM of the galaxy resolution
    - FWHM_tem : FWHM of the template resolution
    Output :
    - templates : the smoothed and rebinned version of the BC03 input
    """
    templates = np.empty(((gal_lambda.size,)+(bc03.shape[1],bc03.shape[2])))
    for i in range(bc03.shape[1]):
        for j in range(bc03.shape[2]):
            ssp = bc03[:,i,j]
            #resize the templates to the log-rebinned galaxy wavelengths
            new_ssp = np.interp(gal_lambda, bc03_lambda, ssp, right=0, left=0)
            #log-rebin the interpolated templates
            #the output velscape should be close to identical to the input
            #Likewise, logLam2 == np.log(gal_lambda) to about 1e-5 
            sspNew, logLam2, velscale = util.log_rebin([gal_lambda[0],gal_lambda[-1]], new_ssp, velscale=velscale)
            #smooth the log-rebinned templates
            FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
            sigma = FWHM_dif/2.355/np.diff(gal_lambda)[0]
            smoothed_ssp = ndimage.gaussian_filter1d(new_ssp, sigma)
            templates[:, i, j] = smoothed_ssp
    return templates


def plot_weights(weights, logAge_grid, metal_grid, nodots=False, colorbar=True, **kwargs):
    xgrid = logAge_grid
    ygrid = metal_grid
    util.plot_weights_2d(xgrid, ygrid, weights,
                             nodots=nodots, colorbar=colorbar, **kwargs)
