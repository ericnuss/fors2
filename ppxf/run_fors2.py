import glob, os
from stellar_libs import setup_miles_library, load_bc03_library
import ppxf_util as util
import  numpy as np
from ppxf import ppxf, reddening_curve
import matplotlib.pyplot as plt
import pickle

c = 299792.458

def load_fors2(match):
    files =  glob.glob(os.path.join("/home/cohen/lsst/photoz/sedfitting/SPEC/edmond_lib/",'*'+match+'*'))
    t=[]
    for filename in files:
        t.append(np.loadtxt(filename, dtype={'names':('wavelength', 'flux'),'formats':('float', 'float')}))
    z=np.zeros(shape=(len(t),)) #spetra have been shifted to z=0
    return t,z, files

def load_brown(match):
    files =  glob.glob(os.path.join("/home/cohen/lsst/photoz/sedfitting/SPEC/An_Atlas_of_Galaxy_SEDs/104/379/",'*'+match+'*'))
    t=[]
    for filename in files:
        t.append(np.loadtxt(filename, dtype={'names':('wavelength', 'flux', 'obs'),'formats':('float', 'float', 'float')}))
    z=np.zeros(shape=(len(t),)) #spetra are of close galaxies
    return t,z, files

def make_plot(pp,
              #metal_grid, logAge_grid,
              title, save=None ):
    plt.figure(figsize=(20,15))
    plt.subplot(211)
    pp.plot()
    plt.title(title)
    plt.text(5500, 0.6,'Mass-weighted <logAge> [Gyr]: %.3g' %(np.sum(pp.weights*pp.grid['logAge'].ravel())/np.sum(pp.weights)))
    plt.text(5500, 0.7,'Mass-weighted <[M/H]>: %.3g' %
          (np.sum(pp.weights*pp.grid['Z'].ravel())/np.sum(pp.weights)))
    if pp.reddening is not None:
        plt.text(5500, 0.8,'Reddening E(B-V): %.3g'%pp.reddening)
    plt.subplot(212)
    s = templates.shape
    weights = pp.weights.reshape(s[1:])/pp.weights.sum()
    plt.imshow(np.rot90(weights), interpolation='nearest', 
               cmap='gist_heat', aspect='auto', origin='upper',
               extent=[pp.grid['logAge'].min(), pp.grid['logAge'].max(), pp.grid['Z'].min(), pp.grid['Z'].max()])
#               extent=[np.log10(1), np.log10(17.7828), -1.9, 0.45])
    plt.colorbar()
    plt.title("Mass Fraction")
    plt.xlabel("log$_{10}$ Age (Gyr)")
    plt.ylabel("[M/H]")
    plt.tight_layout()
    if save is not None:
        plt.savefig(save)

if __name__ == "__main__":
    TYPE = "Brown"
    MASK_MIN = 1000
    MASK_MAX = 15800
#    MASK_MIN = 3200
#    MASK_MAX = 6800
    if TYPE == "Brown":
        specs, zs, files = load_brown('Mrk_1490')
#        specs, zs, files = load_brown('Mrk_1450')
#        specs, zs, files = load_brown('CGCG_049-057')
#        specs, zs, files = load_brown('NGC_4125')
    elif TYPE == "Fors2":
        specs, zs, files = load_fors2('029Abr')
        
    for t,z,filename in zip(specs,zs, files):

        print filename
        #log rebin the galactic spectrum
        mask = True*np.ones(t['flux'].shape, dtype=bool)
        if TYPE == "Brown":
            mask=np.logical_and(t['wavelength']>MASK_MIN,
                                t['wavelength']<MASK_MAX)
        flux=t['flux'][mask]
        wave=t['wavelength'][mask]
        # else:
        #     flux=t['flux']
        #     wave=t['wavelength']
        
        velscale=np.log(wave[1]/wave[0])*c
        print 'velocity scale from the galaxy: ',velscale

        if TYPE == 'Brown':
            #we need to sample evenly the Brown spectrum
            step = np.diff(wave).min()
            x = np.arange(wave[0], wave[-1], step)
            y = np.interp(x, wave, flux)
            flux, log_wave, velscale = util.log_rebin([x[0],x[-1]], y)
        else:    
            flux, log_wave, velscale = util.log_rebin([wave[0],wave[-1]], flux)
        wave = np.exp(log_wave)
                
        galaxy = flux/np.median(flux)

        print 'velocity scale from the galaxy (after rebin): ',velscale
        #assume a constant noise level for spectrum
        noise = galaxy*0 + 0.1

        #now load the stellar library, modifying it based on above info
        #MILES
        # templates, lamRange_temp, logAge_grid, metal_grid = \
        # setup_miles_library(velscale, FWHM_gal)
        #BC03
        if TYPE == 'Brown':
            lam, full, full_data, logAge_grid, metal_grid, templates = \
            load_bc03_library(x, velscale*(0.9999), FWHM_gal=1, version='S')
        else:
            lam, full, full_data, logAge_grid, metal_grid, templates = \
            load_bc03_library(t['wavelength'][mask], velscale, FWHM_gal=1, version='N')
            
        lamRange_temp = [wave[0],wave[-1]]
        
        #normalize templates by median scalar
        norm=np.median(templates)
        templates /= norm
        full /= norm
        
        dv = c*np.log(lamRange_temp[0]/wave[0])  # km/s
        
        goodpixels = util.determine_goodpixels(log_wave, lamRange_temp, z)

        start = [c*np.log(1 + z), 0.01]#3*velscale]
        fixed = [True, True]#, False, False]
        print start
        
        pp = ppxf(templates, galaxy, noise, velscale, start,
          goodpixels=goodpixels, plot=False, moments=2, degree=-1,
          vsyst=dv, clean=False, regul=1,
          mdegree=0, reddening=0.1,lam=wave,
          fixed=fixed
          )

        # noise *= np.sqrt(pp.chi2)
        # pp = ppxf(templates, galaxy, noise, velscale, start,
        #           goodpixels=goodpixels, plot=False, moments=2, degree=-1,
        #           vsyst=dv, clean=False, regul=0.,
        #           mdegree=0, reddening=None,lam=wave
        #           )
        
        pp.grid={'Z':metal_grid,'logAge':logAge_grid}
        title=os.path.basename(filename).split('.')[0]
        outdir='output_bc03s'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        make_plot(pp, title=title,
                  save=os.path.join(outdir,title+'.png')
                  )

        print('Desired Delta Chi^2: %.4g' % np.sqrt(2*goodpixels.size))
        print('Current Delta Chi^2: %.4g' % ((pp.chi2 - 1)*goodpixels.size))

        
        fig = plt.figure(figsize=(20,15))
        ax2 = fig.add_subplot(111)#figure(figsize=(20,15))
        pop_spectrum = np.dot(full.reshape((full.shape[0],full.shape[1]*full.shape[2])),pp.weights)
        ax2.semilogy(lam,pop_spectrum, label='weighted star pop')
        red_curve = reddening_curve(lam, pp.reddening)
        reddened = pop_spectrum*red_curve
        ax2.semilogy(lam, reddened, label='reddened star pop.')
        ax2.semilogy(wave, pp.bestfit, label='bestfit')
        # red_curve = reddening_curve(wave, pp.reddening)
        # dereddened = pp.bestfit/red_curve
        # ax2.semilogx(wave, dereddened, label='dereddened fit')
        ax2.semilogy(t['wavelength'], t['flux']/np.median(flux), label='gal')
        ax2.legend()
        # ax2.plot([MASK_MIN,MASK_MIN],ax2.get_ylim(), 'r')
        # ax2.plot([MASK_MAX,MASK_MAX],ax2.get_ylim(), 'r')
        ax2.set_xlim(MASK_MIN,MASK_MAX)
        ax2.set_ylim(0.01, 2*pop_spectrum.max())
        plt.savefig(os.path.join(outdir,title+'_full.png'))
                
        pickle.dump(pp,open(os.path.join(outdir,title+'.pkl'), 'w'))
        

        # print('Mass-weighted <logAge> [Gyr]: %.3g' %
        #       (np.sum(pp.weights*logAge_grid.ravel())/np.sum(pp.weights)))
        # print('Mass-weighted <[M/H]>: %.3g' %
        #       (np.sum(pp.weights*metal_grid.ravel())/np.sum(pp.weights)))

    plt.show()

    # templates convolues et attenues:
    # plot(np.dot(pp.matrix,pp.weights))

    # templates de depart (need reshape)
    #plot(np.dot(templates,pp.weights))

    # galaxy apres log rebin
    # plot(galaxy)
