import glob
import pyfits
import os
import numpy as np

list=glob.glob("*.fits")
for fitsfile in list:
    try:
        base=os.path.basename(fitsfile).split('.')[0]
        f=pyfits.open(fitsfile)
        hdr=f[0].header
        spectrum=f[0].data
        wavelength=hdr['CRVAL1'] + np.arange(hdr['NAXIS1'])*hdr['CDELT1']
        out=open(base+'.txt','w')
        [out.write("%f %f\n"%(w,s)) for w,s in zip(wavelength,spectrum)]
        out.close()
    except:
        print "%s failed"%fitsfile
