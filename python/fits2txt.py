import glob
import pyfits
import os, sys
import numpy as np

default = "*.fits"
inarg = sys.argv[1]
print inarg
if os.path.exists(inarg):
    if os.path.isfile(inarg):
        default = inarg
    elif os.path.isdir(inarg):
        default = os.path.join(inarg,"*.fits")
else:
    print "wrong input: ", inarg
    exit(0)

print default

fitslist=glob.glob(default)
for fitsfile in fitslist:
    try:
        savedir = os.path.dirname(fitsfile)
        base=os.path.basename(fitsfile).split('.')[0]
        savefile = os.path.join(savedir, base+'.txt')
        f=pyfits.open(fitsfile)
        hdr=f[0].header
        spectrum=f[0].data
        wavelength=hdr['CRVAL1'] + np.arange(hdr['NAXIS1'])*hdr['CDELT1']
        out=open(savefile,'w')
        [out.write("%f %f\n"%(w,s)) for w,s in zip(wavelength,spectrum)]
        out.close()
    except:
        print "%s failed"%fitsfile
