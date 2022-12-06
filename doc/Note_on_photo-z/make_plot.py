from galsed import fits_utils
import numpy as np
from matplotlib import pylab as plt

#FRANZONA
fits_utils.draw_franzona(which="El_cww_fix2.txt",yscale=1.e-3)
#fits_utils.draw_franzona(which="El_B2004a.sed")
#fits_utils.draw_franzona(which="El_B2004a.sed.lsst")
#fits_utils.draw_franzona(which="El_sdss.txt")

#FORS
filename="/home/cohen/WORK/LSST/SPEC/edmond_lib/SPEC&TAB/SPEC029png&tab/SPEC029A3.out.txt"
fits_utils.draw_spectrum(filename)
filename="/home/cohen/WORK/LSST/SPEC/edmond_lib/SPEC&TAB/SPEC043png&tab/SPEC043A3sB.out.txt"
fits_utils.draw_spectrum(filename)
filename="/home/cohen/WORK/LSST/SPEC/edmond_lib/SPEC&TAB/SPEC068png&tab/SPEC058A2.out.txt"
fits_utils.draw_spectrum(filename)
filename="/home/cohen/WORK/LSST/SPEC/edmond_lib/SPEC&TAB/SPEC068png&tab/SPEC063A3.out.txt"
fits_utils.draw_spectrum(filename)
filename="/home/cohen/WORK/LSST/SPEC/edmond_lib/SPEC&TAB/SPEC082png&tab/SPEC082A2next.out.txt"
fits_utils.draw_spectrum(filename)
filename="/home/cohen/WORK/LSST/SPEC/edmond_lib/SPEC&TAB/SPEC082png&tab/SPEC082A3next.out.txt"
fits_utils.draw_spectrum(filename)
filename="/home/cohen/WORK/LSST/SPEC/edmond_lib/SPEC&TAB/SPEC09png&tab/SPEC090A3sn.out.txt"
fits_utils.draw_spectrum(filename)
filename="/home/cohen/WORK/LSST/SPEC/edmond_lib/SPEC&TAB/SPEC099png&tab/SPEC099A3sn.out.txt"
fits_utils.draw_spectrum(filename)
