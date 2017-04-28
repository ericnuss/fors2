import numpy as np

# import matplotlib
# matplotlib.use('GTKAgg',warn=False, force=True)
from matplotlib import pylab as plt

from exceptions import *
from astropy.io import fits
import os, glob

defined_lines={
    'Ly{beta}':1025, 'Ly{alpha}':1215, 'CIV':1549, 'AlIII':1859, 'CIII':1909,
    'CII':2326, 'FeII24':[2344.2,2374.5,2382.8], 'FeII26':[2586.5,2600.2], 'MgII':[2796.4,2803],
    'AlII':3584, 'K':3933.7, '[OII]':3728.2, 'H11':3769.7, 'H10':3797.9, 'H9':3835,
    'NeIII':[3869,3968], 'H8':3889, 'H':3968.5, 'H{epsilon}':3970.07,
    '4000{AA}-break':4000, 'H{delta}':4101.73, 'G':4303, 'H{gamma}':4340.46, 'Fe43':4383, 'Ca44':4455,
    'Fe45':4531, 'NIII':[4634,4640], 'H{beta}':4861.32, '[OIII]':[4958.9,5006.8],
    'HeI':[4388,4471,6,4713,3,4922,5015,5048,5876,6678,7065.], 'Mgb':5175.3, 'E':5270, 'NaD':5892,
    'H{alpha}':6563, 'NII':6584, 'S2':6716}
pretty_labels={'H{alpha}':'H$_\\alpha$','H{beta}':'H$_\\beta$','H{delta}':'H$_\\delta$', 'H{epsilon}':'H$_\\epsilon$', 'H{gamma}':'H$_\\gamma$'}


def plot_spectrum(filename, show=True, save=False, savedir='.'):
    x,y,m = np.loadtxt(filename, unpack=True)
    yy = np.ma.masked_array(y,m)
    comments = open(filename,'r').readlines()[:5]
    obs, slit, z, lines, spec = [ll.split(' ')[-1].strip('\n') for ll in comments]
    lines_list=lines.split(",")
    pretty_plot(x, y, yy, spec, float(z), lines_list, show, save, savedir)

    
def pretty_plot(x, y, yy, spec, z, lines, show=True, save=False, savedir='.'):
    fig=plt.figure(figsize=(20,9))
    plt.plot(x,y)
    plt.plot(x,yy)
    plt.title('SPEC%sn.txt'%spec)
    m=yy.min()
    M=yy.max()
    if m<0:
        m=1.1*m
    else:
        m=0.9*m
    plt.ylim(m, 1.1*M)
    
    for l in lines:
        if l in ['?','weak','broad','(QSO)','QSO','',' ']:
            continue
        val=0.
        if l in defined_lines:
            #turn all into array, in case of line multiplets
            vals=np.array(defined_lines[l])*(1.+z)
            if vals.shape==():
                vals=vals[np.newaxis]
            for val in vals:
                plt.plot([val,val],[m,M],'--')
                if l in pretty_labels:
                    l=pretty_labels[l]
                plt.text(val,M,l)
        else:
            print "%s not found"%l


    if save==True:
        imgdir = os.path.join(savedir,'IMG')
        if not os.path.exists(imgdir):
            os.mkdir(imgdir)
        imgfile = os.path.join(imgdir,'IMG%sn.png'%spec)
        plt.savefig(imgfile)
    if show==True:
        plt.show()
    else:
        plt.close(fig)
    
def get_catalog_info(spec, cat):
    try:
        spec = int(spec)
        if spec in cat.data['ID']:
            catid=(cat.data['ID']==spec)
            z=cat.data[catid]['z'][0]
            lines=cat.data[catid]['Lines'][0]
        else:
             z=-1
             lines='not in catalog'
    except:
        z=-1
        lines='redshift unknown'
    return z,lines

def do_normalize(f,out=None):
    basedir = os.path.dirname(f)
    hdlt=glob.glob(os.path.join(basedir,'HDLT*.txt'))
    if len(hdlt)>1:
        raise Exception('Several HDLT files found!')
    elif len(hdlt)==0:
        hdlt=glob.glob(os.path.join(basedir,'HD*LT*.txt'))
        if len(hdlt)==0:
            raise Exception('HDLT file not found in %s'%basedir)
    hdlt=hdlt[0]
    x,y=np.loadtxt(f,unpack=True)
    xx,ccd=np.loadtxt(hdlt, unpack=True)
    #align wavelength  of HDLT and spectrum
    if len(x)>len(xx):
        try:
            maxid = np.where(x==xx[-1])[0][0]
        except IndexError:
            maxid = np.where(np.isclose(x,xx[-1]))[0][0]
        x=x[:maxid+1]
        y=y[:maxid+1]
    #remove zeron entries of the ccd
    mask=ccd>0
    xx=xx[mask]
    ccd=ccd[mask]
    x=x[mask]
    y=y[mask]
    if not np.allclose(xx,x,1.e-3):
        raise Exception('HDLT.txt has a different wavelength array as the spectrum file %s'%f)
    slit=f.split('SPEC')[1][:2]
    if out is not None:
        np.savetxt(out, np.c_[x,y/ccd], fmt="%.6f")
    return x, y/ccd

def run_on_spec(obs, slit, spec, datadir, show=False, save=True, savedir='.'):
        #get catalog info and discard stars (negative redshifts)
        z,lines=get_catalog_info(spec,cat)
        if z<=0:
            return [], [], z, lines
        
        filename=os.path.join(datadir,"%s"%obs,"SPEC%sn.txt"%slit)

        if not os.path.exists(filename):
            unormalized=os.path.join(datadir,"%s"%obs,"SPEC%s.txt"%slit)
            if os.path.exists(unormalized):
                try:
                    x,y = do_normalize(unormalized)
                except Exception,e:
                    print e, 'for slit %s'%slit
                    return [], [], z, lines
            else:
                print "%s does not exists, not creating %s"%(unormalized, 'SPEC%sn.txt'%spec)
                return [], [], z, lines
        else:
            x,y=np.loadtxt(filename,unpack=True)

        yy=np.ma.masked_equal(y,0)

        try:
            if bands[spec][0]==-1:
                #disregard the obs ban masking.
                spect_bands=np.array(bands[spec])[1:]
            else:
                spect_bands = np.array(bands[obs]+bands[spec])
        except KeyError:
            try:
                spect_bands = np.array(bands[obs])
            except KeyError:
                spect_bands = np.array([])

        if len(spect_bands)!=0:
            for start, stop in zip(spect_bands[::2],spect_bands[1::2]):
                yy.mask[np.where((stop>x)&(x>start))[0]]=True

        #saving:
        np.savetxt(os.path.join(new_dir,'SPEC%sn.txt'%spec), np.c_[x,yy,yy.mask.ravel()], fmt=("%.6f","%.6f","%d"), header="obs: %s\nslit: %s\nz: %.5f\nlines: %s\nspec: %s"%(str(obs),slit,z,lines, str(spec)))

        #plotting
        if yy.shape!=yy.mask.shape:
            yy.mask = yy.mask.ravel()
        lines_list=lines.split(",")
        pretty_plot(x,y,yy,spec,z,lines_list, show, save, savedir)
        
        return x, yy, z, lines_list
        
def run_all(obs, datadir, savedir, show=False, save=True):
    
    slit_ids = obs_dict[obs]['slit']
    spec_ids = obs_dict[obs]['spec']
    
    if len(slit_ids)!=len(spec_ids):
        raise Exception('Something is wrong! Bailing out')
        return

    for i in range(len(slit_ids)):
        slit=str(slit_ids[i])
        spec=str(spec_ids[i])
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)


if __name__=="__main__":
    import yaml, sys
    config_file = 'config.yaml'
    map_table = 'spec_table.yaml'
    if len(sys.argv)>1:
        config_file = sys.argv[1]
    if len(sys.argv)>2:
        map_table = sys.argv[2]
    
    
    #load the config yaml:
    config = yaml.load(open(config_file,'r'))
    datadir=config['datadir']
    workdir=config['workdir']
    catfile = config['fors2catalog']
    cat = fits.open(catfile)[1]
    
    #load the master table of fors2 identifiers : obs/slit_id vs spec_id mappings
    maptable = yaml.load(open(map_table,'r'))
    bands = maptable['mask_bands']
    obs_dict = maptable['fors2_obs']

    obs_list = obs_dict.keys()
    if len(sys.argv)>3:
        obs_list = sys.argv[3].split(",")
        
    print("config: %s\nmapping: %s\nobservatios: %s\n"%(config_file,map_table,obs_list))
        
    for obs in obs_list:
        print "Running on ",obs
        new_dir=os.path.join(workdir,obs)
        if not os.path.exists(new_dir):
            os.mkdir(new_dir)
        run_all(obs, datadir, new_dir, show=False, save=True)
        
