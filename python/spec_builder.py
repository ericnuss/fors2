import numpy as np

# import matplotlib
# matplotlib.use('GTKAgg',warn=False, force=True)
from matplotlib import pylab as plt

from exceptions import *
from astropy.io import fits
import os, glob
from operator import or_

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

def merge_spectra(spec, files, savedir):
    obs_weight = {"255706":1, "255708":1, "255712":1, "255714":1, "255716":1, "255718":1, "255720-25fits":5, "255730-33":1, "255734fits":1, "255754-58fits":1, "255750":1, "255763":1, "255748":1, "255761":1, "255738":1, "255739-43":4, "255752-59-60":3, "255744-45-47":3}
    outfile = os.path.join(savedir,'SPEC%sn.txt'%spec)
    spec_count = 0
    sum_yy=None
    all_list=[]
    for f in files:
        x, y, m, yy, comments = open_spec_file(f)
        print x[0], x[-1],len(x), len(y), len(m), len(yy)
        obs, slit, z, lines_list, spec = parse_comments(comments)
        weight = obs_weight[obs]
        spec_count += weight
        if sum_yy is None:
            sum_yy = yy
        else:
            sum_yy = np.ma.array(sum_yy.data+yy.data,mask=map(or_,sum_yy.mask,yy.mask))
        all_list += lines_list
    all_list = np.array(all_list).unique()
    sum_yy.data = sum_yy.data/spec_count
    np.savetxt(outfile, np.c_[x,sum_yy,sum_yy.mask.ravel()], fmt=("%.6f","%.6f","%d"), header="obs: %s\nslit: %s\nz: %.5f\nlines: %s\nspec: %s"%('merged','NA',z,lines, str(spec)))     

def run_merge_spectra(input_dir='2*/SPEC*',savedir='./merged'):
    if not os.path.exists(savedir):
        os.mkdir(savedir)
    d=glob.glob(input_dir)
    t=[d[i].split('SPEC')[1].strip('n.txt') for i in range(len(d))]
    u=np.unique(t)
    for uu in u:
        indices = np.where(np.array(t)==uu)[0]
        if len(indices)==1:
            continue
        print(uu, indices)
        #merge_spectra(uu, [d[idx] for idx in indices], savedir)
        fig=plt.figure(figsize=(20,9))
        xmin=0
        xmax=0
        s=""
        for idx in indices:
            x, y, m, yy, comments = open_spec_file(d[idx])
            if xmin == 0:
                xsave=np.ma.masked_array(x,m)
                ysave=yy
                xmin=xsave[0]
                xmax=xsave[-1]
                opt=1.
            else:
                x=np.ma.masked_array(x,m)
                y=yy
                # xmin=np.max([xmin,x[0]])
                # xmax=np.min([xmax,x[-1]])
                # mask=np.logical_and(xmin<xsave,xsave<xmax)
                # xx1 = xsave[mask]
                # yy1 = ysave[mask]
                # mask=np.logical_and(xmin<x,x<xmax)
                # xx2 = x[mask]
                # yy2 = y[mask]
                # y2new=np.interp(xx1,xx2,yy2)
                # opt=(y2new*yy1).sum()/(yy1*yy1).sum()
                ynew=np.interp(xsave,x,y)
                opt=(ynew*ysave).sum()/(ysave*ysave).sum()
            plt.plot(x,yy/opt, label=d[idx])
            s += " %s : %f, "%(d[idx],opt)
        print s
        plt.title('SPEC%sn.txt'%uu)
        plt.legend()
        plt.savefig('merged/SPEC%sn.png'%uu)
        
def open_spec_file(filename):
    x, y, m = np.loadtxt(filename, unpack=True)
    yy = np.ma.masked_array(y,m)
    comments = open(filename,'r').readlines()[:5]
    return x, y, m, yy, comments

def parse_comments(comments):
    obs, slit, z, lines, spec = [ll.split(' ')[-1].strip('\n') for ll in comments]
    lines_list=lines.split(",")
    return obs, slit, z, lines_list, spec

def plot_spectrum(filename, show=True, save=False, savedir='.'):
    x, y, m, yy, comments = open_spec_file(filename)
    obs, slit, z, lines_list, spec = parse_comments(comments)
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
        if l in ['?','weak','broad','(QSO)','QSO','',' ','not in catalog']:
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
            print("%s not found"%l)


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
    outstr = "SPEC %s: "%str(spec)
    #get catalog info and discard stars (negative redshifts)
    z,lines=get_catalog_info(spec,cat)
    if z<=0:
        outstr += "z<=0 in cat, with lines %s. Maybe a star?"%lines
        print outstr
        #return [], [], z, lines

    filename=os.path.join(datadir,"%s"%obs,"SPEC%sn.txt"%slit)

    if not os.path.exists(filename):
        unormalized=os.path.join(datadir,"%s"%obs,"SPEC%s.txt"%slit)
        if os.path.exists(unormalized):
            try:
                x,y = do_normalize(unormalized)
            except Exception,e:
                print(e, 'for slit %s'%slit)
                return [], [], z, lines
        else:
            outstr += "%s does not exists, not creating %s"%(unormalized, 'SPEC%sn.txt'%spec)
            print outstr
            return [], [], z, lines
    else:
        x,y=np.loadtxt(filename,unpack=True)
        
    x, yy, lines_list = mask_and_save(x,y,obs,spec,z,lines,show, save, savedir)
    return x, yy, z, lines_list

def mask_and_save(x,y,obs,spec,z,lines,show=False, save=True, savedir='.'):
    yy = np.ma.MaskedArray(y, mask=(y==0))
    # yy=np.ma.masked_values(y,0)
    # if np.isscalar(yy.mask):
    #     yy.mask=np.zeros(y.shape,dtype=bool)
    spec=str(spec)
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
            x=np.array(x)
            yy.mask[np.where((stop>x)&(x>start))[0]]=True
    #saving:
    np.savetxt(os.path.join(savedir,'SPEC%sn.txt'%spec), np.c_[x,yy,yy.mask.ravel()], fmt=("%.6f","%.6f","%d"), header="obs: %s\nslit: %s\nz: %.5f\nlines: %s\nspec: %s"%(str(obs),slit,z,lines, str(spec)))

    #plotting
    if yy.shape!=yy.mask.shape:
        yy.mask = yy.mask.ravel()
    lines_list=lines.split(",")
    pretty_plot(x,y,yy,spec,z,lines_list, show, save, savedir)
    return x, yy, lines_list

def weighted_average(a,b,w1,w2):
    c = (a*b).sum()/(b*b).sum()
    return ( a*w1 + c * b*w2 )/(w1+w2), c
    
def average(x1,y1,x2,y2,w1,w2):
    #print len(x1),len(x2)
    print x1,x2,len(x1),len(x2)
    if len(x1)==len(x2):
        if np.all(np.isclose(x1,x2)):
            yc,coeff = weighted_average(y1,y2,w1,w2)
            return x1,yc, coeff

    i,x_0,x_1=(1,x1[0],x2[0]) if x1[0]<=x2[0] else (2,x2[0],x1[0])
    j,x_2,x_3=(1,x1[-1],x2[-1]) if x1[-1]<=x2[-1] else (2,x2[-1],x1[-1])
    m1=np.logical_and(x1>=x_1, x1<=x_2)
    m2=np.logical_and(x2>=x_1, x2<=x_2)

    xx1=x1[m1]
    xx2=x2[m2]
    yy1=y1[m1]
    yy2=y2[m2]
    x1x2=np.concatenate((xx1,xx2))
    if len(np.unique(x1x2))==len(x1x2):
        yy2 = np.interp(xx1,xx2,yy2)
        yc, coeff = weighted_average(yy1, yy2,w1,w2)
        xc = xx1
    else:
        print "recursing..."
        xc,yc, coeff = average(xx1,yy1,xx2,yy2,w1,w2)
        
    xg = x1[x1<x_1] if i==1 else x2[x2<x_1]
    if i==1:
        print 'toto1'
        idx=np.where(x1<x_1)[0][-1]
        yg = y1[x1<x_1] *yc[0:10].mean()/y1[idx]
        print yc[0:10].mean(),y1[idx], y1[x1<x_1].mean()
    else:
        print 'toto2'
        idx=np.where(x2<x_1)[0][-1]
        yg = y2[x2<x_1]*yc[0:10].mean()/y2[idx]

    xd = x2[x2>x_2] if j==1 else x1[x1>x_2]
    if j==1:
        print 'toto3'
        idx=np.where(x2>x_2)[0][0]
        yd = y2[x2>x_2]*yc[-10:-1].mean()/y2[idx]
    else:
        print 'toto4'
        idx = np.where(x1>x_2)[0][0]
        yd = y1[x1>x_2]*yc[-10:-1].mean()/y1[idx]
    
    if len(xd)==0 and len(yd)==0:
        return xc, yc, None
    elif len(xd)==0:
        return np.concatenate((xg,xc)), np.concatenate((yg,yc)), None
    elif len(xg)==0:
        return np.concatenate((xc,xd)), np.concatenate((yc,yd)), None
    else:
        return np.concatenate((xg,xc,xd)), np.concatenate((yg,yc,yd)), None
        
def merge_files(filenames, obses, weights, spec,z,lines,show=False, save=True, savedir='.', debug=True):
    outstr = "SPEC %s: "%str(spec)
    if len(filenames)>2:
        #outstr +=  'merging '
        #print outstr,filenames[:2]
        while len(filenames)>=2:
            merge_files(filenames[:2], obses[:2],weights[:2], spec,z,lines,show, save, savedir, debug)
            filenames=[os.path.join(savedir,"SPEC%sn.txt"%str(spec))]+filenames[2:]
        return
    else:
        outstr += 'merging '
        print outstr, filenames
        try:
            x1, y01, m1, y1, c1 = open_spec_file(filenames[0])
        except ValueError:
            x,y = np.loadtxt(filenames[0], unpack=True)
            x1, y1, l1 = mask_and_save(x,y,obses[0],spec,z,lines,show,save,savedir)
        try:
            x2, y02, m2, y2, c2 = open_spec_file(filenames[1])
        except ValueError:
            x,y = np.loadtxt(filenames[1], unpack=True)
            x2, y2, l2 = mask_and_save(x,y,obses[1],spec,z,lines,show,save,savedir)
            
        # tmp = np.loadtxt(filenames[0])
        # x1,y1=tmp[:,0],tmp[:,1]#need to do that for the recursion for cases when len(filenames)>2
        # x2,y2 = np.loadtxt(filenames[1]        tmp = np.loadtxt(filenames[0])
        # x1,y1=tmp[:,0],tmp[:,1]#need to do that for the recursion for cases when len(filenames)>2
        # x2,y2 = np.loadtxt(filenames[1], unpack=True)
        w1,w2=[1,1]#weights
        #give up on 0s, too hard to keep track and not useful
        if len(y1)==len(y2):
            m=(y1!=0)&(y2!=0)
            x1=x1[m]
            y1=y1[m]
            x2=x2[m]
            y2=y2[m]
        else:
            x1=x1[y1!=0]
            y1=y1[y1!=0]
            x2=x2[y2!=0]
            y2=y2[y2!=0]
        
        xf,yf, coeff = average(x1,y1,x2,y2,w1,w2)
        if debug:
            debugdir = './merge_debug'
            if not os.path.exists(debugdir):
                os.mkdir(debugdir)
            fig = plt.figure(figsize=(20,9))
            plt.plot(x1,y1,label=os.path.basename(filenames[0])[:-4])
            plt.plot(x2,y2,label=os.path.basename(filenames[1])[:-4])
            plt.plot(xf,yf,label='merged')
            plt.title('%s'%str(spec))
            plt.legend()
            plt.savefig(os.path.join(debugdir,'MERGED%s.png'%str(spec)))
            plt.close(fig)

        mask_and_save(xf,yf,obs,spec,z,lines,show,save,savedir)
        
    # x_f = None
    # y_f = None
    # tmp = None
    # norm = 0
    # if len(filenames)>2:
    #     print obs_slits
    #     length=0
    #     for filename in filenames:
    #         x,y = np.loadtxt(filename, unpack=True)
    #         if length==0:
    #             length = len(x)
    #         elif len(x)!=length:
    #             print filenames
    #             return
    # for filename, obs_slit in zip(filenames, obs_slits):
    #     obs,slit=obs_slit.split(':')
    #     weight=obs_weight[obs]
    #     x,y = np.loadtxt(filename, unpack=True)
    #     if x_f is None:
    #         x_f = x
    #         y_f = weight*y
    #         tmp = y
    #         norm = weight
    #     else:
    #         c = (tmp*y).sum()/(y*y).sum()
    #         y_f += c*y
    #         norm += weight
            
    # mask_and_save(x_f,y_f,obs,spec,z,lines,show,save,savedir)
    # return
    
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

def check_specifics(obs_slits,spec, datadir, show, save, savedir):
    # if '255720-25fits:39' in obs_slits:
    #     print '255720-25fits:39 is ugly compared to 255730-33:38 : using only the latter'
    #     obs,slit='255730-33:38'.split(':')
    #     x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    # elif '255720-25fits:28' in obs_slits:
    #     print '255720-25fits:28 is ugly compared to 255730-33:29 : using only the latter'
    #     obs,slit='255730-33:29'.split(':')
    #     x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    # elif '255720-25fits:25' in obs_slits:
    #     print '255720-25fits:25 is ugly compared to 255730-33:26 : using only the latter'
    #     obs,slit='255730-33:26'.split(':')
    #     x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    # elif '255720-25fits:4' in obs_slits:
    #     print '255720-25fits:4 is ugly compared to 255730-33:6 : using only the latter'
    #     obs,slit='255730-33:6'.split(':')
    #     x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    # elif '255720-25fits:3' in obs_slits:
    #     print '255720-25fits:3 is ugly compared to 255730-33:5 : using only the latter'
    #     obs,slit='255730-33:5'.split(':')
    #     x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    # elif '255720-25fits:2' in obs_slits:
    #     print '255720-25fits:2 is ugly compared to 255730-33:4 : using only the latter'
    #     obs,slit='255730-33:4'.split(':')
    #     x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    if '255728:21' in obs_slits:
        #confirmed
        print '255728:21 is ugly compared to 255761:23 : using only the latter'
        obs,slit='255761:23'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255728:8' in obs_slits:
        #confirmed
        print '255728:8 is ugly compared to 255761:9 : using only the latter'
        obs,slit='255761:9'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255728:24' in obs_slits:
        print '255728:24 et 255763:33 are very problematic, keeping 255728:24 for now'
        obs,slit='255728:24'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    # elif '255720-25fits:32' in obs_slits:
    #     print '255720-25fits:32 is ugly compared to 255738:32 : using only the latter'
    #     obs,slit='255738:32'.split(':')
    #     x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    # elif '255734fits:47' in obs_slits:
    #     print '255734fits:47 is ugly compared to 255720-25fits:46 : using only the latter'
    #     obs,slit='255720-25fits:46'.split(':')
    #     x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255726:6' in obs_slits:
        print '255726:6 is ugly compared to 255752-59-60:22 : using only the latter'
        obs,slit='255752-59-60:22'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255726:5' in obs_slits:
        print '255726:5 is ugly compared to 255734fits:16 : using only the latter'
        obs,slit='255734fits:16'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255744_46_47:34' in obs_slits:
        print '255744_46_47:34 is ugly compared to 255720-25fits:43 : using only the latter'
        obs,slit='255720-25fits:43'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255739-43:36' in obs_slits:
        print '255739-43:36 is ugly compared to 255720-25fits:37 : using only the latter'
        obs,slit='255720-25fits:37'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    # elif '255739-43:30a' in obs_slits:
    #     print '255739-43:30a is ugly compared to 255720-25fits:45 : using only the latter'
    #     obs,slit='255720-25fits:45'.split(':')
    #     x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255726:7' in obs_slits:
        #confirmed
        print '255726:7 is ugly compared to 255750:37 : using only the latter'
        obs,slit='255750:37'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255748:7m' in obs_slits:
        print '255750/SPEC6 unfound but used in 255748/SPEC7_6, renamed SPEC7m'
        obs,slit='255748:7m'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255748:9m' in obs_slits:
        print '255750/SPEC8 unfound but used in 255748/SPEC9_8, renamed SPEC9m'
        obs,slit='255748:9m'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255748:11m' in obs_slits:
        print '255750/SPEC10 unfound but used in 255748/SPEC11_10, renamed SPEC11m'
        obs_slits.remove('255750:10')
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255750:24m' in obs_slits:
        print '255748/SPEC26 unfound but used in 255750/SPEC24m'
        obs,slit='255750:24m'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255748:27m' in obs_slits:
        print '255750/SPEC25 unfound but used in 255748/SPEC27m'
        obs,slit='255748:27m'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255761:38am' in obs_slits:
        print '255763/SPEC36a unfound but used in 255761/SPEC38am'
        obs,slit='255761:38am'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255761:38bm' in obs_slits:
        print '255763/SPEC36b unfound but used in 255761/SPEC38bm'
        obs,slit='255761:38bm'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255744_46_47:13' in obs_slits:
        print '255739-43/SPEC12 unfound, using only 255744_46_47/SPEC13'
        obs,slit='255744_46_47:13'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255728:31' in obs_slits:
        print '255744_46_47:42 unfound, using only 255728:31'
        obs,slit='255728:31'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255754-58fits:5m' in obs_slits:
        print '255752-59-60:5 unfound but used in 255754-58fits:5m'
        obs,slit='255754-58fits:5m'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255754-58fits:16' in obs_slits:
        print '255752-59-60:17 unfound, using only 255754-58fits:16'
        obs,slit='255754-58fits:16'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255752-59-60:13' in obs_slits:
        print '255752-59-60:13 includes 255734fits:29 and 255738:27'
        obs,slit='255752-59-60:13'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255752-59-60:23' in obs_slits:
        print '255752-59-60:23 not found, merging 255754-58fits:24 and 255744_46_47:11'
        obs_slits.remove('255752-59-60:23')
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255752-59-60:31' in obs_slits:
        print '255752-59-60:31 not found, merging 255720-25fits:29 and 255739-43:15'
        obs_slits.remove('255752-59-60:31')
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255752-59-60:10' in obs_slits:
        print '255752-59-60:10 not found, and  255754-58fits:10m includes 255734fits:28, merging 255754-58fits, 255730-33:27, and 255738:26'
        obs_slits.remove('255752-59-60:10')
        obs_slits.remove('255734fits:28')
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255752-59-60:2' in obs_slits:
        print '255752-59-60:2 not found, merging 255763:1 and 255761:1'
        obs_slits.remove('255752-59-60:2')
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255752-59-60:3' in obs_slits:
        print '255752-59-60:3 not found, merging 255752-59-60:3, 255763:2 and 255761:2'
        obs_slits.remove('255752-59-60:3')
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255752-59-60:19' in obs_slits:
        print '255752-59-60:19 not found, merging 255726:12, 255763:17 and 255761:16'
        obs_slits.remove('255752-59-60:19')
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255752-59-60:11' in obs_slits:
        print '255752-59-60:11 not found, merging 255720-25fits:4 and 255730-33:6'
        obs_slits.remove('255752-59-60:11')
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255748:47' in obs_slits:
        print '255748:47 not found, merging 255754-58fits:37 and 255750:45'
        obs_slits.remove('255748:47')
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255726:3' in obs_slits:
        print '255726:3 discarded, merging 255750:3 and 255761:20N'
        obs_slits.remove('255726:3')
        #obs_slits = ['255750:3','255761:20N','255726:3']
        do_merge(obs_slits, spec, datadir, show, save, savedir)
    elif '255750:13' in obs_slits:
        print '255750:13 not found, using only 255748:14'
        obs,slit='255748:14'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255754-58fits:14' in obs_slits:
        obs_slits=['255750:11','255739-43:20','255744_46_47:21','255748:12m','255754-58fits:14']
    elif '255750:25' in obs_slits:
        print '255750:25 not found, using only 255748:27'
        obs,slit='255748:27'.split(':')
        x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show, save, savedir)
    elif '255750:38' in obs_slits:
        print '255750:38 and 255748:40 not found, SPEC id %s results missing'%str(spec)
    else:
        return False
    return True

def obs_cmp(obs_slit1, obs_slit2):
    print obs_slit1, obs_slit2
    obs1,slit1 = obs_slit1.split(':')
    obs2,slit2 = obs_slit2.split(':')
    if obs1=='255720-25fits' and obs2=='255730-33':
        return -1
    if obs1=='255750':
        return 1
    elif obs2=='255750':
        return -1
    if obs1=='255748':
        return 1
    elif obs2=='255748':
        return -1
    if obs1=='255763':
        if obs2 in ['255750','255748']:
            return -1
        else:
            return 1
    if obs1=='255754-58fits':
        if obs2 in ['255750','255748','255763']:
            return -1
        else:
            return 1
    
    if obs1=='255728' and obs2=='255763':
        return -1
    if obs2=='255728' and obs1=='255763':
        return 1
    if obs1=='255754-58fits' or obs1=='255728':
        return -1
    return -1
        
    
def do_merge(obs_slits, spec, datadir, show, save, savedir):
    obs_weight = {"255706":1, "255708":1, "255712":1, "255714":1, "255716":1, "255718":1, "255720-25fits":5, '255726':1,'255728':1,"255730-33":3, "255734fits":1, '255744_46_47':3,"255754-58fits":1, "255750":1, "255763":1, "255748":1, "255761":1, "255738":1, "255739-43":4, "255752-59-60":3, "255744-45-47":3}
    if len(obs_slits)>=2:
        #in case of more than 2 spectra to merge, order may help,
        #so try here to put some order in the chaos.
        print "list before compare: ",obs_slits
        obs_slits = sorted(obs_slits, cmp=obs_cmp, reverse=True)
        print "list after compare: ",obs_slits
    filenames = []
    weights = []
    obses = []
    for obs_slit in obs_slits:
        obs, slit = obs_slit.split(":")
        filename = os.path.join(datadir,obs,'SPEC%sn.txt'%slit)
        if not os.path.exists(filename):
            filename = os.path.join(datadir,obs,'SPEC%s.txt'%slit)
            if not os.path.exists(filename):
                print "%s not found"%filename
            else:
                tmp_n=os.path.join(savedir,'to_be_merged',"%s_%s_%sn.txt"%(spec,obs,slit))
                if not os.path.exists(tmp_n):
                    do_normalize(filename, out=tmp_n)
                filenames.append(tmp_n)
                weights.append(obs_weight[obs])
                obses.append(obs)
        else:
            filenames.append(filename)
            weights.append(obs_weight[obs])
            obses.append(obs)
    print filenames
    merge_files(filenames, obses, weights,spec,z,lines,show, save, savedir)
    return
    
if __name__=="__main__":
    import yaml, sys
    config_file = 'config2.yaml'
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
        
    print("config: %s\nmapping: %s\nobservations: %s\n"%(config_file,map_table,obs_list))
        
    # for obs in obs_list:
    #     print("Running on ",obs)
    #     new_dir=os.path.join(workdir,obs)
    #     if not os.path.exists(new_dir):
    #         os.mkdir(new_dir)
    #     run_all(obs, datadir, new_dir, show=True, save=True)
    

#    for spec in range(1,739):
#    for spec in [272,273,274,290,293,300]:
#    for spec in [103]:
    for spec in range(1,10):
        outstr = "SPEC %s: "%str(spec)
        z, lines = get_catalog_info(spec, cat)
        obs_slits=[]
        for obs in obs_list:
            idx = np.where(np.array(obs_dict[obs]['spec'],dtype=str)==str(spec))[0]
            if len(idx)>0:
                idx=idx[0]
                slit = obs_dict[obs]['slit'][idx]
                obs_slits.append("%s:%s"%(str(obs),str(slit)))
        
        if len(obs_slits)==0:
            outstr += 'no match'
            print outstr
        elif len(obs_slits)==1:
            outstr += obs_slits[0]
            print outstr
            obs,slit=obs_slits[0].split(":")
            x, masked_y, z, lines_list = run_on_spec(obs, slit, spec, datadir, show=False, save=True, savedir=workdir)
        else:
                already_merged = check_specifics(obs_slits,spec, datadir, show=True, save=True, savedir=workdir)
                if not already_merged:
                        do_merge(obs_slits, spec, datadir, show=True, save=True, savedir=workdir)
        
