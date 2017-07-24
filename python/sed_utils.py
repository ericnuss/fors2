from matplotlib import pylab as plt
import lineid_plot
import numpy as np

a_lines={'FeII24_1':2344.2,'FeII24_2':2374.5,'FeII24_3':2382.8,'FeII26_1':2586.5,'FeII26_2':2600.2,'MgII_1':2795,'MgII_2':2802,'Ca(H)':3933.7, 'Ca(K)':3968.5, 'G':4304.4, 'Mg':5175.3, 'Na':5894.0, '4000break':4000.}
e_lines={'[OII]':3727.3,'H9':3835,'H8':3889,'Heps':3970.07,'Hdelta':4101.73,'Hgamma':4340.46,'Hbeta':4861.32,'[OIII]':4959.0,'[OIII]':5006.8, 'Halpha':6562.80,'S2':6716.0, 'NeIII':3869, 'NeIII_2':3968}


def clean(f, sigma=5, out=None):
    x,y=np.loadtxt(f,unpack=True)
    v=(y-y.mean())**2/len(y)
    mask = (y!=0)&(v<=sigma*v.mean())
    if out is not None:
        np.savetxt(out, np.c_[x[mask],y[mask]], fmt="%.6f")
        add_comment(out,["#mask: %s"%mask.tostring()])
    return x[mask], y[mask], mask

def has_normalized(flist,slit):
    res=[]
    for f in flist:
        suf=f.split(str(slit))[1][:-4]
        if 'n' in suf:
            res.append(f)
    return res

def is_normalized(f):
    return ('n' in f)

def do_normalize(f,out=None):
    hdlt=glob.glob('HDLT*.txt')
    if len(hdlt)>1:
        raise Exception('Several HDLT files found!')
    hdlt=hdlt[0]
    x,y=np.loadtxt(f,unpack=True)
    xx,ccd=np.loadtxt(hdlt, unpack=True)
    #align wavelength  of HDLT and spectrum
    if len(x)>len(xx):
        maxid = np.where(x==xx[-1])[0][0]
        x=x[:maxid+1]
        y=y[:maxid+1]
    #remove zeron entries of the ccd
    mask=ccd>0
    xx=xx[mask]
    ccd=ccd[mask]
    x=x[mask]
    y=y[mask]
    if not np.allclose(xx,x,1.e-10):
        raise Exception('HDLT.txt has a different wavelength array as the spectrum file %s'%f)
    slit=f.split('SPEC')[1][:2]
    if out is not None:
        np.savetxt(out, np.c_[x,y/ccd], fmt="%.6f")
    return x, y/ccd

def rebin(f,z,out=None):
    x,y=np.loadtxt(f,unpack=True)
    x/=(1.+z)
    if out is not None:
        np.savetxt(out, np.c_[x,y], fmt="%.6f")
    return x,y

def add_comment(f, comments):
    ff=open(f, 'r')
    lines=ff.readlines()
    ff.close()
    for count, comment  in enumerate(comments):
        lines.insert(count,comment+'\n')
    ff=open(f,'w')
    ff.writelines(lines)
    ff.close()

def plot_spectrum(x,y,z=0,lines=[],xlim=None, ylim=None,**kwargs):
    #fig = plt.figure(1)
    ax1=plt.subplot(111)
    #ax = fig.add_axes([0.1,0.06, 0.85, 0.35])
    # ax.plot(x, y)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    lvals=[]
    for l in lines:
        if l in a_lines:
            lvals.append(a_lines[l])
        elif l in e_lines:
            lvals.append(e_lines[l])
    #print lines, lvals
    lineid_plot.plot_line_ids(x/(1.+z), y*(1.+z), lvals, lines, ax=ax1)
    ax1.plot(x/(1.+z),y*(1.+z), **kwargs)
    if "label" in kwargs:
        ax1.legend()
    #plt.show()
    return ax1

if __name__=='__main__':
    import numpy as np
    import os
    from sed_utils import plot_spectrum
    #lines=['FeII24_1','FeII24_2','FeII24_3','FeII26_1','FeII26_2']
    #f='SPEC43n.txt'
    #zspec=1.303
    lines=['[OII]','MgII_1','MgII_2', 'Heps', 'Hbalmer1', 'Hbalmer2', 'Ca(H)','Ca(K)','NeIII']
    f='SPEC38n.txt'
    zspec=0.7633

    f='SPEC31.txt'
    zspec=0.4665
    x,y=np.loadtxt(f,unpack=True)

    ax=plot_spectrum(x,y, z=zspec, lines=lines, ylim=(-0.05,0.1))
    ax.set(title=f+':'+str(zspec))

    
