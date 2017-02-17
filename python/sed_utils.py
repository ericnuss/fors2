from matplotlib import pylab as plt
import lineid_plot

a_lines={'FeII24_1':2344.2,'FeII24_2':2374.5,'FeII24_3':2382.8,'FeII26_1':2586.5,'FeII26_2':2600.2,'MgII_1':2795,'MgII_2':2802,'Ca(H)':3933.7, 'Ca(K)':3968.5, 'Gband':4304.4, 'Mg':5175.3, 'Na':5894.0, '4000break':4000.}
e_lines={'[OII]':3727.3,'Hbalmer1':3835,'Hbalmer2':3889,'Heps':3969.5,'Hdelta':4102.8,'Hgamma':4340.0,'Hbeta':4861.3,'[OIII]':4959.0,'[OIII]':5006.8, 'Halpha':6562.8,'S2':6716.0, 'NeIII':3869, 'NeIII_2':3968}


def plot_spectrum(x,y,z=0,lines=[],xlim=None, ylim=None):
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
    print lines, lvals
    lineid_plot.plot_line_ids(x/(1.+z), y, lvals, lines, ax=ax1)
    plt.plot(x/(1.+z),y)
    plt.show()
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

    
