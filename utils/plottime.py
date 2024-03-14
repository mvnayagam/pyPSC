import os
import h5py
import numpy as np

from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

from psc.x3Dreadwrite import readh5pymean


def plottime_nD(fpath, h, interval=1):
    
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)

    fig, axs = plt.subplots(1,1,figsize=(7,5.5),frameon=False,sharex=True,sharey=True,
                            constrained_layout=True,subplot_kw={'aspect':'auto'})

    boxpro     = dict(linestyle='--', linewidth=0.1, alpha=1)
    medianpro  = dict(color='k', linestyle='-', linewidth=1.5, alpha=1)
    meanpro    = dict(linestyle='-', linewidth=0.5, color='k')
    whiskerpro = dict(color="red", linewidth = 1.2, alpha=0.7)
    capprop    = dict(color="red", linewidth = 1.2, alpha=0.7) 


    hs = np.arange(2, h+1, interval) #[2, 4, 6, 8]
    mmean = []
    for ci, ls in enumerate(hs):
        fn = os.path.join(fpath,'pnew_%g.h5'%(ls))
        mmean.append(readh5pymean(fn)[1])
            
    cmap = plt.get_cmap("cool")
    norm = cm.colors.LogNorm(vmax=1*np.max(mmean), vmin=1*np.min(mmean))

    emax, t, nsol = [], [], []
    tave, eave, vave = [], [], []

    for ci, i in enumerate(hs):            
        fn  = os.path.join(fpath,'pnew_%g.h5'%(i))
        
        with h5py.File(fn, 'r') as f:
            
            tinf = f.get('timeexe')
            tt   = [np.array(tinf.get(ii)) for ii in np.array(tinf)] ; tt = np.array(tt)
            tti  = [np.max(np.array(tinf.get(ii)))for ii in np.array(tinf)]
            t.append(tt[:,6]); tave.append(tti)
            
            er    = f.get('error')
            err   = [np.array(er.get(i)) for i in np.array(er)]
            emaxi = [np.max(i) for i in err]
            emax.append(emaxi)
            
            ns  = f.get('allsolution')
            nss = [(np.shape(np.array(f.get('allsolution/Pair'+str(ii))))[0])/2 for ii in range(1, len(ns)+1)]
            nsol.append(nss)
            
            vave.append(readh5pymean(fn)[1])

    ##### plot and change properties

    bp = axs.boxplot(t, labels=hs, notch=False, patch_artist=True, boxprops=boxpro, meanprops=meanpro,
                                medianprops=medianpro, meanline=False, showmeans=False, showfliers=False, 
                                whiskerprops=whiskerpro, capprops=capprop)

    for bi, box in enumerate(bp['boxes']):
        box.set(facecolor=cmap(norm(vave[bi])), linewidth=1, alpha=0.65) 
    for bi, box in enumerate(bp['whiskers']):
        if bi%2==0:
            binx = int(bi/2)
        box.set(color=cmap(norm(vave[binx])), linewidth=2, alpha=1)
    for bi, box in enumerate(bp['caps']):
        if bi%2==0:
            binx = int(bi/2)
        box.set(color=cmap(norm(vave[binx])), linewidth=2, alpha=1)


    ####### plot secondary axis

    csel = 'green' ; csels = 'orangered'

    ax_err=axs.twinx() 

    ax_err.set_ylim(-0.01, 0.12)
    ax_err.set_xticklabels([]);
    ax_err.yaxis.label.set_color(csel); ax_err.tick_params(axis='y', colors=csel)

    bp_err = ax_err.boxplot(emax, labels=hs, notch=False, showfliers=False, whis=False, patch_artist=False,
                    showcaps=False, showbox=False,medianprops=dict(linestyle='-', linewidth=0, alpha=1, color='k'),
                    showmeans=True, meanline=False, meanprops=dict(linestyle='-', linewidth=0, color=csel,
                    marker="o",markerfacecolor=csel, markeredgecolor="k", alpha=0.5))


    ax_sol=axs.twinx()

    ax_sol.spines.right.set_position(("axes", 1.4));
    ax_sol.set_ylabel("no. of solutions", fontsize=12, labelpad=10.0)
    ax_err.set_ylabel("max. error",       fontsize=12, labelpad=5.0)
    ax_sol.yaxis.label.set_color(csels) ; ax_sol.tick_params(axis='y', colors=csels)

    bp_sol = ax_sol.boxplot(nsol, labels=hs, notch=False, showfliers=False, whis=False, patch_artist=False,
                showcaps=False, showbox=False,medianprops=dict(linestyle='-', linewidth=0, alpha=1, color='k'),
                showmeans=True, meanline=False, 
                meanprops=dict(linestyle='-', linewidth=0, color=csels, marker="o",markerfacecolor=csels, 
                            markeredgecolor="k", alpha=0.5))
        
    X, Y = [], []
    for m in bp_err['means'][:]:
        lx, ly = m.get_data() ; X.append(lx[0]); Y.append(ly[0])
    ax_err.plot(X,Y,'-',c=csel, lw=1.5, alpha=0.5, zorder=-1)

    X, Y = [], []
    for m in bp_sol['means'][:]:
        lx, ly = m.get_data();  X.append(lx[0]);  Y.append(ly[0])
    ax_sol.plot(X, Y, '-', c=csels, lw=1.5, alpha=0.5, zorder=-1)

    axs.grid(which='major', color='#CCCCCC', lw=0.5, linestyle='--',alpha=0.5)

    axs.set_ylabel (r"time $(s)$",  fontsize=16,labelpad=12.0)
    axs.set_xlabel(r"Reflection order", fontsize=16,labelpad=12.0)
    
    sm   =  ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=axs, orientation='vertical', pad=0.1, aspect=20, shrink=0.95)
    cbar.set_label(r'$\log~ Volume_{\mathrm{~average}}$', size=14, weight='bold', loc='center', labelpad=0.15)
    cbar.ax.tick_params(labelsize=14)
    
    #plt.savefig(os.path.join(os.getcwd(), 'MCres/folder1/folder2/timeinfo.pdf'), dpi=200, bbox_inches='tight')  # when want to save in some other location
    plt.savefig(os.path.join(fpath, 'timeinfo.pdf'), dpi=250, bbox_inches='tight')   # when want to save in the some location of results directory
    
    plt.show()
    
    return
