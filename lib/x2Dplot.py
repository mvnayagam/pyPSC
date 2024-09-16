
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import cm

from psc.g_space import g, hsurf_g, F, hsurf_F, hsurf_F2


def plot_segment(ax, p, cc, lww=1.5, al=0.5):
    for i in p:
        ax.plot(i[:,0], i[:,1],'-',lw=lww, c=cc, alpha=al)
    return

def plotisosurf_EPA(l, h, gi, ax, isos, giso1, giso2, cc, lw=0.12, imax=0.5):
    
    for hi in range(l+1):
        if ( (hi/l <= imax and h%2 !=0) or (hi/l <= imax and h%2 ==0) ):
            if hi == 0:
                ax.plot(isos, giso1 + hi/l, '-' ,lw=lw, c=cc,label=r'$\mathcal{G }\mathrm{( %g,%1.2f)}$'%(l, gi))
                ax.plot(isos, giso2 + hi/l, '--',lw=lw, c=cc)
            else:
                if (hi/l < imax and l%(2*hi) !=0):
                    ax.plot(isos, giso1 + hi/l, '-', lw=lw, c=cc)
                    ax.plot(isos, giso2 + hi/l, '--',lw=lw, c=cc)
                ax.plot(isos, -1*giso1  + hi/l, '-', lw=lw, c=cc)
                ax.plot(isos, -1*giso2  + hi/l, '--',lw=lw, c=cc)
                if (hi/l < imax and l%(2*hi) ==0):
                    ax.plot(isos, giso1    + hi/l, '-', lw=lw, c=cc)
                    ax.plot(isos, giso2    + hi/l, '--',lw=lw, c=cc)
                
            if (l%(2*l) == 0 and l/h <= imax):
                ax.plot(isos, -1*giso1 + (hi+2)/l, '-', lw=lw,c=cc)
                ax.plot(isos, -1*giso2 + (hi+2)/l, '--',lw=lw,c=cc)
    return

def plotisosurf_nEPA(l, h, gi, ax, isos, y1, y2, y3, y4, cc, lw=2, imax=0.5, alp=0.5):
    for hi in range(l+1):
        if ( (hi/l <= imax and h%2 !=0) or (hi/l <= imax and h%2 ==0) ):
            if hi == 0:
                ax.plot(isos, y1 + hi/(l), '-',  c='k', alpha=alp, label='h=%g'%(l))
                ax.plot(isos, y3 + hi/(l), '--', c='b', alpha=alp)
            else:
                if (hi/l <imax and l%(2*hi) !=0):
                    ax.plot(isos, y1 + hi/(l), '-',  c='k', alpha=alp)
                    ax.plot(isos, y3 + hi/(l), '--', c='b', alpha=alp)
                    
                ax.plot(isos, y2 + hi/(l), '-',  c='r', alpha=alp)
                ax.plot(isos, y4 + hi/(l), '--', c='g', alpha=alp)
                
                if (hi/l < imax and l%(2*hi) ==0):
                    ax.plot(isos, y1 + hi/(l), '-',  c='k', alpha=alp)
                    ax.plot(isos, y3 + hi/(l), '--', c='b', alpha=alp)
                
                if (l%(2*l) == 0 and l/h <= imax):
                    ax.plot(isos, y1 + (hi+2)/l, '-',  c='k', alpha=alp)
                    ax.plot(isos, y2 + (hi+2)/l, '-',  c='r', alpha=alp)
    return
