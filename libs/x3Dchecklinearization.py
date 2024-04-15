import sys
import matplotlib.pyplot as plt
import numpy as np
import polytope as pc
from .g_space import g, F, hsurf_g, hsurf_F, hsurf_F2
from .x3Dplot import plot_polytope


def getpoly_mitd( l, normal, distance, scom, dlist, imax=1/6):
    
    polylist = []
    
    gpsc  = np.identity(len(normal))
    Apsc  = np.array(np.vstack([-gpsc, gpsc]))
    bpsc  = np.array([0]*len(normal) + [imax]*len(normal))
    psc   = pc.Polytope(Apsc, bpsc)
    
    aa    = np.array(normal)
    bb    = np.array(distance)
    
    for d in dlist:
        d  = np.array(d)
        oo = np.cos(2*np.pi*l*d)
        if (np.all(np.sign(oo) == 1) or np.all(np.sign(oo) == -1)):
            for i in scom:
                
                A = []
                A.append(-i*aa)
                A.append( i*aa)
                
                if i[len(normal)-1]>0:
                    b=np.array(np.array([-i[len(normal)-1], i[len(normal)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                else:
                    b=np.array(np.array([i[len(normal)-1], -i[len(normal)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                # ---> inner
                iden = np.identity(len(normal))
                for k in range(len(normal)):
                    A=np.vstack([A,-1*iden[k]])
                
                de = d + (i-1)*(1/(4*l))
                b=np.append(b, -de)
                
                # ---> outter
                for k in range(len(normal)):
                    A=np.vstack([A,iden[k]])
                    
                de = d + 1*(i+1)*(1/(4*l))
                b=np.append(b, de)
                
                w=pc.Polytope(np.array(A),np.array(b))
                
                if w.chebXc is not None:
                    if (w.chebXc in psc):
                        polylist.append(w)
                
    return pc.Region(polylist)

def checklinear(l, xcoor, f, normal, distance, j=2, n=100, s=1, testiso=True):
    
    lspace  = np.linspace(0, 1/(2*l), n)
    kj = [lspace]*(len(f)-1)
    kz = np.meshgrid(*kj)
    #gz = np.zeros_like(kz[0])
    
    gi = np.abs(g(l, xcoor, f))
    
    gzp = hsurf_g(l, [*kz], f, gi, j, s=1)
    
    o = getpoly_mitd(l, normal, distance, scom=np.array([[1]*len(f)]), dlist=np.array([[0]*len(f)]), imax=lspace.max())
    
    if testiso:
        kz.extend([np.array(gzp)])
        tz = np.vstack(np.dstack([*kz]))
        
        for ti in tz:
            if np.all(~np.isnan(ti)):
                if ti in o[0]:
                    continue
                else:
                    print("\x1b[1;31m--> Checking the quality of linearization process ",end=" ")
                    print("\n--> Found isosurface outside for the point on the location ",ti)
                    print("--> The point is: ", ti,"isoutside the polytope \n")
                    print("\x1b[1;31m--> Check the linearization step <-- \x1b[0m")
                    print("\x1b[1;32m--> I am quitting hier, BYE <-- \x1b[0m")
                    sys.exit()
        print("\x1b[1;32m--> Polytope contains complete isosurface. Successful Linearization for \x1b[1;31mRO = %g\x1b[0m"%(l))
    
    return 

def checklinearplot(l, xexp, f, normal, distance, j=2,n=100, s=1, testiso=True, plot=True):
    
    lspace  = np.linspace(0, 1/(2*l), n)
    kj = [lspace]*(len(f)-1)
    kz = np.meshgrid(*kj)
    gz = np.zeros_like(kz[0])
    
    gi = np.abs(g(l, xexp, f))
    
    gzp = hsurf_g(l, [*kz], f, gi, j, s=1)
    
    o = getpoly_mitd(l, normal, distance, scom=np.array([[1]*len(f)]), dlist=np.array([[0]*len(f)]), imax=lspace.max())
    
    if testiso:
        kz.extend([np.array(gzp)])
        tz = np.vstack(np.dstack([*kz]))
        
        for ti in tz:
            if np.all(~np.isnan(ti)):
                if ti in o[0]:
                    continue
                else:
                    print("\x1b[1;31m--> Checking the quality of linearization process ",end=" ")
                    print("\n--> Found isosurface outside for the point on the location ",ti)
                    print("--> The point is: ", ti,"isoutside the polytope \n")
                    print("\x1b[1;31m--> Check the linearization step <-- \x1b[0m")
                    print("\x1b[1;32m--> I am quitting hier, BYE <-- \x1b[0m")
                    sys.exit()
        print("\x1b[1;32m--> Polytope contains complete isosurface. successful Linearization :) \n\x1b[0m")
    if plot:
        
        print("\x1b[1;31m---> plotting isosurface with polytope ...\x1b[0m")
        
        fig   = plt.figure(figsize = (6,5),frameon=False)
        ax    = plt.axes(projection='3d')
        
        fig.tight_layout()
        ax.set_xlabel(r'$x_\mathrm{1}$', fontsize=12)
        ax.set_ylabel(r'$x_\mathrm{2}$', fontsize=12)
        ax.set_zlabel(r'$x_\mathrm{3}$', fontsize=12)
        ax.grid(False)
        
        ax.plot_surface(kz[0], kz[1], gzp, color='k', alpha=0.5, antialiased=True,facecolor='r', linewidth=0)
        v=plot_polytope(o[0], ax, alpha=0.15, color ='C0')


def checkisosurf(l, I, f, normal, distance, n=100, s=1):
    
    j = len(f)-1
    lspace  = np.linspace(0, 1/(2*l), n)
    kj = [lspace]*(len(f)-1)
    kz = np.meshgrid(*kj)
    #gz = np.zeros_like(kz[0])
    
    gzp = hsurf_F2(I, l, [*kz], f, j, s=1, s2=1)
    o = getpoly_mitd(l, normal, distance, scom=np.array([[1]*len(f)]), dlist=np.array([[0]*len(f)]), imax=lspace.max())
    
    kz.extend([np.array(gzp)])
    tz = np.vstack(np.dstack([*kz]))
    x=tz[~np.isnan(tz).any(axis=1)]
    
    check=[i in o for i in x]
    #check=[ ti in o    for ti in tz   if ~np.all(np.isnan(ti)) ]
    
    if ~np.all(check):
            index = np.where(~np.array(check))[0]
            dr  = [np.dot(normal,x[inx]) for inx in index]
            return False, [np.min(dr), np.max(dr)]
    else:
        return [True]        
    return  
