import sys
import numpy as np
from .g_space import g, hsurf_F2, hsurf_g

from scipy.optimize import root


def checklinear(l, xcoor, f, normal, distance, j=2, n=200, s=1, testiso=True):
    
    lspace  = np.linspace(0, 1/(2*l), n)
    kj = [lspace]*(len(f)-1)
    kz = np.meshgrid(*kj)
    gz = np.zeros_like(kz[0])
    
    gi = g(l, xcoor, f)**2
    
    gzp= hsurf_F2(gi, l, [*kz], f, j=2, s=1, s2=1)   ## check hsurf_F2 argment position
    
    o  = getpoly_mitd(l, normal, distance, scom=np.array([[1]*len(f)]), dlist=np.array([[0]*len(f)]), imax=lspace.max())
    
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

def extremetangent(h, gi, f, normal, percentage=1):
    
    def func_vec(x0, h, gi, f, normal):
        k  = 2*np.pi*h
        gg = np.sum([ f[i]* np.sqrt(1- ((x0[i]*normal[i])/(k*f[i]))**2 ) for i in range(len(f))])
        return [gi-gg]*len(x0)
    
    def funcjac_vec(x0, h, f, n):
        k     = 2*np.pi*h
        gradg = [ f[i]*x0[i]*((n[i]/(k*f[i]))**2)*(1-((x0[i]*n[i])/(k*f[i]))**2)**(-1/2) for i in range(len(f)) ]
        return gradg
    
    k  = 2*np.pi*h  ; k1 = 1/k
    
    lam_max  = [ k*(f[i]/normal[i]) if (normal[i] >1E-5) else k*f[i]*normal[i] for i in range(len(f))]
    lam_mask = np.ma.masked_equal(lam_max, 0.0, copy=False)
    
    x0 = [(ik - ik%percentage)-1  if ik !=0 else normal[ci]*ik for ci, ik in enumerate(lam_max) ]
    
    res    = root(func_vec, x0, args=(h, gi, f, normal), options={'maxiter':100,'gtol': 1e-16, 'disp': True}).x
    
    midpnt = [ k1*np.arcsin((res[i]*normal[i])/(k*f[i])) for i in range(len(f))]
    
    return [res, midpnt]

def remaining(l, gi, f, n_select, di):
    
    mp =  np.abs(extremetangent(l, gi, f, n_select, percentage=5.)[1])
    do = np.abs(np.dot(n_select, mp))
    
    distance = [di, do]
    
    dok=checklinear(l, xexp, f, n_select, distance, j=2, n=500, s=1)
    
    if not dok[0]: # dok[1] = [di_n, do_n]
        do = dok[1][1] if ~dok[0] and do <dok[1][1] else do
        di = dok[1][0] if ~dok[0] and di >dok[1][0] else di
    
    distance  = [di, do]
    
    return distance

def plane_eq2n(pnt):
    
    p1 = np.array(pnt[0])
    p2 = np.array(pnt[1])
    p3 = np.array(pnt[2])
    
    # These two vectors are in the plane
    v1 = p3 - p2
    v2 = p1 - p2
    
    # the cross product is a vector normal to the plane
    cp = np.cross(v1, v2)
    cp = np.array(cp, dtype=float)
    #print("n before normalization : ", cp)
    cp /= np.linalg.norm(cp)
    
    a, b, c = cp
    
    a = a if ~np.isnan(a) else 0.5/1
    b = b if ~np.isnan(b) else 0.5/1
    c = c if ~np.isnan(c) else 0.5/1
    
    n = np.array([a, b, c], dtype = float)
    d =  np.abs(np.dot(n, p1))
    
    dist1 = np.abs(d)  #/ np.linalg.norm(n)
    
    return np.abs(np.array([a, b, c])), d

def getxyz_opt3(l, f, I, minmax=0):
    
    xyz = []
    
    for jj in [0, 0.5/l]:
        for i in range(len(f)-1):
            t = [0] * len(f)
            t[len(f)-1] = jj
            
            z = hsurf_F2(I*I, l, t, f, j=i, s=1, s2=1) # hsurf_g(l, t, f, I, j=i, s=1)
            
            if not np.isnan(z):
                t[i] = z
            else:
                t[i] = 0.5/l
                z = hsurf_F2(I*I, l, t, f, j=i, s=1, s2=1) # hsurf_g(l, t, f, I, j=i-1, s=1)
                if not np.isnan(z):
                    t[i-1] = z
                else:
                    t[i-1] = 0.5/l
            xyz.append(t)
    
    a = np.array(xyz)
    
    if minmax==0:
        a[:,0] = np.min(a[:,0], axis=0)
    else:
        a[:,0] = np.unique(a)[-2] if np.unique(a)[-1] == 0.5/l else  np.unique(a)[-1]
        
    return a, np.array(xyz)

def getxyz(l, f, gi):
    xyz = []
    
    for i in range(len(f)-1):
        for jj in [0, 1/(2*l)]:
            t = [0]*len(f)
            t[len(f)-1] = jj
            #z = hsurf_g(l, t, f, gi, j=i, s=1)
            z = hsurf_F2(gi, l, t, f, j=i, s=1, s2=1)
            
            if ~np.isnan(z):
                t[i] = z
            else:
                t[i] = 0.5/l
                #z = hsurf_g(l, t, f, gi, j=i, s=1)
                z = hsurf_F2(gi*gi, l, t, f, j=i, s=1, s2=1)
                
                if not np.isnan(z):
                    t[i-1] = z
                else:
                    t[i-1] = 0.5/l
                            
            xyz.append(t)
    
    a = np.copy(xyz)
    a = a[a[:, 1].argsort()][::-1]
    a[0][:2] = a[1][:2] if np.all(a[1][:2] <= a[0][:2]) else a[0][:2]
    a[1][:2] = a[0][:2] if np.all(a[1][:2] >= a[0][:2]) else a[1][:2]
    
    a[2][:2] = a[3][:2] if np.all(a[3][:2] <= a[2][:2]) else a[2][:2]
    a[3][:2] = a[2][:2] if np.all(a[3][:2] >= a[2][:2]) else a[3][:2]
    
    return np.array(a), np.array(xyz)    

def get_pnts(l, f, I):
    tot_coor=[]
    for i in range(len(f)):
        tem_coor    = np.zeros(len(f))
        tem_coor[i] = hsurf_F2(I*I, l, tem_coor, f, j=i, s=1, s2=1) #hsurf_g(l, tem_coor, f, gi, i, s=1)
        tot_coor.append(tem_coor)
        
    return tot_coor

def gt3(l, f, I):
    xyz = []
    
    t = [0] * len(f)
    z = hsurf_F2(I*I, l, t, f, j=0, s=1, s2=1)
    
    if not np.isnan(z):
        t[0] = z
    else:
        t[0] = 0.5/l
    ###
    xyz.append(t)
    
    for i in range(1,len(f)):
        t = [0] * len(f)
        t[i] = 0.5/l
        z = hsurf_F2(I*I, l, t, f, j=0, s=1, s2=1)
        
        if not np.isnan(z):
            t[0] = z
        else:
            t[0] = 0.5/l
        
        xyz.append(t)
    
    t = [0.5/l] * len(f)
    z = hsurf_F2(I*I, l, t, f, j=0, s=1, s2=1)
    
    if not np.isnan(z):
        t[0] = z
    else:
        t[0] = 0.5/l
    
    xyz.append(t)
    
    a = np.array(xyz) ; a[:,0] = np.unique(a)[-1] if np.unique(a)[-1] == 0.5/l else  np.unique(a)[-1]
    return xyz, a

def newfnq(l, f, gi):
    
    pnt = get_pnts(l, f, gi)
    count_Nan = np.count_nonzero(np.isnan(pnt))
    
    
    if   count_Nan == 0:
        f = np.array(f)
        k   = 2*np.pi*l   ;     xp  = (1/k)*np.arccos(gi/np.sum(f))
        p4  = [xp]*len(f)
        
        pnt = np.vstack([np.array(pnt),p4])
        
        ni, di   = plane_eq2n(pnt)
        
        n_select     = ni
        distance_all = remaining(l, gi, f, n_select, di)
        
    elif count_Nan == 1:
        pi, pext    = getxyz(l, f, gi) # getxyz_opt3(l, f, gi, minmax=min) # 
        
        n_ext, dext = plane_eq2n(pext)
        ni, di      = plane_eq2n(pi)
        n_select    = ni #n_ext
        distance_all = remaining(l, gi, f, n_select, di)
        
    elif count_Nan > 1:
        pi, pext    = gt3(l, f, gi)
        
        n_ext, dext = plane_eq2n(pext)
        ni, di      = plane_eq2n(pi)
        n_select    = ni # n_ext
        
        di = di if di !=0 and di<dext else dext       
        distance_all = remaining(l, gi, f, n_select, di)
        
    else:
        print("--> count_Nan ", count_Nan)
        print("--> Wired isosurface check linearization routine ")
        print("--> predicted points are : \n", pnt)
    
    #print("n_select, distance_all : ", n_select, distance_all)
    
    return n_select, distance_all
