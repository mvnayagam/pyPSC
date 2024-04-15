import numpy as np
import polytope as pc
from psc.g_space import g, F, hsurf_g, hsurf_F, hsurf_F2
from psc.x3Dchecklinearization import getpoly_mitd

from scipy.optimize import minimize, minimize_scalar, root, basinhopping

# ============== Modules for EPA ==============

def get_pnts_new4(l, f, gi, xinit=0, imax=0.5):
    
    k = 2*np.pi*l  ; tot_coor = []
    
    for i in range(len(f)):        
        tem_coor    = np.zeros(len(f))
        tem_coor[i] = hsurf_g(l, tem_coor, f, gi, i, s=1)
        
        tot_coor.append(tem_coor)        
    
    if np.all(~np.isnan(np.array(tot_coor))):
        return np.array(tot_coor)
    else:
        
        inx = np.argwhere(np.isnan(tot_coor).any(axis=1)).flatten()
        
        for iw in inx:
            tem_coor    = np.zeros(len(f))
            tem_coor[iw] = xinit/l
            
            mlist       = list([ii for ii in range(len(f)) if ii != iw])
            
            temp_coor = (1/k)*np.arccos( (gi-f[iw]*np.cos(k*tem_coor[iw])) / 
                                         (np.sum([f[ii]*np.cos(k*tem_coor[ii]) for ii in mlist])))
            
            tem_coor[iw+1:] = temp_coor
            tem_coor[:iw]   = temp_coor
            
            tot_coor[iw]=tem_coor
            
    return np.array(tot_coor)

def get_pnts_new(l, f, gi, imax=0.5):
    tot_coor=[]
    for i in range(len(f)):
        tem_coor    = np.zeros(len(f))
        #tem_coor[i] = x[i]
        tem_coor[i] = hsurf_g(l, tem_coor, f, gi, i, s=1)
        
        
        if (~np.isnan(tem_coor[i])):
            tem_coor[i] = tem_coor[i]
        else:
            tem_coor[i] = imax/l
            tem_coor[i] = hsurf_g(l, tem_coor, f, gi, i, s=1)
            
        tot_coor.append(tem_coor)
    
    return tot_coor

def linearizenD_EPA(l:int, f: list, gi: int) ->list:
    
    k = 2*np.pi*l
    #======= 1. finding first three points
    pnt = get_pnts_new(l, f, gi)
    
    #======= 2. finding fourth point : the point on the surface :: At particular point x=y=z or x1*=x2*=x3*. so
    xp  = (1/k)*np.arccos(gi/np.sum(f))
    p4  =[xp]*len(f)
    
    if np.all(~np.isnan(pnt)): #~np.isnan(a):
        centroid = np.mean(pnt, axis=0)
        u, s, v  = np.linalg.svd(pnt-centroid)
        #n = v[-1] if np.all((pnt-centroid)[0]>=0) else -1*v[-1]
        n = np.abs(v[-1])
        
        d_int = np.double(np.sum(np.multiply(n,centroid)))
        d_out = np.double(np.sum(np.multiply(n,p4)))
        
        aps = [d_out/i for i in n]
                        
        for iv in range(len(f)):
            vv = np.zeros(len(f))
            vv [iv] = aps[iv]
            pnt = np.vstack([pnt, vv])
        
        normal = n
        d_all  = [d_int,d_out]
    else:
        
        pnt=[]
        
        pin = get_pnts_new4(l, f, gi, xinit=0)
        pnt = np.vstack([pin])
               
        centroidi   = np.mean(pin, axis=0)
        ui, si, vi  = np.linalg.svd(pin-centroidi)
        
        d_int    = np.abs(np.sum(np.multiply(-1*vi[-1],centroidi)))
                
        pon = get_pnts_new4(l, f, gi, xinit=0.5)
        pnt = np.vstack([pnt, pon])
       
        centroido   = np.mean(pon, axis=0)
        uo, so, vo  = np.linalg.svd(pon-centroido)
        
        d_out1      = np.abs(np.sum(np.multiply(vo[-1],centroido)))
        d_out2      = np.abs(np.dot([xp]*len(f),vo[-1]))
        
        d_all  = [d_int, d_out1] if d_out1 > d_out2 else  [d_int, d_out2]
        
        normal = np.abs(vo[-1])
        
    return normal, d_all




# ============== Modules for non EPA ==============

def linearizenD_nEPA(l, f, gi):
    
    from scipy.optimize import minimize, minimize_scalar, root, basinhopping
    
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
        
        mp = np.abs(extremetangent(l, gi, f, n_select, percentage=5.)[1])
        do = np.abs(np.dot(n_select, mp))
        
        distance = [di, do]
        
        dok=checkisoa(l, gi, f, n_select, distance, n=100, s=1)
        
        if not dok[0]: # dok[1] = [di_n, do_n]
            do = dok[1][1] if ~dok[0] and do <dok[1][1] else do
            di = dok[1][0] if ~dok[0] and di >dok[1][0] else di
        
        distance  = [di, do]
        
        return distance
    
    def checkisoa(l, gi, f, normal, distance, n=100, s=1):
        j = len(f) - 1
        lspace = np.linspace(0, 1 / (2 * l), n)
        kj = [lspace] * (len(f) - 1)
        kz = np.meshgrid(*kj)
        
        gzp = hsurf_F2(gi, l, [*kz], f, j, s=1, s2=1)
        #print(f"gzp shape: {np.shape(gzp)}")
        
        o = getpoly_mitd(l, normal, distance, scom=np.array([[1] * len(f)]), dlist=np.array([[0] * len(f)]), imax=lspace.max())
        
        kz.append(np.array(gzp))
        
        t = [kzi.flatten() for kzi in kz]
        t = list(zip(*t))
        t2 = [list(ti) for ti in t]
        t2 = np.array(t2)
        tz = t2[~np.isnan(t2).any(axis=1)]
        
        #print(f"t2: {np.shape(t2)} tz: {np.shape(tz)}")

        check = [ti in o for ti in tz]

        if np.any(np.logical_not(np.array(check))):
            index = np.where(~np.array(check))[0]
            dr = [np.dot(normal, tz[inx]) for inx in index]
            return False, [np.min(dr), np.max(dr)]
        else:
            return [True]

    def checkiso(l, gi, f, normal, distance, n=100, s=1):
        
        j = len(f)-1
        lspace  = np.linspace(0, 1/(2*l), n)
        kj = [lspace]*(len(f)-1)
        kz = np.meshgrid(*kj)
        gz = np.zeros_like(kz[0])
        
        gzp = hsurf_F2(gi, l, [*kz], f, j, s=1, s2=1)
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
    
    def plane_eq2n(pnt):
        
        centroid = np.mean(pnt, axis=0)
        u, s, v  = np.linalg.svd(pnt-centroid)
        
        d_int = np.double(np.sum(np.multiply(np.abs(v[-1]),centroid)))
                
        return np.abs(v[-1]), d_int #np.abs(np.array([a, b, c])), d

    def getxyz_opt3(l, f, I, minmax=0):
        
        xyz = []
        
        for jj in [0, 0.5/l]:
            for i in range(len(f)-1):
                t = [0] * len(f)
                t[len(f)-1] = jj
                
                z = hsurf_F2(I*I, l, t, f, j=i, s=1, s2=1) # hsurf_g3d(l, t, f, I, j=i, s=1)
                
                if not np.isnan(z):
                    t[i] = z
                else:
                    t[i] = 0.5/l
                    z = hsurf_F2(I*I, l, t, f, j=i, s=1, s2=1) # hsurf_g3d(l, t, f, I, j=i-1, s=1)
                    if not np.isnan(z):
                        t[i-1] = z
                    else:
                        t[i-1] = 0.5/l
                xyz.append(t)
        
        a = np.array(xyz)
        
        if minmax==0:
            #a[:,0:1] = np.min(a[:,0:1])
            a[:,0] = np.min(a[:,0], axis=0)
        else:
            #a[:,0:1] = np.max(a[:,0:1])
            a[:,0] = np.unique(a)[-2] if np.unique(a)[-1] == 0.5/l else  np.unique(a)[-1]
            
        return a, np.array(xyz)
    
    def getxyz(l, f, gi):
        xyz = []
        
        for i in range(len(f)-1):
            for jj in [0, 1/(2*l)]:
                t = [0]*len(f)
                t[len(f)-1] = jj
                z = hsurf_g(l, t, f, gi, j=i, s=1)
                #z = hsurf_F2a(gi, l, t, f, j=i, s=1, s2=1)
                
                if ~np.isnan(z):
                    t[i] = z
                else:
                    t[i] = 0.5/l
                    #z = hsurf_g3d(l, t, f, gi, j=i, s=1)
                    z = hsurf_F2(gi*gi, l, t, f, j=i, s=1, s2=1)
                    
                    if not np.isnan(z):
                        t[i-1] = z
                    else:
                        t[i-1] = 0.5/l
                                
                xyz.append(t)
        
        a = np.copy(xyz)
        a = a[a[:, 1].argsort()][::-1]
        
        #print("\n a - before :: \n", a)
        a[0][:2] = a[1][:2] if np.all(a[1][:2] <= a[0][:2]) else a[0][:2]
        a[1][:2] = a[0][:2] if np.all(a[1][:2] >= a[0][:2]) else a[1][:2]
        
        a[2][:2] = a[3][:2] if np.all(a[3][:2] <= a[2][:2]) else a[2][:2]
        a[3][:2] = a[2][:2] if np.all(a[3][:2] >= a[2][:2]) else a[3][:2]
        
        #print("\n a - after :: \n", a)
        return np.array(a), np.array(xyz)    
    
    def get_pnt_nEPA(l, f, I):
        tot_coor=[]
        for i in range(len(f)):
            tem_coor    = np.zeros(len(f))
            tem_coor[i] = hsurf_F2(I*I, l, tem_coor, f, j=i, s=1, s2=1, nan=False) #hsurf_g3d(l, tem_coor, f, gi, i, s=1)
            tot_coor.append(tem_coor)
            
        return tot_coor
    
    def gt3(l, f, I):
        xyz = []
        
        t = [0] * len(f)
        z = hsurf_F2(I*I, l, t, f, j=0, s=1, s2=1, nan=False)
        
        if not np.isnan(z):
            t[0] = z
        else:
            t[0] = 0.5/l
        ###
        xyz.append(t)
        
        for i in range(1,len(f)):
            t = [0] * len(f)
            t[i] = 0.5/l
            z = hsurf_F2(I*I, l, t, f, j=0, s=1, s2=1, nan=False)
            
            if not np.isnan(z):
                t[0] = z
            else:
                t[0] = 0.5/l
            
            xyz.append(t)
        
        t = [0.5/l] * len(f)
        z = hsurf_F2(I*I, l, t, f, j=0, s=1, s2=1, nan=False)
        
        if not np.isnan(z):
            t[0] = z
        else:
            t[0] = 0.5/l
        ###
        xyz.append(t)
        
        a = np.array(xyz) ; a[:,0] = np.unique(a)[-1] if np.unique(a)[-1] == 0.5/l else  np.unique(a)[-1]
        return xyz, a
    
    
    # ---> Main func starts
    
    pnt = get_pnt_nEPA(l, f, gi)
    count_Nan = np.count_nonzero(np.isnan(pnt))
    
    if   count_Nan == 0:
        f = np.array(f)
        k   = 2*np.pi*l ; xp = (1/k)*np.arccos(gi/np.sum(f))
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
        #pi, pext  = getxyz_opt3(l, f, gi, minmax=1)
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
        
    return n_select, distance_all



# -------------- bak - not working properly
# def linearizenD_nEPA(l, f, gi):
    
#     def extremetangent(h, gi, f, normal, percentage=1):
        
#         def func_vec(x0, h, gi, f, normal):
#             k  = 2*np.pi*h
#             gg = np.sum([ f[i]* np.sqrt(1- ((x0[i]*normal[i])/(k*f[i]))**2 ) for i in range(len(f))])
#             return [gi-gg]*len(x0)
        
#         def funcjac_vec(x0, h, f, n):
#             k     = 2*np.pi*h
#             gradg = [ f[i]*x0[i]*((n[i]/(k*f[i]))**2)*(1-((x0[i]*n[i])/(k*f[i]))**2)**(-1/2) for i in range(len(f)) ]
#             return gradg
        
#         k  = 2*np.pi*h  ; k1 = 1/k
        
#         lam_max  = [ k*(f[i]/normal[i]) if (normal[i] >1E-5) else k*f[i]*normal[i] for i in range(len(f))]
#         lam_mask = np.ma.masked_equal(lam_max, 0.0, copy=False)
        
#         x0 = [(ik - ik%percentage)-1  if ik !=0 else normal[ci]*ik for ci, ik in enumerate(lam_max) ]
        
#         res    = root(func_vec, x0, args=(h, gi, f, normal), options={'maxiter':100,'gtol': 1e-16, 'disp': True}).x
        
#         midpnt = [ k1*np.arcsin((res[i]*normal[i])/(k*f[i])) for i in range(len(f))]
        
#         return [res, midpnt]
    
#     def remaining(l, gi, f, n_select, di):
        
#         mp = np.abs(extremetangent(l, gi, f, n_select, percentage=5.)[1])
#         do = np.abs(np.dot(n_select, mp))
        
#         distance = [di, do]
        
#         dok=checkiso(l, gi, f, n_select, distance, n=500, s=1)
        
#         if not dok[0]: # dok[1] = [di_n, do_n]
#             do = dok[1][1] if ~dok[0] and do <dok[1][1] else do
#             di = dok[1][0] if ~dok[0] and di >dok[1][0] else di
        
#         distance  = [di, do]
        
#         return distance
    
#     def cheiso(l, gi, f, normal, distance, n=100, s=1):
#         j = len(f)-1
#         lspace = np.linspace(0, 1/(2*l), n)
#         gx, gy = np.meshgrid(lspace, lspace)
        
#         gz     = np.zeros_like(gx)
        
#         #gii    = F(l, xexp, f)**2; 
#         gzp = hsurf_F2(gi, l, [gx, gy, gz], f, j, s=1, s2=1)
#         #gii    = np.abs(g3d(l, xexp, f));
#         #gzp = hsurf_g3d(l, [gx, gy, gz], f, gii, j, s=s)
        
#         scom=np.array([[1]*len(f)])  ;  dlist=np.array([[0]*len(f)])
        
#         o = getpoly_mitd(l, normal, distance, scom, dlist, imax=np.max(lspace))
        
#         points = np.column_stack((gx.ravel(), gy.ravel(), gzp.ravel()))
        
#         x=points[~np.isnan(points).any(axis=1)]
#         check=[i in o for i in x]
         
#         if ~np.all(check):
#             inx = np.where(~np.array(check))[0]
#             dr  = [np.dot(normal,x[i]) for i in inx]
            
#             return False, [np.min(dr), np.max(dr)]
        
#         else:
#             return [True]
#         return

#     def checkiso(l, gi, f, normal, distance, n=100, s=1):
        
#         j = len(f)-1
#         lspace  = np.linspace(0, 1/(2*l), n)
#         kj = [lspace]*(len(f)-1)
#         kz = np.meshgrid(*kj)
#         gz = np.zeros_like(kz[0])
        
#         gzp = hsurf_F2(gi, l, [*kz], f, j, s=1, s2=1)
#         o = getpoly_mitd(l, normal, distance, scom=np.array([[1]*len(f)]), dlist=np.array([[0]*len(f)]), imax=lspace.max())
        
#         kz.extend([np.array(gzp)])
#         tz = np.vstack(np.dstack([*kz]))
#         x=tz[~np.isnan(tz).any(axis=1)]
        
#         check=[i in o for i in x]
#         #check=[ ti in o    for ti in tz   if ~np.all(np.isnan(ti)) ]
        
#         if ~np.all(check):
#                 index = np.where(~np.array(check))[0]
#                 dr  = [np.dot(normal,x[inx]) for inx in index]
#                 return False, [np.min(dr), np.max(dr)]
#         else:
#             return [True]        
#         return  
    
#     def plane_eq2n(pnt):
        
#         centroid = np.mean(pnt, axis=0)
#         u, s, v  = np.linalg.svd(pnt-centroid)
        
#         d_int = np.double(np.sum(np.multiply(np.abs(v[-1]),centroid)))
                
#         return np.abs(v[-1]), d_int #np.abs(np.array([a, b, c])), d

#     def getxyz_opt3(l, f, I, minmax=0):
        
#         xyz = []
        
#         for jj in [0, 0.5/l]:
#             for i in range(len(f)-1):
#                 t = [0] * len(f)
#                 t[len(f)-1] = jj
                
#                 z = hsurf_F2(I*I, l, t, f, j=i, s=1, s2=1) # hsurf_g3d(l, t, f, I, j=i, s=1)
                
#                 if not np.isnan(z):
#                     t[i] = z
#                 else:
#                     t[i] = 0.5/l
#                     z = hsurf_F2(I*I, l, t, f, j=i, s=1, s2=1) # hsurf_g3d(l, t, f, I, j=i-1, s=1)
#                     if not np.isnan(z):
#                         t[i-1] = z
#                     else:
#                         t[i-1] = 0.5/l
#                 xyz.append(t)
        
#         a = np.array(xyz)
        
#         if minmax==0:
#             #a[:,0:1] = np.min(a[:,0:1])
#             a[:,0] = np.min(a[:,0], axis=0)
#         else:
#             #a[:,0:1] = np.max(a[:,0:1])
#             a[:,0] = np.unique(a)[-2] if np.unique(a)[-1] == 0.5/l else  np.unique(a)[-1]
            
#         return a, np.array(xyz)
    
#     def getxyz(l, f, gi):
#         xyz = []
        
#         for i in range(len(f)-1):
#             for jj in [0, 1/(2*l)]:
#                 t = [0]*len(f)
#                 t[len(f)-1] = jj
#                 z = hsurf_g(l, t, f, gi, j=i, s=1)
#                 #z = hsurf_F2a(gi, l, t, f, j=i, s=1, s2=1)
                
#                 if ~np.isnan(z):
#                     t[i] = z
#                 else:
#                     t[i] = 0.5/l
#                     #z = hsurf_g3d(l, t, f, gi, j=i, s=1)
#                     z = hsurf_F2(gi*gi, l, t, f, j=i, s=1, s2=1)
                    
#                     if not np.isnan(z):
#                         t[i-1] = z
#                     else:
#                         t[i-1] = 0.5/l
                                
#                 xyz.append(t)
        
#         a = np.copy(xyz)
#         a = a[a[:, 1].argsort()][::-1]
#         #print("\n a - before :: \n", a)
#         a[0][:2] = a[1][:2] if np.all(a[1][:2] <= a[0][:2]) else a[0][:2]
#         a[1][:2] = a[0][:2] if np.all(a[1][:2] >= a[0][:2]) else a[1][:2]
        
#         a[2][:2] = a[3][:2] if np.all(a[3][:2] <= a[2][:2]) else a[2][:2]
#         a[3][:2] = a[2][:2] if np.all(a[3][:2] >= a[2][:2]) else a[3][:2]
        
#         #print("\n a - after :: \n", a)
#         return np.array(a), np.array(xyz)    
    
#     def get_pnt_nEPA(l, f, I):
#         tot_coor=[]
#         for i in range(len(f)):
#             tem_coor    = np.zeros(len(f))
#             tem_coor[i] = hsurf_F2(I*I, l, tem_coor, f, j=i, s=1, s2=1, nan=False) #hsurf_g3d(l, tem_coor, f, gi, i, s=1)
#             tot_coor.append(tem_coor)
            
#         return tot_coor
    
#     def gt3(l, f, I):
#         xyz = []
        
#         t = [0] * len(f)
#         z = hsurf_F2(I*I, l, t, f, j=0, s=1, s2=1, nan=False)
        
#         if not np.isnan(z):
#             t[0] = z
#         else:
#             t[0] = 0.5/l
#         ###
#         xyz.append(t)
        
#         for i in range(1,len(f)):
#             t = [0] * len(f)
#             t[i] = 0.5/l
#             z = hsurf_F2(I*I, l, t, f, j=0, s=1, s2=1, nan=False)
            
#             if not np.isnan(z):
#                 t[0] = z
#             else:
#                 t[0] = 0.5/l
            
#             xyz.append(t)
        
#         t = [0.5/l] * len(f)
#         z = hsurf_F2(I*I, l, t, f, j=0, s=1, s2=1, nan=False)
        
#         if not np.isnan(z):
#             t[0] = z
#         else:
#             t[0] = 0.5/l
#         ###
#         xyz.append(t)
        
#         a = np.array(xyz) ; a[:,0] = np.unique(a)[-1] if np.unique(a)[-1] == 0.5/l else  np.unique(a)[-1]
#         return xyz, a
    
    
#     # ---> Main func starts
    
#     pnt = get_pnt_nEPA(l, f, gi)
#     count_Nan = np.count_nonzero(np.isnan(pnt))
    
#     if   count_Nan == 0:
#         f = np.array(f)
#         k   = 2*np.pi*l ; xp = (1/k)*np.arccos(gi/np.sum(f))
#         p4  = [xp]*len(f)
        
#         pnt = np.vstack([np.array(pnt),p4])
        
#         ni, di   = plane_eq2n(pnt)
        
#         n_select     = ni
#         distance_all = remaining(l, gi, f, n_select, di)
        
#     elif count_Nan == 1:
#         pi, pext    = getxyz(l, f, gi) # getxyz_opt3(l, f, gi, minmax=min) # 
        
#         n_ext, dext = plane_eq2n(pext)
#         ni, di      = plane_eq2n(pi)
#         n_select    = ni #n_ext
#         distance_all = remaining(l, gi, f, n_select, di)
        
#     elif count_Nan > 1:        
#         #pi, pext  = getxyz_opt3(l, f, gi, minmax=1)
#         pi, pext    = gt3(l, f, gi)
        
#         n_ext, dext = plane_eq2n(pext)
#         ni, di      = plane_eq2n(pi)
#         n_select    = ni # n_ext
        
#         di = di if di !=0 and di<dext else dext       
#         distance_all = remaining(l, gi, f, n_select, di)
        
#     else:
#         print("--> count_Nan ", count_Nan)
#         print("--> Wired isosurface check linearization routine ")
#         print("--> predicted points are : \n", pnt)
        
#     return n_select, distance_all
