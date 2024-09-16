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


# def checklinear(l: int, xcoor: list, f: list, normal: list, distance: list, j: int=2, n: int=20, s:int =1, testiso: bool=True):
    
#     """_Checks the quality of linearization. This is for EPA model not for non-EPA_
#     Args:
#         l (int)          : _The reflection order to be processed_
#         xcoor (list)     : _ Given atomic structure_
#         f (list)         : _atomic scattering factors it is actually [1.0]*len(xcoor)_
#         normal (list)    : _Found normal vector of isosurface of l_
#         distance (list)  : _Distance of inner and outer boundaries_
#         j (int, optional): _The atom index along last axis_. Defaults to 2.
#         n (int, optional): _Number of points to create isosurface_. Defaults to 20.
#         s (int, optional): _sign of amplitude_. Defaults to 1.
#         testiso (bool, optional): _Testing the isosurface_. Defaults to True.
#     """
     
#     # Create a linearly spaced array for the polytope's boundary
#     lspace = np.linspace(0, 1 / (2 * l), n)

#     # Generate the meshgrid for k-space dimensions
#     kz = np.meshgrid(*([lspace] * (len(f) - 1))) ; kz = list(kz)

#     # Compute the surface and isosurface points
#     gi = np.abs(g(l, xcoor, f))
#     gzp = hsurf_g(l, kz, f, gi, j=len(f)-1, s=s)

#     # Calculate the polytope
#     o = getpoly_mitd(
#         l, normal, distance,
#         scom=np.ones((1, len(f))),
#         dlist=np.zeros((1, len(f))),
#         imax=lspace.max()
#     )
    
#     # Flatten each 3D array in kz and the gzp array
#     kz_flattened = [kz_i.flatten() for kz_i in kz]
#     gzp_flattened = gzp.flatten()
    
#     # Stack the flattened arrays column-wise to get the desired 2D array
#     iso_grid = np.vstack(kz_flattened + [gzp_flattened]).T
#     valid_iso_grid = iso_grid[~np.isnan(iso_grid).any(axis=1)]
    
#     # Check for isosurface containment within the polytope if required
#     if testiso:
#         # kz.append(np.array(gzp))  # Add isosurface points
#         # tz = np.vstack(np.dstack(kz))  # Stack all arrays along the third axis
#         # Check if any point lies outside the polytope
        
#         outside_points = [ti for ti in valid_iso_grid if ti not in o]
        
#         if outside_points:
#             print(f"\x1b[1;31m--> Checking the quality of linearization process")
#             print(f"\n--> Found isosurface outside for the point at {len(outside_points)} locations")
#             print("\x1b[1;31m--> Check the linearization step <--\x1b[0m")
#             raise ValueError("\x1b[1;32m--> Exiting: Linearization failed as isosurface points are outside the polytope\x1b[0m")
#         else:
#             print("\x1b[1;32m--> Polytope contains complete isosurface. Successful Linearization for \x1b[1;31mRO = %g\x1b[0m" % l)
#     return

def checklinear(l: int, f: list, I: float, normal: list, distance: list, j: int=2, n: int=20, s:int =1, testiso: bool=True):
        
    """_Checks the quality of linearization. This is for EPA model not for non-EPA_
    Args:
        l (int)          : _The reflection order to be processed_
        xcoor (list)     : _ Given atomic structure_
        f (list)         : _atomic scattering factors it is actually [1.0]*len(xcoor)_
        normal (list)    : _Found normal vector of isosurface of l_
        distance (list)  : _Distance of inner and outer boundaries_
        j (int, optional): _The atom index along last axis_. Defaults to 2.
        n (int, optional): _Number of points to create isosurface_. Defaults to 20.
        s (int, optional): _sign of amplitude_. Defaults to 1.
        testiso (bool, optional): _Testing the isosurface_. Defaults to True.
    """
    
    # Create a linearly spaced array for the polytope's boundary
    lspace = np.linspace(0, 1 / (2 * l), n)

    # Generate the meshgrid for k-space dimensions
    kz = np.meshgrid(*([lspace] * (len(f) - 1))) ; kz = list(kz)

    # Compute the surface and isosurface points
    gi = I #np.abs(g(l, xcoor, f))
    gzp = hsurf_g(l, kz, f, gi, j=len(f)-1, s=s)
    
    #; print(f'I={I} f={f} l={l} normal: {normal}')
    
    # Calculate the polytope
    o = getpoly_mitd(l, normal, distance, scom=np.ones((1, len(f))), dlist=np.zeros((1, len(f))), imax=lspace.max())
    
    # Flatten each 3D array in kz and the gzp array
    kz_flattened = [kz_i.flatten() for kz_i in kz]
    gzp_flattened = gzp.flatten()
    
    # Stack the flattened arrays column-wise to get the desired 2D array
    iso_grid = np.vstack(kz_flattened + [gzp_flattened]).T
    valid_iso_grid = iso_grid[~np.isnan(iso_grid).any(axis=1)]
    
    dx=np.dot(valid_iso_grid, normal)
    
    #print(f'from checklineara d_min & d_max : {np.min(dx)} {np.max(dx)}')
    #print(f'from checklineara location are  : {valid_iso_grid[np.where(dx<=np.min(dx))][0]} {valid_iso_grid[np.where(dx==np.max(dx))]}\n--------------------')
    
    # Check for isosurface containment within the polytope if required
    if testiso:
        # Check if any point lies outside the polytope
        outside_points = np.array([ti for ti in valid_iso_grid if ti not in o])
        
        if len(outside_points)!=0:
            #ds=np.dot(outside_points, normal)
            #print(f"\n\x1b[1;31m--> Checking the quality of linearization process")
            #print(f"--> Found isosurface outside for the point at {len(outside_points)} locations.")# ds are ---> {np.min(dx)} {np.max(dx)}")
            #print("\x1b[1;31m--> Check the linearization step <--\x1b[0m")
            dt_corrected = [np.min(dx), np.max(dx)]
            status = False
            #raise ValueError("\x1b[1;32m--> Exiting: Linearization failed as isosurface points are outside the polytope\x1b[0m")
        else:
            print("\x1b[1;32m--> Polytope contains complete isosurface. Successful Linearization for \x1b[1;31mRO = %g\x1b[0m" % l)
            dt_corrected = []
            status = True
    return status, dt_corrected


def checklinearplot(l, xexp, f, normal, distance, j=2,n=100, s=1, testiso=True, plot=True):
    
    lspace  = np.linspace(0, 1/(2*l), n)
    kj = [lspace]*(len(f)-1)
    kz = np.meshgrid(*kj)
    gz = np.zeros_like(kz[0])
    j = len(f)-1
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


def checklinear_I(l, I, f, normal, distance, n=100, s=1):
    
    j = len(f)-1
    lspace  = np.linspace(0, 1/(2*l), n)
    kj = [lspace]*(len(f)-1)
    kz = np.meshgrid(*kj)
    #gz = np.zeros_like(kz[0])
    
    gzp = hsurf_F(I, l, [*kz], f, j, s=1, s2=1)
    o   = getpoly_mitd(l, normal, distance, scom=np.ones((1, len(f))), dlist=np.zeros((1, len(f))), imax=lspace.max() )
    
    # Flatten each 3D array in kz and the gzp array
    kz_flattened = [kz_i.flatten() for kz_i in kz]
    gzp_flattened = gzp.flatten()
    
    # Stack the flattened arrays column-wise to get the desired 2D array
    iso_grid = np.vstack(kz_flattened + [gzp_flattened]).T
    valid_iso_grid = iso_grid[~np.isnan(iso_grid).any(axis=1)]
    
    check=[i in o for i in valid_iso_grid]
        
    if not np.all(check):
        index = np.where(~np.array(check))[0]
        dr  = [np.dot(normal,valid_iso_grid[inx]) for inx in index]
        return False, [np.min(dr), np.max(dr)]
    else:
        return True

# Old function. delte may later
# def checklinear(l: int, xcoor: list, f: list, normal: list, distance: list, j: int=2, n: int=100, s:int =1, testiso: bool=True):
#     """_Checks the quality of linearization. This is for EPA model not for non-EPA_
#     Args:
#         l (int)          : _The reflection order to be processed_
#         xcoor (list)     : _ Given atomic structure_
#         f (list)         : _atomic scattering factors it is actually [1.0]*len(xcoor)_
#         normal (list)    : _Found normal vector of isosurface of l_
#         distance (list)  : _Distance of inner and outer boundaries_
#         j (int, optional): _The atom index along last axis_. Defaults to 2.
#         n (int, optional): _Number of points to create isosurface_. Defaults to 100.
#         s (int, optional): _sign of amplitude_. Defaults to 1.
#         testiso (bool, optional): _Testing the isosurface_. Defaults to True.
#     """    
#     lspace  = np.linspace(0, 1/(2*l), n)
#     kj = [lspace]*(len(f)-1)
#     kz = np.meshgrid(*kj)
#     #gz = np.zeros_like(kz[0])
    
#     gi  = np.abs(g(l, xcoor, f))
#     j   = len(f)-1 
#     gzp = hsurf_g(l, [*kz], f, gi, j, s=1)
    
#     o = getpoly_mitd(l, normal, distance, scom=np.array([[1]*len(f)]), dlist=np.array([[0]*len(f)]), imax=lspace.max())
    
#     if testiso:
#         kz.extend([np.array(gzp)])
#         tz = np.vstack(np.dstack([*kz]))
        
#         for ti in tz:
#             if np.all(~np.isnan(ti)):
#                 if ti in o[0]:
#                     continue
#                 else:
#                     print("\x1b[1;31m--> Checking the quality of linearization process ",end=" ")
#                     print("\n--> Found isosurface outside for the point on the location ",ti)
#                     print("--> The point is: ", ti,"isoutside the polytope \n")
#                     print("\x1b[1;31m--> Check the linearization step <-- \x1b[0m")
#                     print("\x1b[1;32m--> I am quitting hier, BYE <-- \x1b[0m")
#                     sys.exit()
#         print("\x1b[1;32m--> Polytope contains complete isosurface. Successful Linearization for \x1b[1;31mRO = %g\x1b[0m"%(l))
    
#     return 