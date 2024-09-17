import numpy as np
import polytope as pc
from itertools import permutations


# -------------------------------------------------------
# ============== Common Modules ==============
# -------------------------------------------------------

def getmesh(l: int, coordinates: list, imax:int =0.5) -> list:
    
    c = np.linspace(0,imax,int(2*l*imax+1) )
    
    k = [c, c]*len(coordinates)
    k = k[0:len(coordinates)]
    
    j = np.meshgrid(*k)
    
    [*dim] = np.shape(j)
    
    f1=(np.array([j[i].reshape(-1,1) for i in range([*dim][0])]))
    f2=np.hstack([f1[i] for i in range([*dim][0])])
    
    meshlist=np.array(f2)

    plist = []
    for meshid in meshlist:
        oo=np.cos(2*np.pi*l*meshid)
        if (np.all(np.sign(oo) == 1) or np.all(np.sign(oo) == -1)):
            plist.append(meshid)
       
    return np.array(plist)

def getsigncombination(r:int) -> list:
    scom=[]

    for i in range(0, r+1):
        t = [-1]*i+[1]*(r-i)
        w = set(permutations(t))
        for u in w:
            scom.append(u)
    return np.array(scom)

def  getpolytope(l: int, normal: list, distance: list, imax: int =0.5) -> list:
    """ returns a collection of polytope for given reflection l.
        The polytope parameters boundary distance and normal are
        the required inputs. This module is for both EPA and nEPA
        Also this do not assume I or G. So it will return 2*(l**m)
        polytope for given l. m is dimension of PS
    Args:
        l (int): reflection
        normal (list): direction of polytope
        distance (list): boundary distance
        imax (int, optional): Limit of PS. Defaults to 0.5.

    Returns:
        _type_: Region of polytope.
    """    
    
    polylist = []
    
    dlist = getmesh(l, normal, imax=0.5)
    
    scom  = getsigncombination(len(normal))
    scom  = scom[scom[:,len(normal)-1].argsort()][::-1]
    
    gpsc  = np.identity(len(normal))
    Apsc  = np.array(np.vstack([-gpsc, gpsc]))
    bpsc  = np.array([0]*len(normal) + [0.5]*len(normal))
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


# -------------------------------------------------------
# ============== Modules for EPA ==============
# -------------------------------------------------------
def getpolytope_EPA( l, normal, distance, amplitudesign, IorG='amplitude', imax=0.5):
    
    # temp=[]
    # for i in range(len(normal)+1):
    #     zero=np.zeros(len(normal)) 
    #     zero[0:i]=0.5
    #     temp.append(zero)

    temp = np.tril( np.ones(shape=(len(normal), len(normal))) , 0 )
    temp = imax*np.vstack([[0]*len(normal), temp])
    asym = pc.qhull(np.array(temp))
       
    polylist = []
    
    dlist1 = getmesh(l, normal, imax=0.5)
    
    if l==1:
        dlist = np.delete(dlist1, 1, 0)
    else:
        #dinx = [ dlx for cc, dlx in enumerate(dlist1) if dlx in asym]
        #dlist=np.array(dinx)
        dlist=dlist1
    
    scom  = getsigncombination(len(normal))
    scom  = scom[scom[:,len(normal)-1].argsort()][::-1]
    
    gpsc  = np.identity(len(normal))
    Apsc  = np.array(np.vstack([-gpsc, gpsc]))
    bpsc  = np.array([0]*len(normal) + [0.5]*len(normal))
    psc   = pc.Polytope(Apsc, bpsc)
    
    aa    = np.array(normal)
    bb    = np.array(distance)
    
    for d in dlist:
        d  = np.array(d)
        oo = np.cos(2*np.pi*l*d)
        if IorG == 'amplitude':
            if (np.all(np.sign(oo) == amplitudesign)):
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
                        if (w <= psc):
                            polylist.append(w)
                    
        elif IorG == 'intensity':
            
            if( np.all(np.sign(oo) == 1) or np.all(np.sign(oo) == -1) ):
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
                        if (w <= psc):
                            polylist.append(w)
        else:
            print("===> select correct option for IorG")
                        
    return pc.Region(polylist)


# -------------------------------------------------------
# ============== Modules for non EPA ==============
# -------------------------------------------------------

def getpolytope_nEPA( l, normal, distance, amplitudesign, IorG='amplitude', imax=0.5):

    polylist = []

    dlist = getmesh(l, normal, imax=0.5)

    if l==1:
        dlist = np.delete(dlist, 1, 0)
    else:
        pass

    scom  = getsigncombination(len(normal))
    scom  = scom[scom[:,len(normal)-1].argsort()][::-1]

    gpsc  = np.identity(len(normal))
    Apsc  = np.array(np.vstack([-gpsc, gpsc]))
    bpsc  = np.array([0]*len(normal) + [0.5]*len(normal))
    psc   = pc.Polytope(Apsc, bpsc)
        
    aa    = np.array(normal)
    bb    = np.array(distance)

    for d in dlist:
        d  = np.array(d)
        oo = np.cos(2*np.pi*l*d)
        if IorG == 'amplitude':
            if (np.all(np.sign(oo) == amplitudesign)):
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
                        if (w <= psc):
                            polylist.append(w)
                    
        elif IorG == 'intensity':
            
            if( np.all(np.sign(oo) == 1) or np.all(np.sign(oo) == -1) ):
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
                    #print(f"w.chebXc : {w.chebXc}")
                    if w.chebXc is not None:
                        if (w <= psc):
                            polylist.append(w)
        else:
            print("===> select correct option for IorG")
                        
    return pc.Region(polylist)


# ===> I do not know why i wrote this module. but thinking that if coordinates of linearization point
# #      is known then this module can be used
# def repeat(p, d, f, imin, imax):
    
#     pts =[]
#     inx =np.argwhere(d != 0)
#     nz  =np.count_nonzero(d)
    
#     if nz == 0:
#         e1=np.copy(p)
#         pts.append(e1)
        
#     if nz != 0:
#         r,c = np.shape(p)
        
#         if (nz != len(d)):
            
#             if (np.all((d[inx[:,0]]+p[:,inx[:,0]])>=imin) and np.all((d[inx[:,0]]+p[:,inx[:,0]])<=imax)):
                              
#                 if (nz == 1):
#                     e2=np.copy(p)
#                     e2[:,inx[:,0]]=e2[:,inx[:,0]]+d[inx[:,0]]
#                     pts.append(e2)
            
#         if (np.all((d[inx[:,0]]-p[:,inx[:,0]])>=imin) and np.all((d[inx[:,0]]-p[:,inx[:,0]])<=imax)):
            
#             e4=np.copy(p)
#             e4[:,inx[:,0]]=d[inx[:,0]]-p[:,inx[:,0]]
#             pts.append(e4)
            
#             if (nz >1):
#                 for j in f:
#                     e4a=np.copy(p)
#                     e4a=e4a*j
                    
#                     e4a[:,inx[:,0]]=d[inx[:,0]]-e4a[:,inx[:,0]]
                    
#                     if (np.all(e4a>=imin) and np.all(e4a<=imax)):
#                         pts.append(e4a)
        
#     return pts
