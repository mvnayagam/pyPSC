import polytope as pc
import numpy as np
import intvalpy as ip
from itertools import permutations


def fn_mesh_polytope(l, coordinates, imax):
    
    c = np.linspace(0,imax,int(2*l*imax+1) )
    
    k = [c, c]*len(coordinates)
    k = k[0:len(coordinates)]
    
    j = np.meshgrid(*k)
    
    [*dim] = np.shape(j)
    
    f1=(np.array([j[i].reshape(-1,1) for i in range([*dim][0])]))
    f2=np.hstack([f1[i] for i in range([*dim][0])])
    
    meshlist=np.array(f2)
    
    dl = []
    
    for meshid in meshlist:
        oo=np.cos(2*np.pi*l*meshid)
        if (np.all(np.sign(oo) == 1) or np.all(np.sign(oo) == -1)):
            d = np.array(meshid)
            dl.append(meshid)
    
    return np.array(dl)

def fn_signcom_polytope(r):
    f=[]
    
    for i in range(0, r+1):
        t = [-1]*i+[1]*(r-i)
        w = set(permutations(t))
        for u in w:
            f.append(u)
    return np.array(f)

def find_intersection(s, r):
    u=[]
    
    for count, i in enumerate(r):
        v = s & i
        
        if type(v) is pc.Polytope:
            if not pc.is_empty(v):
                u.append(v)
        elif type(v) is pc.Region:
            for k in v:
                if not pc.is_empty(k):
                    u.append(k)
    return pc.Region(u)

def  getpolytope_amplitude( l, normal, distance, amplitudesign, x=[1, 1, 1], imax=0.5):
    
    polylist = []
    
    dlist = fn_mesh_polytope(l, x, imax=0.5)
    scom  = fn_signcom_polytope(len(x))
    scom  = scom[scom[:,len(x)-1].argsort()][::-1]
    
    gpsc  = np.identity(len(x))
    Apsc  = np.array(np.vstack([-gpsc, gpsc]))
    bpsc  = np.array([0]*len(x) + [0.5]*len(x))
    psc   = pc.Polytope(Apsc, bpsc)
    
    aa    = np.array(normal)
    bb    = np.array(distance)
    
    for d in dlist:
        d  = np.array(d)
        oo = np.cos(2*np.pi*l*d)
        if (np.all(np.sign(oo) == amplitudesign)): # or np.all(np.sign(oo) == -1)):
            for i in scom:
                
                A = []
                A.append(-i*aa)
                A.append( i*aa)
                
                if i[len(x)-1]>0:
                    b=np.array(np.array([-i[len(x)-1], i[len(x)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                else:
                    b=np.array(np.array([i[len(x)-1], -i[len(x)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                # ---> inner
                iden = np.identity(len(x))
                for k in range(len(x)):
                    A=np.vstack([A,-1*iden[k]])
                
                de = d + (i-1)*(1/(4*l))
                b=np.append(b, -de)
                
                # ---> outter
                for k in range(len(x)):
                    A=np.vstack([A,iden[k]])
                    
                de = d + 1*(i+1)*(1/(4*l))
                b=np.append(b, de)
                
                w=pc.Polytope(np.array(A),np.array(b))
                
                if w.chebXc is not None:
                    if (w.chebXc in psc):
                        polylist.append(w)
    
    
    return pc.Region(polylist)


def  get_polytope_nd( l, normal, distance, x=[1, 1, 1], imax=0.5):
    
    polylist = []
    
    dlist = fn_mesh_polytope(l, x, imax=0.5)
    scom  = fn_signcom_polytope(len(x))
    scom  = scom[scom[:,len(x)-1].argsort()][::-1]
    
    gpsc  = np.identity(len(x))
    Apsc  = np.array(np.vstack([-gpsc, gpsc]))
    bpsc  = np.array([0]*len(x) + [0.5]*len(x))
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
                
                if i[len(x)-1]>0:
                    b=np.array(np.array([-i[len(x)-1], i[len(x)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                else:
                    b=np.array(np.array([i[len(x)-1], -i[len(x)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                # ---> inner
                iden = np.identity(len(x))
                for k in range(len(x)):
                    A=np.vstack([A,-1*iden[k]])
                
                de = d + (i-1)*(1/(4*l))
                b=np.append(b, -de)
                
                # ---> outter
                for k in range(len(x)):
                    A=np.vstack([A,iden[k]])
                    
                de = d + 1*(i+1)*(1/(4*l))
                b=np.append(b, de)
                
                w=pc.Polytope(np.array(A),np.array(b))
                
                if w.chebXc is not None:
                    if (w.chebXc in psc):
                        polylist.append(w)
                
    return pc.Region(polylist)


def get_mitd( l, normal, distance, scom, dlis, x=[1, 1, 1], imax=1/6):
    
    polylist = []
    
    gpsc  = np.identity(len(x))
    Apsc  = np.array(np.vstack([-gpsc, gpsc]))
    bpsc  = np.array([0]*len(x) + [imax]*len(x))
    psc   = pc.Polytope(Apsc, bpsc)
    
    #v=plot_poly3dnew(psc, axs, alpha=0.15, color ='C0')
    
    aa    = np.array(normal)
    bb    = np.array(distance)
    
    for d in dlis:
        d  = np.array(d)
        oo = np.cos(2*np.pi*l*d)
        if (np.all(np.sign(oo) == 1) or np.all(np.sign(oo) == -1)):
            for i in scom:
                
                A = []
                A.append(-i*aa)
                A.append( i*aa)
                
                if i[len(x)-1]>0:
                    b=np.array(np.array([-i[len(x)-1], i[len(x)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                else:
                    b=np.array(np.array([i[len(x)-1], -i[len(x)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                # ---> inner
                iden = np.identity(len(x))
                for k in range(len(x)):
                    A=np.vstack([A,-1*iden[k]])
                
                de = d + (i-1)*(1/(4*l))
                b=np.append(b, -de)
                
                # ---> outter
                for k in range(len(x)):
                    A=np.vstack([A,iden[k]])
                    
                de = d + 1*(i+1)*(1/(4*l))
                b=np.append(b, de)
                
                w=pc.Polytope(np.array(A),np.array(b))
                
                if w.chebXc is not None:
                    if (w.chebXc in psc):
                        polylist.append(w)
                
    return pc.Region(polylist)


def get_polytope_nd_mitd( l, normal, distance, dlist, x=[1, 1, 1], imax=0.5):
    
    polylist = []
    
    scom  = fn_signcom_polytope(len(x))
    scom  = scom[scom[:,len(x)-1].argsort()][::-1]
    
    gpsc  = np.identity(len(x))
    Apsc  = np.array(np.vstack([-gpsc, gpsc]))
    bpsc  = np.array([0]*len(x) + [0.5]*len(x))
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
                
                if i[len(x)-1]>0:
                    b=np.array(np.array([-i[len(x)-1], i[len(x)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                else:
                    b=np.array(np.array([i[len(x)-1], -i[len(x)-1]])*(bb + np.sum([i[kk]*aa[kk]*d[kk] for kk in range(len(d))])))
                
                # ---> inner
                iden = np.identity(len(x))
                for k in range(len(x)):
                    A=np.vstack([A,-1*iden[k]])
                
                de = d + (i-1)*(1/(4*l))
                b=np.append(b, -de)
                
                # ---> outter
                for k in range(len(x)):
                    A=np.vstack([A,iden[k]])
                    
                de = d + 1*(i+1)*(1/(4*l))
                b=np.append(b, de)
                
                w=pc.Polytope(np.array(A),np.array(b))
                
                if w.chebXc is not None:
                    if (w.chebXc in psc):
                        polylist.append(w)
                
    return pc.Region(polylist)