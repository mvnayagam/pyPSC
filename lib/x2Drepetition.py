import numpy as np
from itertools import permutations, combinations


def repeat2D(p, d, f, imin, imax):
    
    pts =[]
    inx =np.argwhere(d != 0)
    nz  =np.count_nonzero(d)
    
    if nz == 0:
        e1=np.copy(p)
        pts.append(e1)
        
    if nz != 0:
        r,c = np.shape(p)
        
        if (nz != len(d)):
            
            if (np.all((d[inx[:,0]]+p[:,inx[:,0]])>=imin) and np.all((d[inx[:,0]]+p[:,inx[:,0]])<=imax)):
                              
                if (nz == 1):
                    e2=np.copy(p)
                    e2[:,inx[:,0]]=e2[:,inx[:,0]]+d[inx[:,0]]
                    pts.append(e2)
            
        if (np.all((d[inx[:,0]]-p[:,inx[:,0]])>=imin) and np.all((d[inx[:,0]]-p[:,inx[:,0]])<=imax)):
            
            e4=np.copy(p)
            e4[:,inx[:,0]]=d[inx[:,0]]-p[:,inx[:,0]]
            pts.append(e4)
            
            if (nz >1):
                for j in f:
                    e4a=np.copy(p)
                    e4a=e4a*j
                    
                    e4a[:,inx[:,0]]=d[inx[:,0]]-e4a[:,inx[:,0]]
                    
                    if (np.all(e4a>=imin) and np.all(e4a<=imax)):
                        pts.append(e4a)
        
    return pts

def mesh(l, coordinates, imax):
    
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

def getsigncom(r):
    scom=[]

    for i in range(0, r+1):
        t = [-1]*i+[1]*(r-i)
        w = set(permutations(t))
        for u in w:
            scom.append(u)
    return np.array(scom)

def linrep_DS(h, f, pnt, meshgrid, IorG='intensity', signofIorG=1, imin=0, imax=0.5):
    
    plist   = []
    signcom = getsigncom(len(f))
    
    for i in signcom:
        pinner  = pnt*i
        plist.append(pinner)
        
        poutter = np.flip(pinner, axis=1)
        plist.append(poutter)
    
    pfinal = []
    
    for meshid in np.array(meshgrid):
        for ii in plist:
            if IorG == 'amplitude':
                oo=np.cos(2*np.pi*h*meshid)
                if (np.all( np.sign(oo) == signofIorG )):
                    ji = meshid+ii
                    if np.all(ji<=imax) and np.all(ji>=imin):
                        pfinal.append(ji)
            elif IorG == 'intensity':
                oo=np.cos(2*np.pi*h*meshid)
                if (np.all(np.sign(oo) == 1) or np.all(np.sign(oo) == -1)):
                    ji = meshid+ii
                    if np.all(ji<=imax) and np.all(ji>=imin):
                        pfinal.append(ji)
            else:
                print("please select correct option for IorG. It should be either intensity of amplitude")
    return pfinal

def linrep_SS(h, f, pnt, meshgrid, IorG='intensity', signofIorG=1, imin=0, imax=0.5):
        
    plist   = []
    signcom = getsigncom(len(f))
    
    for i in signcom:
        pinner  = pnt*i
        plist.append(pinner)
        
        #poutter = np.flip(pinner, axis=1)
        #plist.append(poutter)
    
    pfinal = []

    for meshid in np.array(meshgrid):
        for ii in plist:
            if IorG == 'amplitude':
                oo=np.cos(2*np.pi*h*meshid)
                if (np.all( np.sign(oo) == signofIorG )):
                    ji = meshid+ii
                    if np.all(ji<=imax) and np.all(ji>=imin):
                        pfinal.append(ji)
            elif IorG == 'intensity':
                oo=np.cos(2*np.pi*h*meshid)
                if (np.all(np.sign(oo) == 1) or np.all(np.sign(oo) == -1)):
                    ji = meshid+ii
                    if np.all(ji<=imax) and np.all(ji>=imin):
                        pfinal.append(ji)
            else:
                print("please select correct option for IorG. It should be either intensity of amplitude")
                    
    # for j in np.array(meshgrid):
        
    #     for ii in plist:
    #         ji = j+ii
    #         if np.all(ji<=imax) and np.all(ji>=imin):
    #             pfinal.append(ji)
    
    return pfinal

def writedata(fn, data):
    
    dimension, r, c = np.shape(data)
    
    for a in data:
        countr=0
        for i in a:
            countc=0
            for j in i:
                if countr <(r-1):
                    fn.write("%2.8f \t"%(j))
                else:
                    if countc <(c-1):
                        fn.write("%2.8f \t"%(j))
                    else:
                        fn.write("%2.8f\n"%(j))
                countc += 1
            countr += 1    
    return()
