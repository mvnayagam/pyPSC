import numpy as np
import sys
import os
from scipy.optimize import fsolve

def find_interception(x,y,m):
    return y-m*x

def findp3x(x, m, g=1.5, h=1):
    return m+np.sin(2*np.pi*h*x)/np.sqrt(1- (g - np.cos(2*np.pi*h*x))**2) 

def findp3y(x,h,g):
    return 1/(2*np.pi*h)*np.arccos(g-np.cos(2*np.pi*h*x))

def intercept_jonas(a,b):
    if (len(a) != len(b)):
        print("a and b have differnet length. Check it")
        sys.exit()
    
    m1 = (a[1]-a[3])/(a[0]-a[2])
    n1 = a[1]-a[0]*m1
    
    m2 = (b[1]-b[3])/(b[0]-b[2])
    n2 =  b[1]-b[0]*m2
    
    xp = (n2-n1)/(m1-m2)
    yp = m2*xp + n2
    return(xp,yp)


def hsurf_g(h,x,f,gi,j,s=1):
    k = 2*np.pi*h
    for i in range(len(x)):
        argm = s*gi/f[j] - (f[i]/f[j])*(np.cos(k*x[i]))
    xj = (np.arccos(argm))/k
    return xj
    

# def jonaspnts(gi,l):
#     #### Jonas Area 
#     k    = 2*np.pi*l
#     gi   = np.abs(gi)
#     #### Finding point p1 and p2  
#     p1x  = (1/k)*np.arccos(gi/2) ;       #  x2* = x1*
#     p1y  = p1x
    
#     p2x  = (1/k)*np.arccos(gi-1)         #  x1* = pnt2 ; x2* = 0
#     p2y  = 0
    
#     m1   = (p2y-p1y)/(p2x-p1x)           # slope of First line
#     n1   = find_interception(p2x,p2y,m1)
        
#     #### Finding point p3, p4 and p5
#     xini = 0.0
#     p5x  = fsolve(findp3x,xini,args=(m1,gi,l), factor=1, epsfcn=1e-16, maxfev=10000000)[0]
#     p5y  = findp3y(p5x,l,gi)
    
#     n2   = find_interception(p5x,p5y,m1)
    
#     p4x  = -n2 / (m1-1)
#     p4y  = p4x
    
#     p3x  = -n2/m1
#     p3y  = 0
    
    
#     pnt  = np.array([[p1x, p1y], [p2x, p2y], [p3x, p3y], [p4x, p4y],[ p5x, p5y]])
    
#     return pnt
def writetofile(fn,dx,dy,p,gmin,gmax):
    
    fmts="%2.8f \t%2.8f \t%2.8f \t%2.8f \t%2.8f \t%2.8f \t%2.8f \t%2.8f \t%2.8f \t%2.8f\n"
    if (np.min(dx-p)>=gmin and np.max(dx-p)<=gmax and np.min(dy-p)>=gmin and np.max(dy-p)<=gmax):
        
        ### Black line
        fn.write(fmts%(dx-p[1],dy-p[0],dx-p[3],dy-p[2],dx-p[5],dy-p[4],dx-p[7],dy-p[6],dx-p[9],dy-p[8]))
        fn.write(fmts%(dy-p[0],dx-p[1],dy-p[2],dx-p[3],dy-p[4],dx-p[5],dy-p[6],dx-p[7],dy-p[8],dx-p[9]))
    
    if (np.min(dx+p)>=gmin and np.max(dx+p)<=gmax and np.min(dy+p)>=gmin and np.max(dy+p)<=gmax):
        
        ### Green line
        fn.write(fmts%(dx+p[0],dy+p[1],dx+p[2],dy+p[3],dx+p[4],dy+p[5],dx+p[6],dy+p[7],dx+p[8],dy+p[9]))
        fn.write(fmts%(dx+p[1],dy+p[0],dx+p[3],dy+p[2],dx+p[5],dy+p[4],dx+p[7],dy+p[6],dx+p[9],dy+p[8]))

    if (np.min(dx+p) >=gmin and np.max(dx+p)<=gmax and np.min(dy-p)>=gmin and np.max(dy-p)<=gmax):
        
        ### Yello line
        fn.write(fmts%(dx+p[1],dy-p[0],dx-p[3],dy-p[2],dx-p[5],dy-p[4],dx+p[7],dy-p[6],dx+p[9],dy-p[8]))
        fn.write(fmts%(dy-p[0],dx+p[1],dy-p[2],dx-p[3],dy-p[4],dx-p[5],dy-p[6],dx+p[7],dy-p[8],dx+p[9]))
    
    if (np.min(dx+p) >=gmin and np.max(dx-p)<gmax and np.min(dy-p) >=gmin and np.max(dy-p)<=gmax):
        
        ### magenta line
        fn.write(fmts%(dx+p[1],dy-p[0],dx+p[2],dy+p[3],dx+p[4],dy+p[5],dx+p[7],dy-p[6],dx+p[8],dy-p[9]))
        fn.write(fmts%(dy-p[0],dx+p[1],dy+p[3],dx+p[2],dy+p[5],dx+p[4],dy-p[6],dx+p[7],dy-p[9],dx+p[8]))

def jonaspnts_error(gi,l,error):
    #### Jonas Area 
    k    = 2*np.pi*l
    gi   = np.abs(gi)
    #### Finding point p1 and p2  
    p1x  = (1/k)*np.arccos(gi*(1+error)/2) ;       #  x2* = x1*
    p1y  = p1x
    
    p2x  = (1/k)*np.arccos(gi*(1+error)-1)         #  x1* = pnt2 ; x2* = 0
    p2y  = 0
    
    m1   = (p2y-p1y)/(p2x-p1x)           # slope of First line
    n1   = find_interception(p2x,p2y,m1)
        
    #### Finding point p3, p4 and p5
    xini = 0.0
    p5x  = fsolve(findp3x,xini,args=(m1,gi*(1-error),l), factor=1, epsfcn=1e-16, maxfev=10000000)[0]
    p5y  = findp3y(p5x,l,gi*(1-error))
    
    n2   = find_interception(p5x,p5y,m1)
    
    p4x  = -n2 / (m1-1)
    p4y  = p4x
    
    p3x  = -n2/m1
    p3y  = 0
    
    pnt  = np.array([p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y, p5x, p5y])
    
    return pnt

# def jonaspnts(gi,l):
#     #### Jonas Area 
#     k    = 2*np.pi*l
#     gi   = np.abs(gi)
#     #### Finding point p1 and p2  
#     p1x  = (1/k)*np.arccos(gi/2) ;       #  x2* = x1*
#     p1y  = p1x
    
#     p2x  = (1/k)*np.arccos(gi-1)         #  x1* = pnt2 ; x2* = 0
#     p2y  = 0
    
#     m1   = (p2y-p1y)/(p2x-p1x)           # slope of First line
#     n1   = find_interception(p2x,p2y,m1)
        
#     #### Finding point p3, p4 and p5
#     xini = 0.0
#     p5x  = fsolve(findp3x,xini,args=(m1,gi,l), factor=1, epsfcn=1e-16, maxfev=10000000)[0]
#     p5y  = findp3y(p5x,l,gi)
    
#     n2   = find_interception(p5x,p5y,m1)
    
#     p4x  = -n2 / (m1-1)
#     p4y  = p4x
    
#     p3x  = -n2/m1
#     p3y  = 0
    
#     pnt  = np.array([p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y, p5x, p5y])
    
#     return pnt



def jonaspnts_newa(gi,l,f):
    #### Jonas Area 
    k    = 2*np.pi*l
    gi   = np.abs(gi)
    #### Finding point p1 and p2  
    #p1x  = (1/k)*np.arccos(gi/2) ;       #  x2* = x1*
    #p1y  = p1x
    
    p1x  = (1/k)*np.arccos(gi-1)         #  x1* = pnt2 ; x2* = 0
    p1y  = 0
    
    p2y  = (1/k)*np.arccos(gi-1)         #  x1* = pnt2 ; x2* = 0
    p2x  = 0
    
    m1   = (p2y-p1y)/(p2x-p1x)           # slope of First line
    n1   = find_interception(p2x,p2y,m1)
        
    #### Finding point p3, p4 and p5
    xini = 0.0
    
    #p5x  = fsolve(findp3x,xini,args=(m1,gi,l), factor=1, epsfcn=1e-16, maxfev=10000000)[0]
    #p5y  = findp3y(p5x,l,gi)
    
    p5x = (1/k)*np.arccos(gi/np.sum(f))
    p5y = p5x #p4  = [xp]*len(f)
    
    n2   = find_interception(p5x,p5y,m1)
    
    p4y  = n2 #/ (m1-1)  # y = m1x+n2
    p4x  = 0 #p4x
    
    p3x  = -n2/m1
    p3y  = 0
    
    pnt  = np.array([[p2x, p2y], [p1x, p1y], [p3x, p3y], [p5x, p5y], [p4x, p4y]])
    
    return pnt


def getploy(h,points,imax=0.5):
    r1, r2 = 0, 0
    a  = []
    for i in range(1,h+1):
        r2=r2+i*i*(4*isos.max())**2
        aa = multistrip(int(r1), int(r2),points)
        
        try:
            aa = unary_union(aa)
        except:
            print("AssertionFailedException occured for RO h=", i, "trying with make_valid")
            aa = make_valid(aa)
            
        a.append(aa)
        r1=np.copy(r2)
        
    return (a)

def isosurfs_newa(h,xexp,f,j,fname):
    dlist=[]
    npts = 500
    isos = np.linspace(0., 0.5, npts)
    grid = [isos, isos]
    
    for l in range(1,h+1):
        gi    = np.abs(g(l, xexp, f))
        giso1 = hsurf_g(l, grid, f, gi, j, s=1)
        giso2 = hsurf_g(l, grid, f, gi, j, s=-1)
    
        pnts = jonaspnts_newa(gi,l, f)
        
        signcom  = fn_signcombination(len(xexp))
        meshlist = fn_mesh(l, xexp, isos.max())
            
        for meshid in meshlist:
            oo=np.cos(2*np.pi*l*meshid)
            if (np.all(np.sign(oo) == 1) or np.all(np.sign(oo) == -1)):
                d = np.array(meshid)
                dlist.append(meshid)
                plist=fn_repeat(pnts, d, signcom, isos.min(), isos.max())                
                
                if plist:
                    
                    fn_write(fname, plist)
                    plistn=np.flip(np.flip(plist, axis=1))
                    fn_write(fname,plistn)
        
    return()


def plotsnew(axsname,p,dx,dy,xmin,xmax,al,cc):
    
    if (np.min(dx-p)>=xmin and np.max(dx-p)<=xmax and np.min(dy-p)>=xmin and np.max(dy-p)<=xmax):
        
        axsname.plot([dx-p[1], dx-p[3]], [dy-p[0], dy-p[2]], '-', c=cc, alpha=al)
        axsname.plot([dx-p[5], dx-p[7]], [dy-p[4], dy-p[6]], '-', c=cc, alpha=al)
        
        axsname.plot([dy-p[0], dy-p[2]], [dx-p[1], dx-p[3]], '-', c=cc, alpha=al)
        axsname.plot([dy-p[4], dy-p[6]], [dx-p[5], dx-p[7]], '-', c=cc, alpha=al)
    
    if (np.min(dx+p)>=xmin and np.max(dx+p)<=xmax and np.min(dy+p)>=xmin and np.max(dy+p)<=xmax):
        axsname.plot([dx+p[0], dx+p[2]], [dy+p[1], dy+p[3]], '-', c=cc, alpha=al)
        axsname.plot([dx+p[4], dx+p[6]], [dy+p[5], dy+p[7]], '-', c=cc, alpha=al)
        
        axsname.plot([dy+p[1], dy+p[3]], [dx+p[0], dx+p[2]], '-', c=cc, alpha=al)
        axsname.plot([dy+p[5], dy+p[7]], [dx+p[4], dx+p[6]], '-', c=cc, alpha=al)
    
    if (np.min(dx+p) >=xmin and np.max(dx-p)<xmax and np.min(dy-p) >=xmin and np.max(dy-p)<=xmax):
        axsname.plot([dx+p[1], dx-p[3]], [dy-p[0], dy-p[2]], '-', c=cc, alpha=al)
        axsname.plot([dx-p[5], dx+p[7]], [dy-p[4], dy-p[6]], '-', c=cc, alpha=al)
        
        axsname.plot([dy-p[0], dy-p[2]], [dx+p[1], dx-p[3]], '-', c=cc, alpha=al)
        axsname.plot([dy-p[4], dy-p[6]], [dx-p[5], dx+p[7]], '-', c=cc, alpha=al)
        
    if (np.min(dx+p) >=xmin and np.max(dx+p)<=xmax and np.min(dy-p)>=xmin and np.max(dy-p)<=xmax):
        axsname.plot([dx+p[1], dx+p[2]], [dy-p[0], dy+p[3]], '-', c=cc, alpha=al)
        axsname.plot([dx+p[7], dx+p[4]], [dy-p[6], dy+p[5]], '-', c=cc, alpha=al)
        
        axsname.plot([dy-p[0], dy+p[3]], [dx+p[1], dx+p[2]], '-', c=cc, alpha=al)
        axsname.plot([dy+p[5], dy-p[6]], [dx+p[4], dx+p[7]], '-', c=cc, alpha=al)