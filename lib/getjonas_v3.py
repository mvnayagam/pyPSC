import numpy as np
import sys, os, re
# import collections
from itertools import permutations, combinations
from pathlib import Path

import shapely
from scipy.optimize import fsolve
from shapely.geometry import Polygon, mapping, Point
from descartes import PolygonPatch
from shapely.validation import make_valid
from shapely.ops import unary_union

def findpy(x,h,gi,f):
    return (1/(2*np.pi*h))*np.arccos(gi/f[1] - (f[0]/f[1])*np.cos(2*np.pi*h*x))

def findpx(y,h,gi,f):
    return (1/(2*np.pi*h))*np.arccos(gi/f[0] - (f[1]/f[0])*np.cos(2*np.pi*h*y))

def find_interception(x,y,m):
    return y-m*x

def findp3x(x, m, gi, h):
    return m+np.sin(2*np.pi*h*x)/np.sqrt(1- (gi - np.cos(2*np.pi*h*x))**2) 

def findp3y(x,h,gi):
    return 1/(2*np.pi*h)*np.arccos(gi-np.cos(2*np.pi*h*x))

def jonaspnts_error(gi, l, xexp, f, error=0, printres=False):#
    
    #### Jonas Area 
    k    = 2*np.pi*l ;  gi   = np.abs(gi)
    
    #### Finding point p1 and p2  
    p1x  = (1/k)*np.arccos(gi*(1+error)/2) ;       #  x2* = x1*
    p1y  = p1x
    
    p2x  = (1/k)*np.arccos(gi*(1+error)-1)         #  x1* = pnt2 ; x2* = 0
    p2y  = 0
    
    m1   = (p2y-p1y)/(p2x-p1x)           # slope of First line
    n1   = find_interception(p2x,p2y,m1)
        
    #### Finding point p3, p4 and p5
    xini = 0.0
    #p5x  = fsolve(findp3x,xini,args=(m1,gi*(1-error),l), factor=1, epsfcn=1e-16, maxfev=10000000)[0]
    #p5y  = findp3y(p5x,l,gi*(1-error))
    
    j = 1
    p5x, p5y = fn_solveforx_v2(l, gi, f, m1, j, xexp) #fn_solveforx_v2(l, gi, f, m1, xexp, j=1)
    
    n2   = find_interception(p5x,p5y,m1)
    
    p4x  = -n2 / (m1-1)
    p4y  = p4x
    
    p3x  = -n2/m1
    p3y  = 0
    
    #pnt  = np.array([p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y, p5x, p5y])
    pnt  = np.array([[p1x, p1y], [p2x, p2y], [p3x, p3y], [p5x, p5y], [p4x, p4y]])
    
    if printres:
        print("m1 is       :: ", m1)
        print("n1 is       :: ", n1)
        print("n2 is       :: ", n2)
        #print("p1x and p1y :: ", p1x, p1x)
        #print("p4x and p4y :: ", p4x, p4x)
        #print("p4x-p1x     :: ", p4x-p1x)
        #print("p4y-p1y     :: ", p4y-p1y)
        print("-1/m is     :: ",-1/m1)
        #print("p2x-p3x     :: ", p2x-p3x)
        #print("n2-n1       :: ", n2-n1)
        print("(p2y-p1y)/(p2x-p1x) :: ", (p2y-p1y), (p2x-p1x))
    return pnt


def fn_solveforx_v2(l, g, f, m, j, x):
    
    k = 2 * np.pi * l
    
    i = list(range(j)) + list(range(j+1,len(x)))
    
    a = (1-m*m)/(f[j]*f[j])
    b = 2 * m*m * g /(f[j]*f[j])
    c = m*m*( 1 - (g*g) / (f[j]*f[j]) ) - np.array([ (f[ii]*f[ii]) / (f[j]*f[j]) for ii in i]).sum(axis = 0) 
    
    
    if a != 0:
        
        z1 = (-b + np.sqrt(b*b - 4*a*c))/(2*a)
        z2 = (-b - np.sqrt(b*b - 4*a*c))/(2*a)
        
        r1=(1/k)*np.arccos(z1/f[0])
        r2=(1/k)*np.arccos(z2/f[0])
        
        if ~np.isnan(r1) and np.isnan(r2):
            ry=findpy(r1,l,g,f)
            return np.array([r1, ry])
        
        elif np.isnan(r1) and ~np.isnan(r2):
            ry=findpy(r2,l,g,f)
            return np.array([r2, ry])
        
        elif np.isnan(r1) and np.isnan(r2) :
            
            return np.array([float("nan")])
        
        elif ~np.isnan(r1) and ~np.isnan(r2):
            
            r1y=findpy(r1,l,g,f)
            r2y=findpy(r2,l,g,f)
            
            print("Multiple solution: r1, r1y= ", r1,r1y, " r2, r2y = ", r2, r2y)
            return np.array([ [r1, r1y], [r2, r2y] ])
        
        else:
            print("I obtained weird solution. I do not know what to return r1= ", r1," r2= ", r2)
            return np.array([r1, r2])
    
    else:        
        z1 = (-1*c/b)
        
        r1=(1/k)*np.arccos(z1/f[0])
        
        if ~np.isnan(r1):
            ry=findpy(r1,l,g,f)
            return np.array([r1, ry])
        elif np.isnan(r1):
            prx = (1/k)*np.arccos(g*(1+error)/np.sum(f))
            pry = prx
            return np.array([prx, pry])
        else:
            prx = (1/k)*np.arccos(g*(1+error)/np.sum(f))
            pry = prx
            return np.array([prx, pry])


def jonopt_error_nEPA(g, l, f, xcoor, j=1, error=0):
    pnt = []
    k   = 2*np.pi*l
    g   = np.abs(g) 
    
    p1x = (1/k)*np.arccos(g*(1+error)/np.sum(f))
    p1y = p1x
    pnt.append([p1x, p1y])
    
    p2x = (1/k)*np.arccos((g*(1+error)-f[1])/f[0])
    p2y = 0
    
    if np.isnan(p2x):
        p2x = 0.5/l
    
    pnt.append([p2x, p2y])
    
    m1  = (p2y-p1y)/(p2x-p1x)
    n1  = find_interception(p2x,p2y,m1)
    
    p5   = fn_solveforx_v2(l, g, f, m1, j, xcoor)
    
    if ~np.all(np.isnan(p5)):
        
        if len(np.shape(p5)) > 1 :
            
            pnt = []
            
            if np.floor(p5[0,1]/p5[0,0]) <= 1 and np.floor(p5[1,1]/p5[1,0]) > 1:
                pLB=p5[0]
                pUB=p5[1]
            else:
                pLB=p5[1]
                pUB=p5[0]
            
            ### Lower Boundary (LB)
            nLB = pLB[1] - m1*pLB[0]
            
            pLB1y = 0.5/l
            pLB1x = (pLB1y - nLB)/m1
            
            pLB3x = np.abs(nLB/m1)
            pLB3y = 0
            #pLB3x = (1/k)*np.arccos((g*(1+error)-f[0])/f[1])
                        
            ### Upper Boundary (UB)
            
            nUB = pUB[1] - m1*pUB[0]
            pUB4x = np.abs(nUB/m1)
            pUB4y = 0
            
            pUB6y = 0.5/l
            pUB6x = (pUB6y - nUB)/m1
            
            ### Colloecting point in order
            pnt.append([pLB1x,  pLB1y])
            pnt.append([pLB[0], pLB[1]])
            pnt.append([pLB3x,  pLB3y])
            
            pnt.append([pUB4x,  pUB4y])
            pnt.append([pUB[0], pUB[1]])
            pnt.append([pUB6x,  pUB6y])
                        
            return np.array(pnt)
        
        elif ~np.all(np.isnan(p5)) and len(np.shape(p5)) == 1 :
            p5y = p5[1]
            p5x = p5[0]
        else:
            p5y = 0.5/l
            p5x  = (1/k)*np.arccos( g/f[0] - (f[1]/f[0])*np.cos(k*p5y) )
    
    n2  = find_interception(p5x,p5y,m1)
    
    p4x = -n2 / (m1-1)
    p4y = p4x
    #print("===> p4 is ", p4x, p4y)
    
    p3x = -n2 / m1
    p3y = 0
    
    pnt.append([p3x, p3y])
    pnt.append([p5x, p5y])
    
    
    #---> Part 2 linearization 
    
    p6x  = 0
    p6y  = (1/k)*np.arccos((g*(1+error)-f[0])/f[1])
    
    if np.isnan(p6y):
        p6y = 0.5/l
        
    m3   = (-p6y+p1y)/(-p6x+p1x)
    
    p7 = fn_solveforx_v2(l, g, f, m3, j, xexp)
    
    if ~np.all(np.isnan(p7)):
        p7x = p7[0]
        p7y = p7[1]
        
        n4   = find_interception(p7x,p7y,m3)
        p8x  = 0
        if ~np.isnan(n4) and n4<=0.5/l:
            p8y  = n4 
        elif ~np.isnan(n4) and n4>=0.5/l:
            p8y  = 0.5/l
        else:
            print("--> from def jonopt_error_v5: do not know what to do for p8y ")
        
    else:
        p7y = 0.5/l
        p7x  = (1/k)*np.arccos( g/f[0] - (f[1]/f[0])*np.cos(k*p7y) )  # === p7x  = findpx(p7y,l,g,f)
        
        n4   = find_interception(p7x,p7y,m3)
        p8x  = 0
        p8y  = 0.5/l
        
        
    n4   = find_interception(p7x,p7y,m3)
    
    p9x = n4/(1-m3)
    p9y = p9x
    
    if p9y == p4y and p9x == p4x :
        pnt.append([p4x, p4y])
        pnt.append([p9x, p9y])
    else:
        if p9y > p4y and p9x > p4x :
                    
            p9x = (n2-n4) / (m3-m1)
            p9y = m3 * p9x + n4
            
            pnt.append([p4x, p4y])
            pnt.append([p9x, p9y])
            
            m49   = (p9y-p4y)/(p9x-p4x)
            n49   = p1y - m49*p1x
            #pnewx = fn_solveforx(l, g, f, m49, j, xexp) #p1x + ( (0.5/l) - p1y ) / m49
            pnewy = p9y #(1/k)*np.arccos( g/f[1] - (f[0]/f[1])*np.cos(k*pnewx) )
            pnewx = (pnewy - n49)/( m49 )
            
            #pnt.append([pnewx, pnewy])
            
        else:
            
            p9x = (n2-n4) / (m3-m1)
            p9y = m3 * p9x + n4
            
            #pnt.append([p4x, p4y])
            pnt.append([p9x, p9y])
        
            m49   = (p9y-p4y)/(p9x-p4x)
            n49   = p1y - m49*p1x
            pnewy = p9y
            pnewx = (pnewy - n49)/( m49 )
            
            #pnt.append([pnewx, pnewy])
                
    pnt.append([p7x, p7y])
    pnt.append([p8x, p8y])
    pnt.append([p6x, p6y])
    pnt.append([p1x, p1y])
    
    
    return np.array(pnt)



def jonas4ring(gi,l,f,error=0):
    #### Jonas Area 
    k    = 2*np.pi*l
    gi   = np.abs(gi)
    #### Finding point p1 and p2  
    #p1x  = (1/k)*np.arccos(gi/2) ;       #  x2* = x1*
    #p1y  = p1x
    
    p1x  = (1/k)*np.arccos(gi*(1+error)-1)         #  x1* = pnt2 ; x2* = 0
    p1y  = 0
    
    p2y  = (1/k)*np.arccos(gi*(1+error)-1)         #  x1* = pnt2 ; x2* = 0
    p2x  = 0
    
    m1   = (p2y-p1y)/(p2x-p1x)           # slope of First line
    n1   = find_interception(p2x,p2y,m1)
        
    #### Finding point p3, p4 and p5
    xini = 0.0
    
    #p5x  = fsolve(findp3x,xini,args=(m1,gi,l), factor=1, epsfcn=1e-16, maxfev=10000000)[0]
    #p5y  = findp3y(p5x,l,gi)
    
    p5x = (1/k)*np.arccos(gi*(1-error)/np.sum(f))
    p5y = p5x #p4  = [xp]*len(f)
    
    n2   = find_interception(p5x,p5y,m1)
    
    p4y  = n2 #/ (m1-1)  # y = m1x+n2
    p4x  = 0 #p4x
    
    p3x  = -n2/m1
    p3y  = 0
    
    pnt  = np.array([[p2x, p2y], [p1x, p1y], [p3x, p3y], [p5x, p5y], [p4x, p4y]])
    
    return pnt


def fn_isosurfs(h,xexp,f,j,fname):
    dlist=[]
    npts = 500
    isos = np.linspace(0., 0.5, npts)
    grid = [isos, isos]
    
    for l in range(1,h+1):
        gi    = np.abs(g(l, xexp, f))
        giso1 = hsurf_g(l, grid, f, gi, j, s=1)
        giso2 = hsurf_g(l, grid, f, gi, j, s=-1)
    
        pnts = jonaspnts_error(gi,l)
                
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

def fn_isosurfs_4ring(h,xexp,f,j,fname):
    dlist=[]
    npts = 500
    isos = np.linspace(0., 0.5, npts)
    grid = [isos, isos]
    
    for l in range(1,h+1):
        gi    = np.abs(g(l, xexp, f))
        giso1 = hsurf_g(l, grid, f, gi, j, s=1)
        giso2 = hsurf_g(l, grid, f, gi, j, s=-1)
    
        pnts = jonas4ring(gi,l,f)
        
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

def multistrip(a,b,pnts):
    d=[]
    for i in range(a,b,1):
        df=int(len(pnts[i]))
        d1=Polygon([(pnts[i][j], pnts[i][j+1]) for j in range(0,df,2)])
        d.append(d1)
    return shapely.geometry.MultiPolygon([poly for poly in d])


def fn_getploygons(h,points,imax=0.5):
    
    r1, r2 = 0, 0
    a  = []
    
    for i in range(1,h+1):
        #r1 = r1+4*(i-1)**2
        #r2 = r2+4*i*i
        
        r2=r2+i*i*(4*imax)**2
        aa = multistrip(int(r1), int(r2),points)
        
        try:
            aa = unary_union(aa)
        except:
            print("AssertionFailedException occured for RO h=", i, "trying with make_valid")
            aa = make_valid(aa)
            
        a.append(aa)
        r1=np.copy(r2)
        
    return (a)

def fn_getintersections(h,a,xexp,fname,count):
    
    s  = []
    
    for j in range(h-1):
        #print("Doing for j's upto :: ", j+1," with j = ",j+2)
        try:
            if j == 0:
                ss = a[j].intersection(a[j+1])
            else:
                ss = s[-1].intersection(a[j+1])
        except:
            fname.write('Pair-{} : TopologyException error for x1 = {:2.4} and x2 = {:2.4} at h = {}\n'.format(count,xexp[0], xexp[1], (j+1)))
            continue
        
        if not ss:
            #print("===> ss is empty for j = ", j+2)
            ss=s[-1]
        
        s.append(ss)
        
    return (s, j)


def fn_writepolygons(fname, polys):
    
    for i in polys:
        x, y = i.exterior.coords.xy
        
        for xl in range(len(x)):
            #fname.write('{:10.10}\t\t{:10.10}\t\t'.format(x[xl],y[xl]))
            fname.write("%2.12f\t\t%2.12f\t\t"%(x[xl],y[xl]))
        fname.write("\n")
    
    return ()


def isInside(p, v1=np.array([0.0, 0.0]), v2=np.array([0.5,0.0]), v3=np.array([0.25, 0.25])):
    
    def get_area(vert1, vert2, vert3):
        veca = vert2-vert1
        vecb = vert3-vert1
        return 0.5*np.abs(np.cross(veca, vecb))
    
    A = get_area (v1, v2, v3)
    A1 = get_area (p, v2, v3)
    A2 = get_area (v1, p, v3)
    A3 = get_area (v1, v2, p)
    
    if(A >= A1 + A2 + A3):
        return True
    else:
        return False

def get_error_v3a(d):
    xlist=[d[i] for i in range(3,len(d), 2)]
    ylist=[d[i] for i in range(4,len(d), 2)]
    
    x_min=np.min(xlist)
    x_max=np.max(xlist)
    y_min=np.min(ylist)
    y_max=np.max(ylist)
    dx = (x_min-x_max)/2
    dy = (y_min-y_max)/2
    
    return (dx, dy)

def get_error_v3(d):
    xlist=d[0]
    ylist=d[1]
    
    x_min=np.min(xlist)
    x_max=np.max(xlist)
    y_min=np.min(ylist)
    y_max=np.max(ylist)
    
    dx = (x_min-x_max)/2
    dy = (y_min-y_max)/2
    
    return (dx, dy)


def fn_pseudosolution(x,y,fnpoly):  
    for xl in range(len(x)):        
        if xl == len(x)-1:
            fnpoly.write("%2.12f\t %2.12f\n" %(x[xl], y[xl]))
        else:
            fnpoly.write("%2.12f\t %2.12f\t" %(x[xl],y[xl]))
    return

def fn_realsolution(x,y,fcoor):
    for xl in range(len(x)):
        if xl == len(x)-1:            
            fcoor.write("%2.12f\t %2.12f\n"%(x[xl], y[xl]))
        else:
            fcoor.write("%2.12f\t %2.12f\t"%(x[xl], y[xl]))
    return


def fn_write(fn, data):
    
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


def fn_signcombination(r):
    f=[]
    
    for i in range(1, r+1):
        t = [-1]*i+[1]*(r-i)
        w = set(permutations(t))
        for u in w:
            f.append(u)
    return np.array(f)

def fn_repeat(p, d, f, imin, imax):
    
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


def fn_mesh(l, coordinates, imax):
    
    c = np.linspace(0,imax,int(2*l*imax+1) )
    
    k = [c, c]*len(coordinates)
    k = k[0:len(coordinates)]
    
    j = np.meshgrid(*k)
    
    [*dim] = np.shape(j)
    
    f1=(np.array([j[i].reshape(-1,1) for i in range([*dim][0])]))
    f2=np.hstack([f1[i] for i in range([*dim][0])])
    
    meshlist=np.array(f2)
    
    return meshlist




# def fn_solveforx_v2(l, g, f, m, j, x):
    
#     k = 2 * np.pi * l
    
#     i = list(range(j)) + list(range(j+1,len(x)))
    
#     a = (1-m*m)/(f[j]*f[j])
#     b = 2 * m*m * g /(f[j]*f[j])
#     c = m*m*( 1 - (g*g) / (f[j]*f[j]) ) - np.array([ (f[ii]*f[ii]) / (f[j]*f[j]) for ii in i]).sum(axis = 0) 
    
    
#     if a != 0:
        
#         z1 = (-b + np.sqrt(b*b - 4*a*c))/(2*a)
#         z2 = (-b - np.sqrt(b*b - 4*a*c))/(2*a)
        
#         r1=(1/k)*np.arccos(z1/f[0])
#         r2=(1/k)*np.arccos(z2/f[0])
        
#         if ~np.isnan(r1) and np.isnan(r2):
#             ry=findpy(r1,l,g,f)
#             return np.array([r1, ry])
        
#         elif np.isnan(r1) and ~np.isnan(r2):
#             ry=findpy(r2,l,g,f)
#             return np.array([r2, ry])
        
#         elif np.isnan(r1) and np.isnan(r2) :
            
#             return np.array([float("nan")])
        
#         elif ~np.isnan(r1) and ~np.isnan(r2):
            
#             r1y=findpy(r1,l,g,f)
#             r2y=findpy(r2,l,g,f)
            
#             print("Multiple solution: r1, r1y= ", r1,r1y, " r2, r2y = ", r2, r2y)
#             return np.array([ [r1, r1y], [r2, r2y] ])
        
#         else:
#             print("I obtained weird solution. I do not know what to return r1= ", r1," r2= ", r2)
#             return np.array([r1, r2])
    
#     else:
#         #print("===> printing a from else : ", a)
#         #print("===> printing b from else : ", b)
#         #print("===> printing c from else : ", c)
#         #print("===> printing m from else : ", m)
#         #print("===> printing g from else : ", g)
        
#         z1 = (-1*c/b)
        
#         r1=(1/k)*np.arccos(z1/f[0])
        
#         if ~np.isnan(r1):
#             ry=findpy(r1,l,g,f)
#             return np.array([r1, ry])
#         elif np.isnan(r1):
#             prx = (1/k)*np.arccos(g*(1+error)/np.sum(f))
#             pry = prx
#             return np.array([prx, pry])
#         else:
#             prx = (1/k)*np.arccos(g*(1+error)/np.sum(f))
#             pry = prx
#             return np.array([prx, pry])


####
# def jonaspnts_error(gi,l,error=0):
#     #### Jonas Area 
#     k    = 2*np.pi*l
#     gi   = np.abs(gi)
#     #### Finding point p1 and p2  
#     p1x  = (1/k)*np.arccos(gi*(1+error)/2) ;       #  x2* = x1*
#     p1y  = p1x
    
#     p2x  = (1/k)*np.arccos(gi*(1+error)-1)         #  x1* = pnt2 ; x2* = 0
#     p2y  = 0
    
#     m1   = (p2y-p1y)/(p2x-p1x)           # slope of First line
#     n1   = find_interception(p2x,p2y,m1)
        
#     #### Finding point p3, p4 and p5
#     xini = 0.0
#     p5x  = fsolve(findp3x,xini,args=(m1,gi*(1-error),l), factor=1, epsfcn=1e-16, maxfev=10000000)[0]
#     p5y  = findp3y(p5x,l,gi*(1-error))
    
#     n2   = find_interception(p5x,p5y,m1)
    
#     p4x  = -n2 / (m1-1)
#     p4y  = p4x
    
#     p3x  = -n2/m1
#     p3y  = 0
    
#     #pnt  = np.array([p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y, p5x, p5y])
#     pnt  = np.array([[p1x, p1y], [p2x, p2y], [p3x, p3y], [p4x, p4y],[ p5x, p5y]])
#     return pnt
