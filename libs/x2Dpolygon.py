import numpy as np
from shapely.geometry import Polygon, mapping, Point

def multistrip(a, b, pnts):
    d=[]
    for i in range(a,b,1):
        df=int(len(pnts[i]))
        d1=Polygon([(pnts[i][j], pnts[i][j+1]) for j in range(0,df,2)])
        d.append(d1)
    return shapely.geometry.MultiPolygon([poly for poly in d])

def getploygons(h, points, imax=0.5):
    
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

def getintersections(h,a,xexp,fname,count):
    
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

def writepolygons(fname, polys):
    
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

def get_error_v2a(d):
    xlist=[d[i] for i in range(3,len(d), 2)]
    ylist=[d[i] for i in range(4,len(d), 2)]
    
    x_min=np.min(xlist)
    #x_min_inx=np.where(xlist == x_min)[0]
    #y_min=ylist[x_min_inx[0]]
    
    x_max=np.max(xlist)
    #x_max_inx=np.where(xlist == x_max)[0]
    #y_max=ylist[x_max_inx[0]]
    
    y_min=np.min(ylist)
    y_max=np.max(ylist)
    
    dx = np.abs(x_min-x_max)/2
    dy = np.abs(y_min-y_max)/2
    
    return (dx, dy)

def fn_get_error_v2(d):
    xlist=d[0]
    ylist=d[1]
    
    x_min=np.min(xlist)
    x_max=np.max(xlist)
    
    #x_min_inx=np.where(xlist == x_min)[0]
    #y_min=ylist[x_min_inx[0]]
    #x_max_inx=np.where(xlist == x_max)[0]
    #y_max=ylist[x_max_inx[0]]
    
    y_min=np.min(ylist)
    y_max=np.max(ylist)
    
    dx = np.abs(x_min-x_max)/2
    dy = np.abs(y_min-y_max)/2
    
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