import numpy as np
import intvalpy as ip


# -------------------- EPA linearization

def find_interception(x,y,m):
    return y-m*x

def findp3x(x, m, g, h):
    return m+np.sin(2*np.pi*h*x)/np.sqrt(1- (g - np.cos(2*np.pi*h*x))**2) 

def findp3y(x,h,g):
    return 1/(2*np.pi*h)*np.arccos(g-np.cos(2*np.pi*h*x))

def  double_segment_EPA(gi, l, xexp, f, error=0):
    
    #### Jonas Area 
    k    = 2*np.pi*l ;  gi   = np.abs(gi)
    
    #### Finding point p1 and p2  
    p1x  = (1/k)*np.arccos(gi*(1+error)/np.sum(f)) ;       #  x2* = x1*
    p1y  = p1x
    
    p2x  = (1/k)*np.arccos(gi*(1+error)-1)         #  x1* = pnt2 ; x2* = 0
    p2y  = 0
    
    m1   = (p2y-p1y)/(p2x-p1x)           # slope of First line
    n1   = find_interception(p2x,p2y,m1)
        
    #### Finding point p3, p4 and p5
    
    j = 1
    p5x, p5y = fn_solveforx_v2(l, gi, f, m1, j, xexp)
    
    n2   = find_interception(p5x,p5y,m1)
    
    p4x  = -n2 / (m1-1)
    p4y  = p4x
    
    p3x  = -n2/m1
    p3y  = 0
    
    #pnt  = np.array([p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y, p5x, p5y])
    pnt  = np.array([[p1x, p1y], [p2x, p2y], [p3x, p3y], [p5x, p5y], [p4x, p4y]])
    
    return pnt

def single_segment_EPA(gi, l, xexp, f, error=0):
    #### Jonas Area 
    k    = 2*np.pi*l
    gi   = np.abs(gi)
    
    #### Finding point p1 and p2  
    
    p1x  = (1/k)*np.arccos(gi*(1+error)-1)         #  x1* = pnt2 ; x2* = 0
    p1y  = 0
    
    p2y  = (1/k)*np.arccos(gi*(1+error)-1)         #  x1* = pnt2 ; x2* = 0
    p2x  = 0
    
    m1   = (p2y-p1y)/(p2x-p1x)           # slope of First line
    n1   = find_interception(p2x,p2y,m1)
        
    #### Finding point p3, p4 and p5
    
    p5x = (1/k)*np.arccos(gi*(1-error)/np.sum(f))
    p5y = p5x #p4  = [xp]*len(f)
    
    n2   = find_interception(p5x,p5y,m1)
    
    p4y  = n2 #/ (m1-1)  # y = m1x+n2
    p4x  = 0 #p4x
    
    p3x  = -n2/m1
    p3y  = 0
    
    pnt  = np.array([[p2x, p2y], [p1x, p1y], [p3x, p3y], [p5x, p5y], [p4x, p4y]])
    
    return pnt


# -------------------- nEPA linearization

def single_segment_nEPA(gi, l, f, xcoor, j=1, error=0):
    pnt = []
    
    k   = 2*np.pi*l
    gi   = np.abs(gi)
      
    p1x  = 0
    p1y  = (1/k)*np.arccos((gi*(1+error)-f[0])/f[1])
    
    if np.isnan(p1y):
        p1y = 0.5/l
        p1x = findpx(p1y, l, gi, f)
        pnt.append([p1x, p1y])
        
    else:
        pnt.append([p1x, p1y])        
     
    
    p2x  = (1/k)*np.arccos((gi*(1+error)-f[1])/f[0])
    p2y  = 0
    pnt.append([p2x, p2y]) 
    
    
    m1   = (p2y-p1y)/(p2x-p1x)
    n1   = find_interception(p2x,p2y,m1)
    
    p4   = fn_solveforx_v2(l, gi, f, m1, j, xcoor)
            
    if ~np.all(np.isnan(p4)):
        
        if len(np.shape(p4)) > 1:
            pnt = []
            
            if np.floor(p4[0,1]/p4[0,0]) <= 1 and np.floor(p4[1,1]/p4[1,0]) > 1:
                pLB=p4[0]
                pUB=p4[1]
            else:
                pLB=p4[1]
                pUB=p4[0]
            
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
        
        elif ~np.all(np.isnan(p4)):
            p4y = p4[1]
            p4x = p4[0]
        else:
            
            p4y = 0.5/l
            p4x  = (1/k)*np.arccos( gi/f[0] - (f[1]/f[0])*np.cos(k*p4y) )
    
    n2   = find_interception(p4x,p4y,m1)
    p3x  = -n2/m1
    p3y  = 0
    
    pnt.append([p3x, p3y])
    pnt.append([p4x, p4y])
    
    p5x  = 0
    p5y  = n2
    
    if p5y > 0.5/l:
        p5y  = 0.5/l
        p5x = (p5y-n2)/m1
        pnt.append([p5x, p5y])
        
        pextra_x = 0
        pextra_y = 0.5/l
        pnt.append([pextra_x, pextra_y])
    else:
        pnt.append([p5x, p5y])
        
    return np.array(pnt)

def double_segment_nEPA(gi, l, f, xcoor, j=1, error=0):
    pnt = []
    k   = 2*np.pi*l
    gi   = np.abs(gi) 
    
    p1x = (1/k)*np.arccos(gi*(1+error)/np.sum(f))
    p1y = p1x
    pnt.append([p1x, p1y])
    
    p2x = (1/k)*np.arccos((gi*(1+error)-f[1])/f[0])
    p2y = 0
    
    if np.isnan(p2x):
        p2x = 0.5/l
    
    pnt.append([p2x, p2y])
    
    m1  = (p2y-p1y)/(p2x-p1x)
    n1  = find_interception(p2x,p2y,m1)
    
    p5   = fn_solveforx_v2(l, gi, f, m1, j, xcoor)
    
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
            p5x  = (1/k)*np.arccos( gi/f[0] - (f[1]/f[0])*np.cos(k*p5y) )
    
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
    p6y  = (1/k)*np.arccos((gi*(1+error)-f[0])/f[1])
    
    if np.isnan(p6y):
        p6y = 0.5/l
        
    m3   = (-p6y+p1y)/(-p6x+p1x)
    
    p7 = fn_solveforx_v2(l, gi, f, m3, j, xcoor)
    
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
        p7x  = (1/k)*np.arccos( gi/f[0] - (f[1]/f[0])*np.cos(k*p7y) )  # === p7x  = findpx(p7y,l,g,f)
        
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

def findpy(x,h,gi,f):
    return (1/(2*np.pi*h))*np.arccos(gi/f[1] - (f[0]/f[1])*np.cos(2*np.pi*h*x))

def findpx(y,h,gi,f):
    return (1/(2*np.pi*h))*np.arccos(gi/f[0] - (f[1]/f[0])*np.cos(2*np.pi*h*y))

def fn_solveforx_v2(l, gi, f, m, j, x):
    
    k = 2 * np.pi * l
    
    i = list(range(j)) + list(range(j+1,len(x)))
    
    a = (1-m*m)/(f[j]*f[j])
    b = 2 * m*m * gi /(f[j]*f[j])
    c = m*m*( 1 - (gi*gi) / (f[j]*f[j]) ) - np.array([ (f[ii]*f[ii]) / (f[j]*f[j]) for ii in i]).sum(axis = 0) 
    
    
    if a != 0:
        
        z1 = (-b + np.sqrt(b*b - 4*a*c))/(2*a)
        z2 = (-b - np.sqrt(b*b - 4*a*c))/(2*a)
        
        r1=(1/k)*np.arccos(z1/f[0])
        r2=(1/k)*np.arccos(z2/f[0])
        
        if ~np.isnan(r1) and np.isnan(r2):
            ry=findpy(r1, l, gi, f)
            return np.array([r1, ry])
        
        elif np.isnan(r1) and ~np.isnan(r2):
            ry=findpy(r2, l, gi, f)
            return np.array([r2, ry])
        
        elif np.isnan(r1) and np.isnan(r2) :
            return np.array([float("nan")])
        
        elif ~np.isnan(r1) and ~np.isnan(r2):
            r1y=findpy(r1, l, gi, f)
            r2y=findpy(r2, l, gi, f)
            return np.array([ [r1, r1y], [r2, r2y] ])
        
        else:
            return np.array([r1, r2])
    
    else:
        z1 = (-1*c/b)
        
        r1=(1/k)*np.arccos(z1/f[0])
        
        if ~np.isnan(r1):
            ry=findpy(r1, l, gi, f)
            return np.array([r1, ry])
        elif np.isnan(r1):
            prx = (1/k)*np.arccos(gi*(1+error)/np.sum(f))
            pry = prx
            return np.array([prx, pry])
        else:
            prx = (1/k)*np.arccos(gi*(1+error)/np.sum(f))
            pry = prx
            return np.array([prx, pry])
