import numpy as np

def intersect_list(f, g):
    """
    x:
        x-values (in order to make sure that the same x-axis are used)
    f,g:
        two functions (represented by discrete points on the same grid x)
    
    returns:
        ni:
            array indices
    """ 
    #x:
    #    x-values (in order to make sure that the same x-axis are used)
    
    d  = f-g
    ni = np.argwhere((d[1:]*d[:-1]) < 0).flatten() # quick and dirty
    return np.array(ni) 

def intersect_interp(x, f, g): 
    """
    f,g:
        two functions (represented by discrete points on the same grid x)
    
    returns:
            xi:
                x-values
    """
    inds = intersect_list(f, g)
    
    xi, yi = [], []
    for i in iter(inds):
        imat = np.linalg.inv([x[i-1:i+2]**2, x[i-1:i+2], np.ones_like(x[i-1:i+2])])
        a, b, c = np.dot(imat, f[i-1:i+2]-g[i-1:i+2]) 
        #a, b, c = np.polyfit(x[i-1:i+2], f[i-1:i+2]-g[i-1:i+2], deg = 2) # squared interpolation

        x12 = [-b/(2*a)+np.sqrt(b*b/(4*a*a)-c/a), -b/(2*a)-np.sqrt(b*b/(4*a*a)-c/a)] # zeros of squared interpolation
        # print(a, b, c)
        # select the reasonable zero -> x
        if (x12[1] > x[i-1]) and (x12[1] < x[i+1]):
            xi.append(x12[1])
        else:
            xi.append(x12[0])
            
    return np.array(xi)

def angle(x, y):
    """
    intersection angle between two n-dim. vectors.
    x, y:
        array-like objects (will internally be casted to arrays)
    """
    return np.arccos(np.dot(x, y)/(np.linalg.norm(x)*np.linalg.norm(y)))

#def angle(x, y):
#    return np.arccos(np.dot(x, y)/(np.linalg.norm(x)*np.linalg.norm(y)))