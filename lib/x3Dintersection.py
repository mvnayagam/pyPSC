import polytope as pc
import numpy as np

'''
This func is defined to find the intersection of successive / any pair of polytope regions
It is written in serial method. Parallelization could save time. However it is yet to be tested

future goal : make it in parallel code
'''
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

def find_intersection_updated(s, r):
    u=[]
    
    if type(s) is pc.Polytope:
        #print("from if: ", np.shape(s.A)[1])
        dim=np.shape(s.A)[1]
    elif type(s) is pc.Region:
        print("from elif: ", np.shape(s[0].A)[1])
        dim=np.shape(s[0].A)[1]
    else:
        print("type is not found")
    
    gp  = np.identity(dim)
    A  = np.array(np.vstack([-gp, gp]))
    b  = np.array([0]*dim + [0.5]*dim)
    PS = pc.Polytope(A, b)
    
    for count, i in enumerate(r):
        v = s & i
        
        if type(v) is pc.Polytope:
            if not pc.is_empty(v):
                u.append(PS & v)
        elif type(v) is pc.Region:
            for k in v:
                if not pc.is_empty(k):
                    u.append( PS & k)
    return pc.Region(u)