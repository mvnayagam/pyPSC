import polytope as pc


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
