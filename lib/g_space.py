import numpy as np
import sys

"""
naming convenctions: 
g:    scattering amplitude with fixed scattering strengths f
F:    scattering amplitude with atom-dependent scattering strengths f
F2:   sloppy form of F.F* = I (intensity)
_fixed:  including fixed scatterers beyond the variable structure
hsurf    hypersurface
grad:    gradient at a given position
"""

def g(h: int, x: list, f: list) -> float:
    """returns scattered amplitude with fixed scattering strength 
       f for given atomic positions x. Centrosymmetric case.
    Args:
        h: reflection index
        f: atomic scattering factor, scalar quantity
        x: array of atom positions
    Returns:
        float: scattered amplitude
    """
    
    if len(x) != len(f):
        print("\x1b[1;31;43m===> len(x) and len(f) are not same. I am exciting")
        sys.exit()

    return ( sum([f[i]*np.cos(2*np.pi*h*x[i]) for i in range(len(x))] ))

def F(h: int, x: list, f: list) -> float:
    """returns scattered amplitude with individual scattering strengths 
       f(i) for given atomic positions x(i), i - atom index. Centrosymmetric case.
       
    Args:
        h (int): array of possible atomic coordinates
        x (list): array of atom positions
        f (list): array of scattering factors. simply set to 1.
    Returns:
        float: scattered amplitude
    """
    
    return (sum([f[i]*np.cos(2*np.pi*h*x[i]) for i in range(len(f))]))

def hsurf_g(h: int, x: list, f: list, gi: float, j: int, s: int=1) -> list:
    """returns hyper isosurface for given amplitude gi
    
    Args:
        h (int): array of possible atomic coordinates
        x (list): array of atom positions
        f (list): array of scattering factors. simply set to 1.
        gi (float): the amplitude of intensity
        j (int): for which atoms the hsurf is to be solved
        s (int, optional): sign of amplitude. always +/- 1

    Returns:
        list: returns hyper isosurface for given amplitude gi
    """
    
    ilist = list(range(j)) + list(range(j+1,len(x))) # select parameter indices except index j
    
    k = 2*np.pi*h
    argm = s*gi/f[j] - np.array([(f[i]/f[j])*np.cos(k*x[i]) for i in ilist]).sum(axis = 0)
    xj = (np.arccos(argm))/k
    
    return xj

def hsurf_F(h: int, x: list, f: list, I: float, j: int, s: int =1, s2: int = 1) -> list:
    """ list of possible coordinates along j direction
    
    Args:
        h (int): reflection index
        x (list): array of atom positions
        f (list): array of atomic scattering factors
        I (float): intensity
        j (int): selected coordinate
        s (int, optional): sign of amplitude (1: positive, -1: negative)
        s2 (int, optional): sign of ordinate (1: positive, -1: negative)
        nan (str, optional): Defaults to True.

    Returns:
        list: possible coordinates along j direction
    """
    
    s, s2  = np.sign(s),np.sign(s2) 
    inds = list(range(j)) + list(range(j+1,len(f))) # select parameter indices except index j
    #same as inds = np.delete(np.arange(len(f)), j)
    
    gj = [f[i]*np.cos(2*np.pi*h*x[i]) for i in inds]
    
    return s2*np.arccos((s*np.sqrt(I) - sum(gj))/f[j])/(2*np.pi*h)

def hsurf_F2(h: int, x: list, f: list, I: float, j: int, s: int =1, s2: int = 1, nan: str = True) -> list:
    """_summary_

    Args:
        h (int): reflection index
        x (list): array of atom positions
        f (list): array of atomic scattering factors
        I (float): intensity
        j (int): selected coordinate
        s (int, optional): sign of amplitude (1: positive, -1: negative)
        s2 (int, optional): sign of ordinate (1: positive, -1: negative)
        nan (str, optional): _description_. Defaults to True.

    Returns:
        list: possible coordinates along j direction
    """
    
    s, s2  = np.sign([s, s2])
    inds = list(range(j)) + list(range(j+1,len(f))) # select parameter indices except index j
    #same as inds = np.delete(np.arange(len(f)), j)
    
    # double sum over all indices except index j
    gj = [f[i]*np.conj(f[k])*np.cos(2*np.pi*h*x[i])*np.cos(2*np.pi*h*x[k]) for i in inds for k in inds]
    cj = np.array(gj).sum(axis = 0) - I
    
    # single sum over all indices except index j
    b1j = np.array([f[i]*np.cos(2*np.pi*h*x[i]) for i in inds]).sum(axis = 0)*np.conj(f[j])
    b2j = np.array([np.conj(f[i])*np.cos(2*np.pi*h*x[i]) for i in inds]).sum(axis = 0)*f[j]
    bj = b1j + b2j
    
    # contribution from index j
    aj = f[j]*np.conj(f[j])
    
    # print(cj, bj, aj)
    zj = (-bj + s*np.sqrt(bj**2 - 4*aj*cj))/(2*aj)
    xj = np.real(s2*np.arccos(zj)/(2*np.pi*h))

    if nan:
        ifalse = np.where((bj**2 < 4*aj*cj) + (np.abs(zj) > 1.0))
        xj[ifalse] = np.nan

    return xj

def grad_g(h: int, x: list, f:list) -> list:
    """Returns the gradient vector at point x
    
    Args:
        h (int): represent reflection
        x (list): test structure 
        f (list/array): list of atomic scattering factors
        normalize (bool, optional): Defaults to True. This normalizes gradient vector
        
    Returns:
        _array_: Returns the gradient vector as array
    """    
    k = 2*np.pi*h
    
    if len(x) > 1:
        gradG = np.array([-f[i]*k*np.sin(k*x[i]) for i in range(len(f))])
        return gradG/np.linalg.norm(gradG)
    else:
        gradG = np.array([-f[i]*k*np.sin(k*x) for i in range(len(f))])
        return gradG/np.linalg.norm(gradG)

def grad_F(h: int, x: list, f:list, normalize=True) -> list:
    """Returns the gradient vector at point x

    Args:
        h (int): represent reflection
        x (list): test structure 
        f (list/array): list of atomic scattering factors
        normalize (bool, optional): Defaults to True. This normalizes gradient vector
        
    Returns:
        _array_: Returns the gradient vector as array
    """      
    k = 2*np.pi*h
    
    gradF = np.array( [-k*f[i]*np.sin(k*x[i]) for i in range(len(f))])
    
    return gradF/np.linalg.norm(gradF)