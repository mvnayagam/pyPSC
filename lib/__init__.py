#   --- importing all functions for general usage

from .g_space import g, F, grad_F, grad_g, hsurf_F, hsurf_F2, hsurf_g

from .x3Dlinearizazion import linearizenD_EPA, linearizenD_nEPA
from .x3Drepetition import getmesh, getsigncombination, getpolytope, getpolytope_EPA, getpolytope_nEPA

from .x3Dchecklinearization import checklinear, checklinear_I, checklinearplot, getpoly_mitd
from .x3Dintersection import find_intersection
from .x3Dreadwrite import wrtdata, wrtcoor, wrtvolume, wrtallsolution, readoldsolution, readh5file, readh5file_v2


from .x2Dlinearize import find_interception, findpx, findpy, fn_solveforx_v2
from .x2Dlinearize import double_segment_EPA, single_segment_EPA, single_segment_nEPA, single_segment_nEPA, double_segment_nEPA
from .x2Dpolygon import multistrip, getploygons_EPA_SS, getploygons_EPA_DS, getploygons_nEPA, polyintersect, polyintersect_MC

from .x2Drepetition import repeat2D, linrep_DS, linrep_SS, writedata
from .x2Dwritesolution import writepolygons, isInside, get_error, get_error_v3a, pseudosolution, realsolution, analyzesolution


# ----------------------------------------------------------
# Call some tools even in lib for structure solving 
# ----------------------------------------------------------

from .solvestructureNDPS import solvenDstrucutre_EPA, solvenDstrucutre_nEPA
from .solvestructure2DPS import solve2Dstrucutre_nEPA, solve2Dstrucutre_EPA