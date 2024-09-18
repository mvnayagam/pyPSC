# ----------------------------------------------------------
# Dt: 18.09.2024    by: Muthu
# importing all functions for general usage
# ----------------------------------------------------------
from .createfolder import createmcdir

# ----------------------------------------------------------
# Call some tools from lib for structure solving 
# ----------------------------------------------------------
from .MCin2DPS_EPA  import MC2DPS_EPA
from .MCin2DPS_nEPA import MC2DPS_nEPA
from .MCinNDPS_EPA  import isosurfs_EPA, MCNDPS_EPA

# ----------------------------------------------------------
# Main library for psc structure solving
# ----------------------------------------------------------
from .g_space import g, F, grad_F, grad_g, hsurf_F, hsurf_F2, hsurf_g

from .x2Dlinearization import find_interception, findpx, findpy, fn_solveforx_v2
from .x2Dlinearization import double_segment_EPA, single_segment_EPA, single_segment_nEPA, single_segment_nEPA, double_segment_nEPA
from .x2Dpolygon       import multistrip, getploygons_EPA_SS, getploygons_EPA_DS, getploygons_nEPA, polyintersect, polyintersect_MC
from .x2Drepetition    import repeat2D, linrep_DS, linrep_SS, writedata
from .x2Dwritesolution import writepolygons, isInside, get_error, get_error_v3a, pseudosolution, realsolution, analyzesolution

from .x3Dchecklinearization import checklinear, checklinear_I, checklinearplot, getpoly_mitd
from .x3Dintersection       import find_intersection
from .x3DlinearizationEPA   import linearizenD_EPA, linearizenD_nEPA
from .x3Dreadwrite          import wrtdata, wrtcoor, wrtvolume, wrtallsolution, readoldsolution, readh5file, readh5file_v2
from .x3Drepetition         import getmesh, getsigncombination, getpolytope, getpolytope_EPA, getpolytope_nEPA

from .xlinearizationtools   import radius_from_volume
