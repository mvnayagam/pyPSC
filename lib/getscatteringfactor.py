"""Returns Scattering Factor f=f0+f'+if'' of given Chemical Symbol

   Parameter: 	Diffraction Angle Theta in Degrees
				Chemical Symbol (Default = "Si")
				Energy of Incident X-Ray Beam in eV (Default Cu Ka1 = 8047.8227)
				
   adds 	f0 from "/atominfo/scatteringfactorcoeff.dat"
			f' from "/f1f2/chemsymb.dat"
			f''from "/f1f2/chemsymb.dat"			
"""

# from scipy.xplt import *									# modul for graphic output
#from IPython.Shell import IPythonShellEmbed				# modul for debug shell while program runs
#ipshell = IPythonShellEmbed()								# debug shell


import numpy as np

def getfs(l, d, E, symbs = 'H'):
    """
    l:
        l' scattering index (= 2*l for non-mirrorsymmetric structures)
    d:
        lattice constant [m]
    E:
        energy [eV]
    symbs:
        array of symbols, such as ['H', 'O']
    
    returns:
        array of scattering strengths, e.g. [fH = fH' + j*fH'', fO = fO' + j*fO'']
    """
    
    h, c, elec_charge = 6.626069e-34, 299792458, 1.60217646e-19
    angle = np.arcsin(l/d*h*c/(E*elec_charge))*180/np.pi
    fs = []
    
    for s in symbs:
        fs.append(getf(angle, s, energy = E)[3])
    
    return fs

def getf(theta, chemsymb="Si", energy=8047.8227):
	"""Returns Scattering Factor f=f0+f'+if'' of given Chemical Symbol

	   Parameter: 	Diffraction Angle Theta in Degrees
					Chemical Symbol (Default = "Si")
					Energy of Incident X-Ray Beam in eV (Default Cu Ka1 = 8047.8227)

	   adds 	f0 from "/atominfo/scatteringfactorcoeff.dat"
				f' from "/f1f2/chemsymb.dat"
				f''from "/f1f2/chemsymb.dat"
	"""
	
	from math import sin,exp,pi

	theta = float(theta)/360*2*pi		# Diffraction Angle
	energy = float(energy)						# Energy of beam
	datnamef1f2 = str(chemsymb)+".dat"			# Filename of f1f2 table
	
	
	if open("f1f2/"+datnamef1f2,"r"): 		# if table exists for specified Chem Symbol
		f = open("f1f2/"+datnamef1f2,"r")	# open file
		list = f.readlines()					# write lines to list
		f.close									# close file
		
		i = 4									# start with fifth line
		if energy < float(list[i].split()[0]):					
			print("Energy smaller than tabled range!")
		elif energy > float(list[len(list)-3].split()[0]):
			print("Energy higher than tabled range!")
		else:
			while float(list[i].split()[0]) < energy:		
				i = i+1
			ena, f1a, f2a = float(list[i-1].split()[0]),  float(list[i-1].split()[1]),  float(list[i-1].split()[2])
			enb, f1b, f2b = float(list[i].split()[0]),  float(list[i].split()[1]),  float(list[i].split()[2])
			frac = (energy-ena) / (enb-ena)                 	# linear interpolation for f' and f'' between neighboring energies
			f1 = f1a + frac*(f1b-f1a)
			f2 = f2a + frac*(f2b-f2a)
	else:
		print("Chemical Symbol is not listed in Table f' f'' !")
	
	f = open("atominfo/scatteringfactorcoeff.dat","r")		# table of coeff. for f0
	list = f.readlines()
	f.close
	i = 0
	
	for i in range(len(list)-3):
		if '"'+chemsymb+'",' in list[i]:						# get coeff
			a1, a2 = float(list[i].split()[3].strip(",")), float(list[i].split()[4].strip(","))
			a3, a4 = float(list[i].split()[5].strip(",")), float(list[i].split()[6].strip(","))
			b1, b2 = float(list[i+1].split()[1].strip(",")), float(list[i+1].split()[2].strip(","))
			b3, b4 = float(list[i+1].split()[3].strip(",")), float(list[i+1].split()[4].strip(","))
			c	   = float(list[i+2].split()[0].strip(","))
			
			# print a1,a2,a3,a4,b1,b2,b3,b4,c
			# raw_input()
			
			
			lamda = 2*pi*1973.271096 / energy	# wavelength in Angstroem
			q2 = (sin(theta)/lamda)**2							# argument of f0, squared in old convention IT (without 4*pi)
			f0 = a1 * exp(-b1*q2) + \
				 a2 * exp(-b2*q2) + \
				 a3 * exp(-b3*q2) + \
				 a4 * exp(-b4*q2) + \
				 c												# approx for f0 with retriefed coff according to intern. tables 6.1.1.4
			break
		if i == len(list)-3:
			print("Chemical Symbol is not listed in Table f0!")
			break
		i = i+1
	return f0,f1,f2,complex(f0+f1,f2)
	

# test program for routine given above

# res = 200											# set resolution of calculation
# angles = 1.0*arange(res)*90/float(res)				# choose range for diffraction angles
# atom = "Si"											# choose Sort of Atoms
# energy = 8047.8227									# choose energy of incident beam

# fs = zeros((res,4),dtype=complex)					# initialize array for scatteringfactors
# i=0													# indicator numeric calculations	
# for i in range(res):
	# fs[i,:] = getf(angles[i],atom,energy)			# call above routine for angles, sort of atom, energy
	# i=i+1
	# #print i

# window(0)                                  					# plotproperties
# title("Atomic Scattering Factor of "+atom,fontsize=24) 
# xlabel("theta"); ylabel("f")                        
# pldefault(marks=0)                         
# limits(min(angles),max(angles),min(min(fs[:,1].real),min(fs[:,2].real)-2),max(fs[:,0].real)+5) 

# plg(fs[:,0].real,angles,color="black",width=4)             	# plot f0
# plg(fs[:,1].real,angles,color="red")                       	# plot f'
# plg(fs[:,2].real,angles,color="blue")                      	# plot f'' 
# plt("f0",0.2,0.77,color="black")      
# plt("f'",0.2,0.81,color="red")            
# plt("f''",0.2,0.79,color="blue")    
         
# print "EXIT to leave imbedded shell"
# ipshell()
# raw_input()
# winkill(0)
