
#! /usr/bin/env python

################################################################################
# Ruth Lunt 2013                                                               #
################################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pylab import *
from optparse import OptionParser


######################## Set up optional arguments #############################

# specify input file
parser = OptionParser()
parser.add_option("-f", "--file",
                  action = "store", type = "string", dest = "file", default = "GULP.gin",
                  help = "Path to input file [default: ./geometry.in]")
(options, args) = parser.parse_args()
file = options.file


########################### Begin main program #################################

print "A program to calculate the effective spring constants, C6 Lennard-Jones parameter and bsm values for oxygen in a mixed metal oxide. As well as plotting the potentials for the metal oxide."
print "Ruth Lunt 2013\nDate last edited: 20/08/2013"

# Open the file and split into Lines
f = open(file,"r")
lines = f.readlines()
f.close()

# lists to be appended
buck = []
coul = []
lenn_fix = []
lenn_mix = []
spring_fix = []
spring_mix = []
bsm_mix = []
p_fix = []
copy_lenn_fix = 0
copy_lenn_mix = 0
copy_spring_fix = 0
copy_spring_mix = 0
copy_bsm = 0
copy_p = 0

# determination of metal oxide composition
composition = raw_input("What is the metal composition of the mixed metal oxide?")
elements = composition.split()

# start reading lines
for element in elements:
	for line in lines:
		 inp = line.split()
		 if inp == []:
		  continue
		 if len(inp) == 12 and inp[0] == element:
			buck.append(inp)
		 
		 if len(inp) == 3 and inp[0] == element and (inp[1] == "core" or inp[1] == "shel"):
			coul.append(inp)
		
		 if copy_spring_fix == 1 and inp[0] == element:
			spring_fix.append(inp)
			copy_spring_fix = 0
		 if len(inp) == 2 and (inp[0] == "spring" and inp[1] == element):
			copy_spring_fix = 1	
		
		 if copy_spring_mix == 1 and inp[0] == "O":
			copy_spring_mix = 0
			# addition of 'element label' allows for dinstinction between multiple oxygen spring constants
			inp = element.split() + inp
			spring_mix.append(inp)
		 if len(inp) ==2 and (inp[0] == "spring" and inp[1] == element):
			copy_spring_mix = 1
		
		 # copy_lenn_fix differnet as needs to read in multiple lines for lennard 
		 if copy_lenn_fix <= 2 and copy_lenn_fix > 0 and (inp[0] == element or inp[2] == element):
			copy_lenn_fix = copy_lenn_fix + 1
			lenn_fix.append(inp)
		 if len(inp) == 2 and (inp[0] == "lennard" and inp[1] == element):
			copy_lenn_fix = 1
			
		 if copy_lenn_mix == 1 and inp[0] == "O" and inp[2] == "O":
			copy_lenn_mix = 0
			# addition of 'element label' allows for dinstinction between multiple oxygen C6 values
			inp = element.split() + inp
			lenn_mix.append(inp)
		 if len(inp) == 2 and (inp[0] == 'lennard' and inp[1] == element):
			copy_lenn_mix = 1

	 	 if copy_bsm == 1 and inp[0] == "O":
			copy_bsm = 0
			# addition of 'element label' allows for dinstinction between multiple oxygen bsm values
			inp = element.split() + inp
			bsm_mix.append(inp) 
		 if len(inp) == 2 and (inp[0] == "bsm" and inp[1] == element):
			copy_bsm = 1
			
		 if copy_p == 1 and len(inp) == 1:
			copy_p = 0
			# allows for distinction between mulitple oxygen P values (electron number - effective number of electrons contributing to polrizability) 
			# these values are obtained from the static dipole of polarisability (alpha) and C6 values for the pecific metal oxides
			inp = element.split() + inp
			p_fix.append(inp)
		 if len(inp) == 2 and (inp[0] == "P" and inp[1] == element):
			copy_p = 1
			
					
# getting the k2 values for oxygen
for spring in spring_mix:
	if spring[0] == "Zn":
		k2_Zn = spring[2]
	
	elif spring[0] == "Sn":
		k2_Sn = spring[2]
		
	elif spring[0] == "In":
		k2_In = spring[2]
		
# getting the k4 values for oxygen
for spring in spring_mix:
	if spring[0] == "Zn":
		k4_Zn = spring[3]
		
	elif spring[0] == "Sn":
		k4_Sn = spring[3]
		
	elif spring[0] == "In":
		k4_In = spring[3]
		
# getting bsm values for oxygen
for bsm in bsm_mix:
	if bsm[0] == "Zn":
		bsm_Zn = bsm[3]
		
	elif bsm[0] == "Sn":
		bsm_Sn = bsm[3]
		
	elif bsm[0] == "In":
		bsm_In = bsm[3]
		
# getting p values for oxygen
for p in p_fix:
	if p[0] == "Zn":
		P_Zn = p[1]
		
	elif p[0] == "Sn":
		P_Sn = p[1]
		
	elif p[0] == "In":
		P_In = p[1]
		

# determination of metal oxide formula
formula = raw_input("What is the elemental ratio of the mixed metal oxide (same order as composition)?")
ratio = formula.split()

#values for interpoltion formula
x = 0
y = 0
z = 0

# on addition of other elements into the code this section must be modified
if elements[0] == "Zn":
	x = ratio[0]
	if elements[1] == "Sn":
		z = ratio[1]
	elif elements[2] == "Sn":
		z = ratio[2]
	if elements[1] == "In":
		y = ratio[1]
	elif elements [2] == "In":
		y = ratio[2]
	
elif elements[0] == "In":
	y = ratio[0]
	if elements[1] == "Zn":
		x = ratio[1]
	elif elements[2] == "Zn":
		x = ratio[2]
	if elements[1] == "Sn":
		z = ratio[1]
	elif elements [2] == "Sn":
		z = ratio[2]

elif elements[0] == "Sn":
	z = ratio[0]
	if elements[1] == "Zn":
		x = ratio[1]
	elif elements[2] == "Zn":
		 x = ratio[2]
	if elements[1] == "In":
		y = ratio[1]
	elif elements [2] == "In":
		y = ratio[2]
				
# to avoid zero point errors if a metal is not present in the oxide the k2, k4, C6, P and bsm values for that metal oxide are set to 1, this will still give a zero value for that term within the interpolation formula
if x == 0:
	k2_Zn = 1
if y == 0:
	k2_In = 1
if z == 0:
	k2_Sn = 1
	
if x == 0:
	k4_Zn = 1
if y == 0:
	k4_In = 1
if z == 0:
	k4_Sn = 1
	
if x == 0:
	bsm_Zn = 1
if y == 0:
	bsm_In = 1
if z == 0:
	bsm_Sn = 1
	
if x == 0:
	P_Zn = 1
if y == 0:
	P_In = 1
if z == 0:
	P_Sn = 1
	

# calculatuon of effective oxygen k2 spring constant
def spring_k2(x, y, z, k2_Zn, k2_Sn, k2_In):
	k2_eff = (float(x) + float(y) + float(z))/((float(x)/float(k2_Zn)) + (float(y)/float(k2_In)) + (float(z)/float(k2_Sn)))
	return k2_eff

k2 = spring_k2(x, y, z, k2_Zn, k2_Sn, k2_In)
	
# calculation of effective oxygen k4 spring constant
def spring_k4(x, y, z, k4_Zn, k4_Sn, k4_In):
	k4_eff = (float(x) + float(y) + float(z))/((float(x)/float(k4_Zn)) + (float(y)/float(k4_In)) + (float(z)/float(k4_Sn)))
	return k4_eff

k4 = spring_k2(x, y, z, k4_Zn, k4_Sn, k4_In)
	
# calculation of effective oxygen bsm
def bsm(x, y, z, bsm_Zn, bsm_Sn, bsm_In):
	bsm_eff = (float(x) + float(y) + float(z))/((float(x)/float(bsm_Zn)) + (float(y)/float(bsm_In)) + (float(z)/float(bsm_Sn)))
	return bsm_eff

bsm = bsm(x, y, z, bsm_Zn, bsm_Sn, bsm_In)

print "effective oxygen k2 = " + "{0:6.4f}".format(k2)
print "effective oxygen k4 = " + "{0:6.4f}".format(k4)
print "effective oxygen bsm = " + "{0:6.4f}".format(bsm)

# oxygen shell charge
for i in coul:
	if i[0] == "O" and i[1] == "shel":
		Y = float(i[2])

# calculation of alpha (static dipole of polarizability) for metal oxide (to later calculate C6)
def alpha_eff(Y, k2):
	alpha = (Y**2)/k2
	return alpha

alpha = alpha_eff(Y, k2)

print "effective oxygen alpha = " + "{0:6.4f}".format(alpha)

# calculation of P for metal oxide
def P_eff(x, y, z, P_Zn, P_Sn, P_In):
	P = (float(x) + float(y) + float(z))/((float(x)/float(P_Zn)) + (float(y)/float(P_In)) + (float(z)/float(P_Sn)))
	return P

P = P_eff(x, y, z, P_Zn, P_Sn, P_In)
	
print "effective oxygen P = " + "{0:6.4f}".format(P)

# calculation of C6 for mixed metal oxide
def C6_eff(alpha, P):
		C6 = (0.75) * (alpha**(1.5) * (P**0.5))
		return C6 

C6 = C6_eff(alpha, P)
	
print "effective oxygen C6 = " + "{0:6.4f}".format(C6)


# addition of the effective oxygen k2 and k4 values to the spring_fix list
# only need one O spring term - can remove the rest from spring_mix list
spring_mix = spring_mix[0]
spring_mix[2] = k2
spring_mix[3] = k4
# to remove 'element label' from list
spring_mix = spring_mix[1:6]
spring_fix.append(spring_mix)

#Addition of effective oxygen bsm value
bsm_mix = bsm_mix[0]
bsm_mix[3] = bsm
bsm_mix = bsm_mix[1:7]

# Addition of oxygen C6 value to lenn_fix list
lenn_mix = lenn_mix[0] 
lenn_mix[6] = C6
lenn_mix = lenn_mix[1:11]
lenn_fix.append(lenn_mix)


# Buckingham potential - for plot
def buck_pot(a, rho, c, cut):
	array_size = int(cut/0.01)
	buck = np.zeros(shape= (array_size,2))
	for x in range (1, array_size):
		buck[x, 0] = x*0.01
		buck[x, 1] = a*np.exp(-buck[x, 0]/rho) - c**6/buck[x, 0]**6
	return buck
		
# coulomb potential - for plot
def coulomb_pot(q1, q2, cut):
	array_size = int(cut/0.01)
	coulomb = np.zeros(shape = (array_size, 2))
	for x in range (1, array_size):
		coulomb[x, 0] = x*0.01
		#covert Hartree to eV and Amstrom to Bohr
		coulomb[x, 1] = 27.211396132*(q1*q2)/(x*0.01*1.889725989)
	return coulomb		

# Lennard-Jones potential - for plot
def lennard_pot(A, B, cut):
	array_size = int(cut/0.01)
	lenn = np.zeros(shape = (array_size, 2))
	for x in range (1, array_size):
		lenn[x, 0] = x*0.01
		lenn[x, 1] = A/((x*0.01)**12) - B/((x*0.01)**6)
	return lenn		


# print out potentials
print "species"
coul_input = '\n'.join(str(i) for i in coul)
print coul_input

print "buck"
buckingham = '\n'.join(str(i) for i in buck)
print buckingham

print "lennard"
lennard = '\n'.join(str(i) for i in lenn_fix)
print lennard

print "spring"
spring = '\n'.join(str(i) for i in spring_fix)
print spring

print "bsm"
print bsm_mix


# Lennard-Jones paramters set to incase there are no values in input file
A = 0
B = 0
# Getting values to plot potentials
for i in buck:
	buck = buck_pot(float(i[4]), float(i[5]), float(i[6]), float(i[8])) 
	for j in coul:
		if j[0] == i[0] and j[1] == i[1]:
			q1 = j[2]
		if j[0] == i[2] and j[1] == i[3]:
			q2 = j[2]
	coulomb = coulomb_pot(float(q1), float(q2), float(i[8]))
	for k in lenn_fix:
		if k[0] == i[0] and k[2] == i[2]:
			A = k[4]
			B = k[5]
		# in case element order is reversed in input file
		elif k[2] == i[0] and k[0] == i[2]:
			A = k[4]
			B = k[5]
        lenn_j = lennard_pot(float(A), float(B), float(i[8]))
        #total_pot = np.add(buck,coul,lenn_j)
	total_pot = buck + coulomb + lenn_j
	
	plt.plot(buck[1:, 0],buck[1:, 1], label = 'Buckingham potential')
	plt.plot(coulomb[1:, 0], coulomb[1:, 1], label = 'Coulombic interaction')
	plt.plot(lenn_j[1:, 0], lenn_j[1:, 1], label = 'Lennard potential')
	plt.plot(buck[1:, 0], total_pot[1:, 1], label = 'Total potential')
	plt.legend(('Buckingham potential', 'Coulombic interaction', 'Lennard potential', 'Total potential'))
	plt.axis([0.00, 5.00, -150, 150])
	plt.xlabel('Interatomic distance, r', fontsize = 16)
	plt.ylabel('Potential Energy, eV', fontsize = 16)
	plt.title("%s" % (str(i[0] + " (" + i[1] + ")" + " - " + i[2] + " (" + i[1] + ")")), fontsize = 18)
	xticklines = getp(gca(), 'xticklines')
	yticklines = getp(gca(), 'yticklines')
	xgridlines = getp(gca(), 'xgridlines')
	ygridlines = getp(gca(), 'ygridlines')
	xticklabels = getp(gca(), 'xticklabels')
	ygridlines = getp(gca(), 'ygridlines')
	xticklabels = getp(gca(), 'xticklabels')
	yticklabels = getp(gca(), 'yticklabels')

	setp(xticklines, 'linewidth', 3)
	setp(yticklines, 'linewidth', 3)
	#setp(xgridlines, 'linestyle', '-')
	#setp(ygridlines, 'linestyle', '-')
	setp(yticklabels, 'color', 'Black', fontsize='medium')
	setp(xticklabels, 'color', 'Black', fontsize='medium')
	plt.grid(True)
	plt.savefig('%s.eps' % (str(i[0] + i[2])))
	plt.show()

	
		
		