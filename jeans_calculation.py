import numpy as np
import matplotlib.pyplot as plt
import math
from astropy import units as u

mH = 1.67e-24 * u.g
G = 6.67259e-8 *(u.cm**3) / u.g / (u.s**2)
k = 1.380658e-16 * u.erg / u.K
gamma = 1.6667

#def sound_speed(temp):
#	T = 10.**temp * u.K
#	return  np.sqrt(gamma*k*T/mH)
	
def sound_speed(temp):
	T = 10.**temp * u.K
	m = 0.59*mH
	return np.sqrt((5.*k*T)/(3*m))

def jeans_length(cs,hden):
	hden = 10.**hden * (u.cm**-3)
	numerator = math.pi * cs**2
	denom = G * mH * hden
	return np.sqrt(numerator/denom)

#def jeans_mass(temp,hden):
#	T = 10.**temp * u.K
#	hden = 10.**hden * u.cm**-3
#	factor1 = ((math.pi*k*T)/(G*1.2))**1.5
#	factor2 = (hden*mH))**-0.5
#	return (1./8.)*factor1*factor2

def jeans_mass(j_length,hden):
	hden = 10.**hden * u.cm**-3
	factor = (j_length/2.)**3.
	return (4.*math.pi/3.)*hden*mH*factor

## want to make sure I'm getting reasonable column densities...

def equilibrium_T(Ta,na,nb):
	factor = ((10.**Ta)*(10.**na))/(10.**nb)
	return np.log10(factor) 

