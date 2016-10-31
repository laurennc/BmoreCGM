import yt
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from cgs import *
from astropy import units as u

gamma = 1.6667

## Method for calculating the sound speed in the idealized Cloud Crusher simulations
## Assumes pressure equilibrium for an ideal gas
def sound_speed(Tc,delta_cloud):
	prefactor = (gamma*kb)/(mu*mp)
	return np.sqrt(prefactor*Tc*delta_cloud)*u.cm/u.s  #cm/s


