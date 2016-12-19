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

def plot_shocks(fout):
	dens_c = [100.,1000.,1e4]
	mach = [1.,5.,10.,15.]
	temp_c = [1e3,1e4]

	dc_arr,T_arr,mach_arr = np.meshgrid(dens_c,temp_c,mach)
	pts = np.array((dc_arr.ravel(),T_arr.ravel(),mach_arr.ravel() )).T
	
	dc, Tc, mc = pts[:,0],pts[:,1],pts[:,2]
	cs = sound_speed(Tc,dc)
	
	shockV = (cs*mc).to('AU/yr')
	shockD = shockV*(400.*u.yr)
	
	shockD = shockD.value
	shockD = np.log10(shockD)
	shockV = shockV.value
	shockV = np.log10(shockV)

	fig,ax = plt.subplots(1,2)
	ax[0].scatter(np.log10(dc_arr),shockV,c=Tc,alpha=0.5)
	ax[0].set_xlabel('Density Constrast')
	ax[0].set_ylabel('log(Shock Velocity) (AU/yr)')
	
	ax[1].scatter(np.log10(dc_arr),shockD,c=Tc,alpha=0.5)
	ax[1].set_xlabel('Density Contrast')
	ax[1].set_ylabel('log(Distance Traveled) (AU)')
	
	plt.savefig(fout)
	return "made plots!"

fout = 'shocks_in_cycle.png'
plot_shocks(fout)



