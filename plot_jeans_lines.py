import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from jeans_calculation import *

Tas = [5.,5.5,6.,6.5]
nas = [-5,-4.5,-4.,-3.5,-3.]
nbs = [-3.5,-3.,-2.5,-2.,-1.5,-1.]

use_log = True
ylen,xlen = 1,3
fig,ax = plt.subplots(ylen,xlen)
fig.set_size_inches(11,4.5)

for Ta in Tas:
	for na in nas:
		Tbs = []
		jls = []
		css = []
		jms = []
		
		
		for nb in nbs:
			Tb = equilibrium_T(Ta,na,nb)
			Tbs = np.append(Tbs,Tb)

			cs = sound_speed(Tb).to(u.km/u.s)
			jl = jeans_length(cs,nb).to(u.kpc)
			jm = jeans_mass(jl,nb).to(u.solMass)

			css = np.append(css,cs)
			jls = np.append(jls,jl)
			jms = np.append(jms,jm)
		
		if use_log == True:	
				ax[0].plot(nbs,np.log10(css),'o')
				ax[0].set_ylabel('Sound Speed [log(km/s)]')

				ax[1].plot(nbs,np.log10(jls),'o')
				ax[1].set_ylabel('Jeans Length [log(kpc)]')
				ax[1].set_xlabel('Cloud H Density [log(cm^-3)]')
				ax[1].axhline(np.log10(0.7),lw=1.5,color='red')

				ax[2].plot(nbs,np.log10(jms),'o')
				ax[2].set_ylabel('Jeans Mass [log(Msun)]')

		if use_log == False:
				ax[0].plot(nbs,css,'o')
				ax[0].set_ylabel('Sound Speed [km/s]')

				ax[1].plot(nbs,jls,'o')
				ax[1].set_xlabel('Cloud H Density [cm^-3]')
				ax[1].axhline(0.7,lw=1.5,color='red')

				ax[2].plot(nbs,jms,'o')
				ax[2].set_ylabel('Jeans Mass [Msun]')



plt.savefig('explore_jeans.png')
plt.close()


