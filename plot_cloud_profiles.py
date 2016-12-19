import numpy as np
import matplotlib.pyplot as plt

CloudRadius      =  0.0625 *1.6e17
CloudRatio       =  0.7
CloudInnerRadius =  CloudRatio * CloudRadius
UniformDensity   = 1.67e-24
CloudDensity     =  1000.*UniformDensity

r = np.linspace(0.03,0.07,1e4)
r = r * 1.6e17

sigma = np.sqrt( (-1*(CloudRadius-CloudInnerRadius)**2.0) / (2.0*np.log10(UniformDensity/CloudDensity) ) )

#ct == 0
i = 0
dens0 = np.zeros(len(r)) + UniformDensity

for i in range(len(r)):
	if r[i] < CloudRadius:
		dens0[i] = CloudDensity

#ct == 1
i = 0
dens1 = np.zeros(len(r)) + UniformDensity

for i in range(len(r)):
	if r[i] < CloudRadius:
		if r[i] <= CloudInnerRadius:
			dens1[i] = CloudDensity
		if r[i] > CloudInnerRadius:
			dens1[i] = (UniformDensity - CloudDensity)/( (1.0/CloudRadius) - (1./CloudInnerRadius) )*( (1./r[i]) - (1./CloudInnerRadius) ) + CloudDensity


#ct == 2
i = 0
dens2 = np.zeros(len(r)) + UniformDensity

for i in range(len(r)):
	if r[i] < CloudRadius:
		if r[i] <= CloudInnerRadius:
			dens2[i] = CloudDensity
		if r[i] > CloudInnerRadius:
			dens2[i] = CloudDensity * np.exp( (-1.0*(r[i]-CloudInnerRadius)**2.0 ) / (2.*sigma**2.0) )

#plt.plot(r,np.log10(dens0),'r',lw=2.0,label='Uniform')
#plt.plot(r,np.log10(dens1),'b',lw=2.0,label='Power Law')
#plt.plot(r,np.log10(dens2),'g',lw=2.0,label='Gaussian')
plt.plot(r,dens0,'r',lw=2.0,label='Uniform')
plt.plot(r,dens1,'b',lw=2.0,label='Power Law')
plt.plot(r,dens2,'g',lw=2.0,label='Gaussian')
plt.xlabel('Radius')
plt.ylabel("Density")
plt.legend(loc=3)
plt.savefig('cloud_profiles.png')
