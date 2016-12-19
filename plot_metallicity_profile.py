import numpy as np
import matplotlib.pyplot as plt

CloudRadius      =  0.0625 *1.6e17
CloudRatio       =  0.7
CloudInnerRadius =  CloudRatio * CloudRadius
UniformDensity   = 1.67e-24
CloudDensity     =  1000.*UniformDensity
 
r = np.linspace(0.03,0.07,1e4)
r = r * 1.6e17

UniformMetallicityFraction = 0.001
CloudMetallicityFraction   = 0.01

metal_density = np.zeros(len(r)) + UniformMetallicityFraction * UniformDensity 

for i in range(len(r)):
  if r[i] < CloudRadius:
    if r[i] <= CloudInnerRadius:
      metal_density[i] = CloudDensity * CloudMetallicityFraction
    if r[i] > CloudInnerRadius:
      metal_density[i] = ((UniformMetallicityFraction*UniformDensity - CloudMetallicityFraction*CloudDensity)/((1.0/CloudRadius) - (1.0/CloudInnerRadius))) * ((1.0/r[i]) - (1.0/CloudInnerRadius)) + CloudMetallicityFraction * CloudDensity


sigmaZ = np.sqrt((-1.0*(CloudRadius - CloudInnerRadius)**2.0 )/(2.0*np.log10( (UniformDensity*UniformMetallicityFraction) / (CloudDensity * CloudMetallicityFraction) ) ))


for i in range(len(r)):
   if r[i] < CloudRadius:
     if r[i] <= CloudInnerRadius:
       metal_density[i] = CloudDensity * CloudMetallicityFraction
     if r[i] > CloudInnerRadius:
       metal_density[i] = CloudDensity * CloudMetallicityFraction * np.exp((-1.0*(r[i] - CloudInnerRadius)**2.0 )/(2.0*sigmaZ**2.0))
 


plt.plot(r,metal_density,'b',lw=2.0,label='Power Law')
plt.savefig('metallicity_profiles.png')

