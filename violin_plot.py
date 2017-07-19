import yt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_style("whitegrid", {'axes.grid' : False})

## For the z=2 simulations!
ds = yt.load('RD0020/RD0020')
center = np.array([0.493970873304, 0.488697442131, 0.502242657694])
sp = ds.sphere(center,(45.,'kpc'))


dx = np.log10(sp['dx'].in_units('pc'))
grid = sp[('index','grid_level')]
mass = np.log10(sp['cell_mass'].in_units('Msun'))

temp = np.log10(sp['Temperature'])
idx = np.where(temp < 5.)[0]
temp_binary = np.zeros(len(temp))
temp_binary[idx] = 1.0

data = {'cell_size':dx,
        'cell_mass':mass,
        'grid_level':grid,
        'temp_binary':temp_binary}
df = pd.DataFrame(data=data)

ax = sns.violinplot(x=df['grid_level'],y=df['cell_mass'],scale='count') #  ,hue=df['temp_binary'],split=True)
ax.set_ylabel('Cell Mass [Msun]')
ax.set_xlabel('Grid Level')
ax.invert_xaxis()
plt.savefig('violin_resolution_mass_count_highres_RD0020.png')
plt.close()



##############################################################################
## FOR TWO SIMULATIONS
## For the z=2 simulations!
ds1 = yt.load('nref11/RD0020/RD0020')
ds2 = yt.load('nref11_refine200kpc_z4to2/RD0020/RD0020')
center = np.array([0.493970873304, 0.488697442131, 0.502242657694])
sp1 = ds1.sphere(center,(45.,'kpc'))
sp2 = ds2.sphere(center,(45.,'kpc'))

dx1 = np.log10(sp1['dx'].in_units('pc'))
grid1 = sp1[('index','grid_level')]
mass1 = np.log10(sp1['cell_mass'].in_units('Msun'))
sim_ids1 = np.zeros(len(dx1))

dx2 = np.log10(sp2['dx'].in_units('pc'))
grid2 = sp2[('index','grid_level')]
mass2 = np.log10(sp2['cell_mass'].in_units('Msun'))
sim_ids2 = np.zeros(len(dx2))
sim_ids2 = sim_ids2 + 1

dx      = np.concatenate((dx1,dx2))
grid    = np.concatenate((grid1,grid2))
mass    = np.concatenate((mass1,mass2))
sim_ids = np.concatenate((sim_ids1,sim_ids2))

data = {'cell_size':dx,
        'cell_mass':mass,
        'grid_level':grid,
        'sim_ids':sim_ids}
        #'temp_binary':temp_binary}
df = pd.DataFrame(data=data)

##############################################################################

ax = sns.violinplot(x=df['grid_level'],y=df['cell_mass'],scale='count',
                    hue=df['sim_ids'],palette='husl')#,split=True)
ax.set_ylabel('Cell Mass [Msun]')
ax.set_xlabel('Grid Level')
ax.invert_xaxis()
plt.savefig('violin_resolution_mass_RD0020.png')
plt.close()
