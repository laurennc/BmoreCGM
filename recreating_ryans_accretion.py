import yt
import numpy as np
import get_halo_center as gc
from yt.funcs import mylog
import matplotlib.pyplot as plt
mylog.setLevel(40) # This sets the log level to "ERROR"

## Load in the data
fn = '/Users/dalek/data/Ryan/r0058_l10/redshift0058'
ds = yt.load(fn, file_style="%s.grid.cpu%%04i")
val,pos = ds.find_max('Density')


# Figure out the cell-size Ryan was using and set the width of my box
sp = ds.sphere(pos,(275.,'kpc'))
cell_size = np.unique(sp['dx'].in_units('kpc'))[2]
width = (cell_size*550).in_units('code_length')

# Set up the FRB. You def have to set the center for all of the
# following radial quantities
xL,xR = pos[0]-width/2.,pos[0]+width/2.
yL,yR = pos[1]-width/2.,pos[1]+width/2.
zL,zR = pos[2]-width/2.,pos[2]+width/2.
box = ds.r[xL:xR:551j,yL:yR:551j,zL:zR:551j]
box.set_field_parameter("center",pos)

# Everything works better if you just load the data.
# I think I was having memory issues because of how I was
# accessing the cube structures with where
print 'radius'
cell_rad = box['radius'].in_units('kpc').flatten()
print 'temperature'
cell_temp = box['temperature'].flatten()
print 'radial velocity'
cell_vr   = box['radial_velocity'].in_units('kpc/yr').flatten()
print 'cell mass'
cell_mass = box['cell_mass'].in_units('Msun').flatten()

## Start the loop through of all of the different masses that I want
radii = np.linspace(10,200.,91)
i = 0
cold_masses = np.zeros(len(radii))
warm_masses = np.zeros(len(radii))
hot_masses  = np.zeros(len(radii))

while i < len(radii)-1:
    maxr = radii[i+1]
    minr = radii[i]
    idC = np.where((cell_rad >= minr) & (cell_rad < maxr)  &
                   (box['temperature'] < 1e5))[0]
    idW = np.where((cell_rad >= minr) & cell_rad < maxr)  &
                   (temp < 1e6) &
                   (temp >=1e5 ))[0]
    idH = np.where((cell_rad >= minr) & (cell_rad < maxr)  &
                   (box['temperature'] >= 1e6))[0]
    cold_masses[i] = np.sum(cell_mass[idC]*cell_vr[idC])
    warm_masses[i] = np.sum(cell_mass[idW]*cell_vr[idW])
    hot_masses[i]  = np.sum(cell_mass[idH]*cell_vr[idH])
    i = i + 1

masses = {}


dr = radii[1]-radii[0]
#print dr
plt.plot(radii,-1.*cold_masses/(dr),color='b',label='cold')
plt.plot(radii,-1.*warm_masses/(dr),color='g',label='warm')
plt.plot(radii,-1.*hot_masses/(dr),color='r',label='hot')
plt.legend()
#plt.ylim(0,20)
plt.show()
## Now let's create some plots like Ryan's. That's a little difficult
## because of the nice way he's stacked and shaded the plots.
