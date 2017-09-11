import yt
import numpy as np
import matplotlib.pyplot as plt

fn = "PATH/r0054/redshift0054"
ds = yt.load(fn,file_style="%s.grid.cpu%%04i")

##The easiest way to find the high-resolution galaxy is to look for the point of highest density in the simulation

val,pos = ds.find_max('Density')

#from previous code, I know that
rvir = 294.55 #kpc

#Create the Projection object for a Density Projection
dp = yt.ProjectionPlot(ds,'x',('gas','density'),center=pos,width=(1000,"kpc"))
#set the color map
dp.set_cmap(('gas','density'),'jet')
#set the color map limits
dp.set_zlim(('gas','density'),5e-6,1e-2)
#annotate a sphere at the virial radius of the halo
dp.annotate_sphere(pos, radius=(rvir,'kpc'), circle_args={'color':'white', 'ls':'--', 'lw':2})
dp.hide_axes()
dp.save('density_z02')


##Same deal for the temperature but now the field is weighted by the gas density
tp = yt.ProjectionPlot(ds,'x',('gas','temperature'),weight_field=('gas','density'),center=pos,width=(1000,"kpc"))
cmap = sns.blend_palette(("black","#d73027","darkorange","#ffe34d"), n_colors=50, as_cmap=True)
tp.set_cmap("temperature",cmap)
tp.set_zlim("temperature",1e4,5e6)
tp.annotate_sphere(pos, radius=(rvir,'kpc'), circle_args={'color':'white', 'ls':'--', 'lw':2})
tp.hide_axes()
tp.save('temp_z02')

## Metallicity is also density-weighted
mp = yt.ProjectionPlot(ds, 'x', ('gas','metallicity'), weight_field=('gas','density'), center=pos, width=(1000,"kpc"))
mp.set_zlim(('gas','metallicity'),1e-4,1e-2)
cmap = sns.blend_palette(("black","#984ea3","#4575b4","#4daf4a","#ffe34d","darkorange"), as_cmap=True)
mp.set_cmap(('gas','metallicity'),cmap)
mp.annotate_sphere(pos, radius=(rvir,'kpc'), circle_args={'color':'white', 'ls':'--', 'lw':2})
mp.hide_axes()
mp.save('metal_z1')


