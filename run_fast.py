import yt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import AxesGrid


print 'finding most dense'

rvir0,rvir02,rvir05,rvir1 = 316.86,294.55,226.32,162.66
fn = "/Volumes/TOSHIBA EXT/r0058_l10/redshift0058"
ds = yt.load(fn,file_style="%s.grid.cpu%%04i")
val, pos = ds.find_max('Density')
rvir = rvir1

print 'density projection'

dp = yt.ProjectionPlot(ds,'x',('gas','density'),center=pos,width=(1000,"kpc"))
cmap = sns.blend_palette(("black","#984ea3","#d73027","darkorange","#ffe34d","#4daf4a","white"), n_colors=60, as_cmap=True)
dp.set_cmap(('gas','density'),cmap)
dp.set_zlim(('gas','density'),5e-6,1e-2)
dp.annotate_sphere(pos, radius=(rvir,'kpc'), circle_args={'color':'white', 'ls':'--', 'lw':2})
dp.hide_axes()
dp.save('density_z1')

print 'temperature projection'

tp = yt.ProjectionPlot(ds,'x',('gas','temperature'),weight_field=('gas','density'),center=pos,width=(1000,"kpc"))
cmap = sns.blend_palette(("black","#d73027","darkorange","#ffe34d"), n_colors=50, as_cmap=True)
tp.set_cmap("temperature",cmap)
tp.set_zlim("temperature",1e4,5e6)
tp.annotate_sphere(pos, radius=(rvir,'kpc'), circle_args={'color':'white', 'ls':'--', 'lw':2})
tp.hide_axes()
tp.save('temp_z1')


print 'metal high proj'

mp = yt.ProjectionPlot(ds, 'x', ('gas','metallicity'), weight_field=('gas','density'), center=pos, width=(1000,"kpc"))
mp.set_zlim(('gas','metallicity'),1e-2,1)
cmap = sns.blend_palette(("black","#984ea3","#4575b4","#4daf4a","#ffe34d","darkorange"), as_cmap=True)
mp.set_cmap(('gas','metallicity'),cmap)
mp.annotate_sphere(pos, radius=(rvir,'kpc'), circle_args={'color':'white', 'ls':'--', 'lw':2})
mp.hide_axes()
mp.save('metal_z1_higher')


print 'metal proj'

mp = yt.ProjectionPlot(ds, 'x', ('gas','metallicity'), weight_field=('gas','density'), center=pos, width=(1000,"kpc"))
mp.set_zlim(('gas','metallicity'),1e-4,1e-2)
cmap = sns.blend_palette(("black","#984ea3","#4575b4","#4daf4a","#ffe34d","darkorange"), as_cmap=True)
mp.set_cmap(('gas','metallicity'),cmap)
mp.annotate_sphere(pos, radius=(rvir,'kpc'), circle_args={'color':'white', 'ls':'--', 'lw':2})
mp.hide_axes()
mp.save('metal_z1')



