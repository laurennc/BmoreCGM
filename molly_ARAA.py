import yt
import matplotlib.pyplot as plt
import numpy as np


fn = "/Volumes/TOSHIBA EXT/r0058_l10/redshift0058"
ds = yt.load(fn,file_style="%s.grid.cpu%%04i")

pos = [0.40328598,0.47176743,0.46131516]
rvir = 316.86

sph = ds.sphere(pos,(rvir,'kpc'))

nhT = yt.PhasePlot(sph,('gas','H_number_density'),('gas','Temperature'),['cell_mass'],weight_field=None,x_bins=256,y_bins=256,figure_size=12,fontsize=22)

grey = sns.blend_palette(("#fefefe","#eeeeee","#dddddd","#aaaaaa","#888888","#666666","black"), as_cmap=True)

nhT.set_cmap("cell_mass",grey)
nhT.set_xlim(10**-6.5,10**2.5)
nhT.set_ylim(10**3.75,10**6.8)
nhT.set_zlim("cell_mass",1e37,1e41)
bottom_color="#ffffff"
mpl_plot = nhT.plots['cell_mass']
mpl_plot.axes.set_axis_bgcolor(bottom_color)
nhT.save('phaseplot_cellmass_256bins')

nhT = yt.PhasePlot(sph,('gas','H_number_density'),('gas','Temperature'),['cell_mass'],weight_field=None,x_bins=64,y_bins=64,figure_size=12,fontsize=22)
nhT.set_cmap("cell_mass",grey)
nhT.set_xlim(10**-6.5,10**2.5)
nhT.set_ylim(10**3.75,10**6.8)
nhT.set_zlim("cell_mass",1e37,1e41)
bottom_color="#ffffff"
mpl_plot = nhT.plots['cell_mass']
mpl_plot.axes.set_axis_bgcolor(bottom_color)
nhT.save('phaseplot_cellmass_64bins')

nhT = yt.PhasePlot(sph,('gas', 'H_number_density'),'Temperature',('gas', 'H_p0_fraction'),weight_field=None, x_bins=256, y_bins=256, figure_size=12)
nhT.set_cmap("H_p0_fraction",grey)
nhT.set_xlim(10**-6.5,10**2.5)
nhT.set_ylim(10**3.75,10**6.8)
nhT.set_zlim("H_p0_fraction",1e-4,1)
bottom_color="#ffffff"
mpl_plot = nhT.plots['H_p0_fraction']
mpl_plot.axes.set_axis_bgcolor(bottom_color)
nhT.save('phaseplot_H0frac_256bins')


nhT = yt.PhasePlot(sph,('gas', 'H_number_density'),'Temperature',('gas', 'H_p0_fraction'),weight_field=None, x_bins=64, y_bins=64, figure_size=12)
nhT.set_cmap("H_p0_fraction",grey)
nhT.set_xlim(10**-6.5,10**2.5)
nhT.set_ylim(10**3.75,10**6.8)
nhT.set_zlim("H_p0_fraction",1e-4,1)
bottom_color="#ffffff"
mpl_plot = nhT.plots['H_p0_fraction']
mpl_plot.axes.set_axis_bgcolor(bottom_color)
nhT.save('phaseplot_H0frac_64bins')

