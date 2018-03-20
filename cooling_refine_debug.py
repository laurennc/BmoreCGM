import matplotlib
matplotlib.use('Agg')

import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import cPickle
from emission_functions import *

def _cooling_criteria(field,data):
    return data['cooling_time'] / ((data['dx']/data['sound_speed']).in_units('s'))

yt.add_field(("gas","cooling_criteria"),function=_cooling_criteria,units=None)

#ds = yt.load('DD0056/DD0056')
#track = Table.read('nref11_track_box',format='ascii')
ds = yt.load('DD0131/DD0131')
track = Table.read('halo_track',format='ascii')
track.sort('col1')
rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)
rb_width = ds.arr(rb_width,'code_length').in_units('kpc')

track_small = Table.read('nref11_cool_small_region',format='ascii')
track_small.sort('col1')
rb_small,rb_center_small,rb_width_small = get_refine_box(ds,ds.current_redshift,track_small)

diff = track['col1'] - ds.current_redshift
this_loc = track[np.where(diff == np.min(diff[np.where(diff > 1.e-6)]))]
xL,yL,zL = this_loc['col2'][0],this_loc['col3'][0],this_loc['col4'][0]
xR,yR,zR = this_loc['col5'][0],this_loc['col6'][0],this_loc['col7'][0]

diff = track_small['col1'] - ds.current_redshift
this_loc = track_small[np.where(diff == np.min(diff[np.where(diff > 1.e-6)]))]
xLs,yLs,zLs = this_loc['col2'][0],this_loc['col3'][0],this_loc['col4'][0]
xRs,yRs,zRs = this_loc['col5'][0],this_loc['col6'][0],this_loc['col7'][0]

zc = rb_center[2]
ray1 = ds.ray([xL,yL,zc],[xR,yL,zc])
ray2 = ds.ray([xR,yL,zc],[xR,yR,zc])
ray3 = ds.ray([xR,yR,zc],[xL,yR,zc])
ray4 = ds.ray([xL,yR,zc],[xL,yL,zc])

ray1s = ds.ray([xLs,yLs,zLs],[xRs,yLs,zLs])
ray2s = ds.ray([xRs,yLs,zLs],[xRs,yRs,zLs])
ray3s = ds.ray([xRs,yRs,zLs],[xLs,yRs,zLs])
ray4s = ds.ray([xLs,yRs,zLs],[xLs,yLs,zLs])


proj = yt.ProjectionPlot(ds,'z',('index','grid_level'),center=rb_center,width=(100.,'kpc'),method='mip')
proj.set_log(('index','grid_level'),False)
proj.annotate_ray(ray1)
proj.annotate_ray(ray2)
proj.annotate_ray(ray3)
proj.annotate_ray(ray4)
#proj.annotate_ray(ray1s)
#proj.annotate_ray(ray2s)
#proj.annotate_ray(ray3s)
#proj.annotate_ray(ray4s)
proj.save()


## cooling criteria ##
proj = yt.SlicePlot(ds,'z',('gas','cooling_criteria'),center=rb_center,width=(100.,'kpc'))
proj.set_cmap(('gas','cooling_criteria'),'bwr')
proj.set_zlim(('gas','cooling_criteria'),1e-4,1e4)
#proj.set_zlim(('gas','cooling_criteria'),1e-4,1)
proj.annotate_ray(ray1,plot_args={'color':'k'})
proj.annotate_ray(ray2,plot_args={'color':'k'})
proj.annotate_ray(ray3,plot_args={'color':'k'})
proj.annotate_ray(ray4,plot_args={'color':'k'})
#proj.annotate_ray(ray1s)
#proj.annotate_ray(ray2s)
#proj.annotate_ray(ray3s)
#proj.annotate_ray(ray4s)
#proj.save('DD0131_cool2')
proj.save()

proj = yt.SlicePlot(ds,'z',('index','grid_level'),center=rb_center,width=(100.,'kpc'))
proj.set_log(('index','grid_level'),False)
proj.annotate_ray(ray1,plot_args={'color':'k'})
proj.annotate_ray(ray2,plot_args={'color':'k'})
proj.annotate_ray(ray3,plot_args={'color':'k'})
proj.annotate_ray(ray4,plot_args={'color':'k'})
#proj.annotate_ray(ray1s)
#proj.annotate_ray(ray2s)
#proj.annotate_ray(ray3s)
#proj.annotate_ray(ray4s)
proj.save()

### What volume is refined to the max level versus the (max level - 1)
total_volume = 374693.882872 kpc**3
nref9_volume =
nref10_volume =

nref9_numcell = 2665120
nref10_numcell =
