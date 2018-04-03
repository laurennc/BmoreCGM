import matplotlib
matplotlib.use('Agg')

import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import cPickle

def get_refine_box(ds, zsnap, track):
    ## find closest output, modulo not updating before printout
    #diff = track['col1'] - zsnap
    #this_loc = track[np.where(diff == np.min(diff[np.where(diff > 1.e-6)]))]
    #print "using this loc:", this_loc
    x_left = np.interp(zsnap, track['col1'], track['col2'])
    y_left = np.interp(zsnap, track['col1'], track['col3'])
    z_left = np.interp(zsnap, track['col1'], track['col4'])
    x_right = np.interp(zsnap, track['col1'], track['col5'])
    y_right = np.interp(zsnap, track['col1'], track['col6'])
    z_right = np.interp(zsnap, track['col1'], track['col7'])
    refine_box_center = [0.5*(x_left+x_right), 0.5*(y_left+y_right), 0.5*(z_left+z_right)]
    refine_box = ds.r[x_left:x_right, y_left:y_right, z_left:z_right]
    refine_width = np.abs(x_right - x_left)
    return refine_box, refine_box_center, refine_width

def _cooling_criteria(field,data):
    return data['cooling_time'] / ((data['dx']/data['sound_speed']).in_units('s'))

yt.add_field(("gas","cooling_criteria"),function=_cooling_criteria,units=None)

#ds = yt.load('RD0016/RD0016')
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


proj = yt.SlicePlot(ds,'z',('index','grid_level'),center=rb_center,width=(100.,'kpc'))
proj.set_log(('index','grid_level'),False)
proj.annotate_ray(ray1,plot_args={'color':'k'})
proj.annotate_ray(ray2,plot_args={'color':'k'})
proj.annotate_ray(ray3,plot_args={'color':'k'})
proj.annotate_ray(ray4,plot_args={'color':'k'})
proj.annotate_grids()
#proj.annotate_ray(ray1s)
#proj.annotate_ray(ray2s)
#proj.annotate_ray(ray3s)
#proj.annotate_ray(ray4s)
proj.save()

## cooling criteria ##
proj = yt.SlicePlot(ds,'x',('gas','cooling_criteria'),center=rb_center,width=(100.,'kpc'))
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


### What volume is refined to the max level versus the (max level - 1)
total_volume = 374693.882872 kpc**3
nref9_volume = 55261.0643739 kpc**3  -- 15%
nref10_volume = 319432.818498 kpc**3 -- 85%

nref9_numcell =  2665120   - 2%
nref10_numcell = 123244357 - 98%



proj = yt.ProjectionPlot(ds,'x',('index','grid_level'),center=rb_center,width=(100.,'kpc'),method='mip')
proj.set_log(('index','grid_level'),False)
proj.save()

proj = yt.SlicePlot(ds,'x',('index','grid_level'),center=rb_center,width=(100.,'kpc'))
proj.annotate_grids()
proj.set_log(('index','grid_level'),False)
proj.save()

proj = yt.SlicePlot(ds,'z','dz',center=rb_center,width=(100.,'kpc'))
proj.set_log(('index','grid_level'),False)
proj.save()
