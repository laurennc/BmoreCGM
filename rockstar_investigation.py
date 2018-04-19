import matplotlib as mpl
mpl.use('agg')

import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import cPickle
import matplotlib as mpl
#import trident
from yt.data_objects.particle_filters import add_particle_filter

def Stars(pfilter, data):
      return data[("all", "particle_type")] == 2

add_particle_filter("stars", function=Stars, filtered_type='all',
                    requires=["particle_type"])

def get_refine_box(ds, zsnap, track):
    ## find closest output, modulo not updating before printout
    diff = track['col1'] - zsnap
    this_loc = track[np.where(diff == np.min(diff[np.where(diff > 1.e-6)]))]
    print "using this loc:", this_loc
    x_left = this_loc['col2'][0]
    y_left = this_loc['col3'][0]
    z_left = this_loc['col4'][0]
    x_right = this_loc['col5'][0]
    y_right = this_loc['col6'][0]
    z_right = this_loc['col7'][0]
    LE = [x_left,y_left,z_left]
    RE = [x_right,y_right,z_right]
    refine_box_center = [0.5*(x_left+x_right), 0.5*(y_left+y_right), 0.5*(z_left+z_right)]
    refine_box = ds.r[x_left:x_right, y_left:y_right, z_left:z_right]
    refine_width = np.abs(x_right - x_left)
    return refine_box, refine_box_center, refine_width,LE,RE

def within_refine_region(x,y,z,LE,RE):
  if ((x < RE[0]) & (x>LE[0])):
    pass
  else:
    return False
  if ((y < RE[1]) & (y>LE[1])):
    pass
  else:
    return False
  if ((z < RE[2]) & (z>LE[2])):
    return True
  else:
    return False

fn= 'RD0018/RD0018'
track_name = 'halo_track'
ds = yt.load(fn)
ds.add_particle_filter('stars')
track = Table.read(track_name, format='ascii')
track.sort('col1')
rb,rb_center,rb_width,LE,RE = get_refine_box(ds,ds.current_redshift,track)
redshift = ds.current_redshift
box_width = ds.arr(rb_width,'code_length').in_units('kpc')
halos = np.genfromtxt('rockstar_halos/out_4.list',skip_header=16)

xL,yL,zc = LE[0],LE[1],0.5
xR,yR,zc = RE[0],RE[1],0.5
ray1 = ds.ray([xL,yL,zc],[xR,yL,zc])
ray2 = ds.ray([xR,yL,zc],[xR,yR,zc])
ray3 = ds.ray([xR,yR,zc],[xL,yR,zc])
ray4 = ds.ray([xL,yR,zc],[xL,yL,zc])


## halo positions are in Mpc/h and virial radii are in kpc/h
## so need to convert all of these to proper units for plotting!
## starting only with the most massive halos

idm = np.where(np.log10(halos[:,2]) > 7.)[0]

#proj = yt.ProjectionPlot(ds,'z','Density',width=(10,'Mpc'))#,center=rb_center,width=(float(box_width.value),'kpc'))
proj = yt.ProjectionPlot(ds,'z','Density',center=rb_center,width=(20.,'kpc')) #width=(float(box_width.value),'kpc'))
proj.set_cmap('Density','jet')

for i in idm:
  #x = float(ds.quan(halos[i,8],'Mpccm').in_units('code_length').value)
  #y = float(ds.quan(halos[i,9],'Mpccm').in_units('code_length').value)
  #z = float(ds.quan(halos[i,10],'Mpccm').in_units('code_length').value)
  x = halos[i,8]/100.
  y = halos[i,9]/100.
  z = halos[i,10]/100.
  halo_center = [x,y,z]
  virial_radius = float(ds.quan(halos[i,5],'kpccm').in_units('kpc').value)
  proj.annotate_sphere(halo_center,radius=(virial_radius,'kpc'))#,circle_args={'color':'black'})
  #proj.annotate_text(halo_center,str(int(halos[i,0])))#,text_args={'color':'black'})

proj.annotate_ray(ray1,plot_args={'color':'red'})
proj.annotate_ray(ray2,plot_args={'color':'red'})
proj.annotate_ray(ray3,plot_args={'color':'red'})
proj.annotate_ray(ray4,plot_args={'color':'red'})
proj.save('confirming_rockstar_RD0018_startest_zoom')

## so this just looking for the massive halos thing totally didn't work
## SO
## let's find the halos that are actually within the refine region
for i in range(len(halos[:,0])):
  #x = float(ds.quan(halos[i,8],'Mpccm').in_units('code_length').value)
  #y = float(ds.quan(halos[i,9],'Mpccm').in_units('code_length').value)
  #z = float(ds.quan(halos[i,10],'Mpccm').in_units('code_length').value)

  if within_refine_region(x,y,z,LE,RE):
    print int(halos[i,0])
  else:
    continue
