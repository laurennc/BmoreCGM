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

import seaborn as sns
sns.set_style("white", {'axes.grid' : False})

def Stars(pfilter, data):
      return data[("all", "particle_type")] == 2



add_particle_filter("stars", function=Stars, filtered_type='all',
                    requires=["particle_type"])

def Starslow(pfilter,data):
    #idx = (np.where(data["all","particle_type"]) == 2) & (np.where(data["all","particle_mass"].in_units('Msun') < 5e4))
    return ((data["all","particle_type"] == 2) & (data["all","particle_mass"].in_units('Msun') < 5e4 ))

add_particle_filter("stars_low",function=Starslow,filtered_type='all',requires=["particle_type"])

h1_color_map = sns.blend_palette(("white","#ababab","#565656","black","#4575b4","#984ea3","#d73027","darkorange","#ffe34d"), as_cmap=True)
h1_proj_min = 1.e12
h1_proj_max = 1.e24
metal_color_map = sns.blend_palette(("black","#984ea3","#4575b4","#4daf4a","#ffe34d","darkorange"), as_cmap=True)
metal_min = 1.e-4
metal_max = 2.

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
    refine_box_center = [0.5*(x_left+x_right), 0.5*(y_left+y_right), 0.5*(z_left+z_right)]
    refine_box = ds.r[x_left:x_right, y_left:y_right, z_left:z_right]
    refine_width = np.abs(x_right - x_left)
    return refine_box, refine_box_center, refine_width


#base = "/lou/s2m/mpeeples/halo_008508/nref11n/nref11f_refine200kpc_z4to2"
#base = "/lou/s2m/mpeeples/halo_008508/nref11n/nref11n_nref10f_refine200kpc_z4to2"
base = "nref11n/natural/"
#fn = base+"/RD0017/RD0017"
track_name = "/lou/s2m/lcorlies/halo_track"
args = fn.split('/')

fn= 'RD0042/RD0042'
fn = 'DD0020/DD0020'
#track_name = 'halo_track'
ds = yt.load(fn)
ds.add_particle_filter('stars')
ds.add_particle_filter('stars_low')
track = Table.read(track_name, format='ascii')
track.sort('col1')
rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)
redshift = ds.current_redshift

box_width = ds.arr(rb_width,'code_length').in_units('kpc')

proj = yt.ProjectionPlot(ds,'z','H_number_density',center=rb_center,width=(float(box_width.value),'kpc'))
proj.set_cmap('H_number_density',h1_color_map)
proj.set_zlim('H_number_density',h1_proj_min,h1_proj_max)
proj.annotate_particles((500,'kpc'),col='MediumSeaGreen',ptype='stars',alpha=0.4)
proj.annotate_particles((500,'kpc'),col='#33F3FF',ptype='stars_low')
proj.save('nref11n_DD0020_HI_w_stars2')


proj = yt.ProjectionPlot(ds,'y','dark_matter_density',center=rb_center,width=(float(box_width.value),'kpc'))
proj.set_unit('dark_matter_density','Msun/kpc**2')
proj.save('nref11n_10f_RD0042')


proj = yt.ProjectionPlot(ds,'y','metallicity',center=rb_center,width=(float(box_width.value),'kpc'),weight_field='Density')
proj.set_cmap('metallicity',metal_color_map)
proj.set_zlim('metallicity',metal_min,metal_max)
proj.annotate_particles((float(box_width.value)*2,'kpc'),col='MediumSeaGreen',ptype='stars)
proj.save(FILEOUT)

metal_color_map = sns.blend_palette(("black","#984ea3","#4575b4","#4daf4a","#ffe34d","darkorange"), as_cmap=True)
metal_min = 1.e-4
metal_max = 2.
