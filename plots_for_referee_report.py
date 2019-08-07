import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import cPickle
import matplotlib as mpl
import trident
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from get_halo_center import *

from matplotlib.colorbar import Colorbar
from emission_functions import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.cosmology import WMAP9 as cosmo

import holoviews as hv
import pandas as pd
import datashader as dshade
from holoviews.operation.datashader import datashade, aggregate
from holoviews import Store
hv.extension('matplotlib')

#base = "/Volumes/sonic/halo_008508/nref11n/nref11f"
base = "/Users/lcorlies/data/Molly/nref11n_nref10f_refine200kpc_z4to2"
#base = "/Users/lcorlies/data/Molly/natural/nref11"
fn = base+"/RD0016/RD0016"
lines = ['OVI','CIV','CIII_977','SiIV','HAlpha']
lines2 = ['O VI','C IV','C III _977','Si IV','H Alpha']
track_name = base+"/halo_track"
args = fn.split('/')

def _cooling_criteria(field,data):
    """
    Calculates criteria used in enzo CellFlaggingMethod = 7
    """
    return -1*data['cooling_time'] / ((data['dx']/data['sound_speed']).in_units('s'))

# Add field to any ds that gets loaded
yt.add_field(("gas","cooling_criteria"),function=_cooling_criteria,units=None)

ds = yt.load(fn)
track = Table.read(track_name, format='ascii')
track.sort('col1')
rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)
redshift = ds.current_redshift

box_width = ds.arr(rb_width,'code_length').in_units('kpc')

def make_holoviews_radial_profile():
    dens = np.log10(rb['H_nuclei_density'])
    temp = np.log10(rb['Temperature'])
    Zgas = np.log10(rb['metallicity'])
    cell_mass = rb['cell_mass'].in_units('Msun')
    cell_volume = rb['cell_volume'].in_units('kpc**3')
    x = rb['x']
    y = rb['y']
    z = rb['z']
    vz = rb['z-velocity'].in_units('km/s')

    halo_center = ds.arr(rb_center,'code_length')
    dist = np.sqrt((halo_center[0]-rb['x'])**2.+(halo_center[1]-rb['y'])**2.+(halo_center[2]-rb['z'])**2.).in_units('kpc')

    df = pd.DataFrame({'temp':temp, 'dens':dens, 'Zgas':Zgas,'cell_volume':cell_volume,
                        'x':x,'y':y,'z':z,'dist':dist,'cell_mass':cell_mass,'vz':vz})

    temp_dist = hv.Scatter(df,kdims=['dist'],vdims=['temp'],label="Temperature ")
    dens_dist = hv.Scatter(df,kdims=['dist'],vdims=['dens'],label='Hydrogen Number Density')
    metal_dist = hv.Scatter(df,kdims=['dist'],vdims=['Zgas'],label='Metallicity')
    vz_dist    = hv.Scatter(df,kdims=['dist'],vdims=['vz'],label='z Velocity')

    #    if weight_by == 'cell_mass':
    temp_shade = aggregate(hv.Scatter(df,['dist','temp']),y_range=(2,8.4),aggregator=dshade.sum('cell_mass'))
    temp_shade = temp_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.Reds))
    dens_shade = aggregate(hv.Scatter(df,['dist','dens']),y_range=(-7,2.5),aggregator=dshade.sum('cell_mass'))
    dens_shade = dens_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.Blues))
    metal_shade = aggregate(hv.Scatter(df,['dist','Zgas']),y_range=(-7,2.5),aggregator=dshade.sum('cell_mass'))
    metal_shade = metal_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.BuGn))
    vz_shade    = aggregate(hv.Scatter(df,['dist','vz']),y_range=(-500.,500.),aggregator=dshade.sum('cell_mass'))#,dynamic=False)
    vz_shade    = vz_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.Purples))

    dist_plots = (temp_shade + dens_shade + metal_shade + vz_shade)
    fileout = 'basic_profile_cell_mass_'+args[-3]+'_'+args[-1]

    renderer = Store.renderers['matplotlib'].instance(fig='pdf', holomap='gif')
    renderer.save(dist_plots, fileout)
    return 'Finished radial profile!'

def plot_cooling_criteria():
    plt.hist(np.log10(rb[('gas','cooling_criteria')]*-1),alpha=0.65,label='all',range=[-3,6])
    idx = np.where(np.log10(rb['gas','cooling_criteria']*-1) < 0.)[0]
    title_phrase = '# of cells to trigger refinement: '+str(len(idx))+'/'+str(len(rb[('gas','cooling_criteria')]))
    title_phrase += ' = '+str(round(float(len(idx)) / len(rb['gas','cooling_criteria']), 2))
    #plt.hist(np.log10(rb[('gas','cooling_criteria')][idx]*-1),alpha=0.65,label='log(nH) > 0',range=[-3,6])
    plt.axvline(x=0.0,ls='--',lw=2.0,color='#5F6A6A')
    # plt.legend()
    plt.xlabel('log(cooling time / sound crossing time)')
    plt.title(title_phrase)
    fileout = 'cooling_criteria_histogram_'+args[-3]+'_'+args[-1]
    plt.savefig(fileout+'.pdf')
    plt.close()
    print 'Finished Cooling Criteria Histogram!'
    return

#make_holoviews_radial_profile()
plot_cooling_criteria()
