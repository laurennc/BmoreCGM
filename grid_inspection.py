import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

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


fn = '/nobackupnfs2/lcorlies/cooltest/nref10f_sym50kpc/mustrefine/DD0101/DD0101'
track_name = '/nobackupnfs2/lcorlies/cooltest/nref10f_sym50kpc/nref11_track_box'
ds = yt.load(fn)
track = Table.read(track_name, format='ascii')
track.sort('col1')
rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)

proj = yt.ProjectionPlot(ds,'x',('index','grid_level'),data_source=rb,method='mip',
                         center=rb_center,width=(rb_width,'code_length'))
proj.set_log(('index','grid_level'),False)
proj.save('mustrefine_mip_gridlevel')

sl = yt.SlicePlot(ds,'x','Density',data_source=rb,
                         center=rb_center,width=(rb_width,'code_length'))
sl.annotate_cells()
sl.save('mustrefine_dens_cells')
