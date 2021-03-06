import yt
import numpy as np
import cPickle
from emission_instruments import *

base = "/Users/dalek/data/Molly/nref11n_nref10f_refine200kpc_z4to2"
#base = "/Users/dalek/data/Molly/natural/nref11"
fn = base+"/RD0016/RD0016"
track_name = base+"/halo_track"
ds = yt.load(fn)
track = Table.read(track_name, format='ascii')
track.sort('col1')

lines = ['OVI','CIV','CIII_977','SiIV','HAlpha']
#lines = ['HAlpha']

for line in lines:
  #for instrument in instr_dict.keys():
    #inst = instr_dict[instrument]
    instrument = 'MUSE'
    inst = instr_dict['MUSE']
    #instrument = 'KCWI'
    #inst = instr_dict['KCWI']
    for mode in inst.keys():
      emap = EmissionMap(ds,track,instrument,mode)
      frb = emap.make_frb(line,'z','instruments/emis_z_RD0016REF_'+instrument+'_'+mode+'_'+line+'.cpkl')
      frb = frb['gas','Emission_'+line]
      emap.plot_frb_talks(line,frb,instrument,'instruments/green_obs/emis_z_RD0016REF_'+instrument+'_'+mode+'_'+line+'.pdf')
