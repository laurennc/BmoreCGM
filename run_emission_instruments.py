import yt
import numpy as np
import cPickle
from emission_instruments import *

base = "/Users/dalek/data/Molly/natural/nref11"
fn = base+"/RD0020/RD0020"
track_name = base+"/halo_track"
ds = yt.load(fn)
track = Table.read(track_name, format='ascii')
track.sort('col1')

lines = ['OVI','CIV','CIII_977','SiIV','HAlpha']

for line in lines:
  for instrument in instr_dict.keys():
    inst = instr_dict[instrument]
    for mode in inst.keys():
      emap = EmissionMap(ds,track,instrument,mode)
      frb = emap.make_frb(line,'x','instruments/emis_RD0020NAT_'+instrument+'_'+mode+'_'+line+'.cpkl')
      frb = frb['gas','Emission_'+line]
      emap.plot_frb('SiIV',frb,'instruments/emis_RD0020NAT_'+instrument+'_'+mode+'_'+line+'.pdf')
