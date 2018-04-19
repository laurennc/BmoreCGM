import yt
import numpy as np
import matplotlib.pyplot as plt
import cPickle
import trident
from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
from emission_functions import *

# ok so general structure is that I'm going to have to make an FRB
## for each instrument for each line

KCWI = {'ang_res':??,'FOV':??,'bandpass':(??,??),'SBlim':??}
MUSE = {'ang_res':??,'FOV':??,'bandpass':(??,??),'SBlim':??}

line_energies = {'CIII_977':4.*np.pi*2.03e-11,
                 'CIV':4.*np.pi*1.28e-11,
                 'OVI':4.*np.pi*1.92e-11,
                 'SiIV':4.*np.pi*1.42e-11,
                 'HAlpha':4.*np.pi*3.03e-12,
                 'LyAlpha':4.*np.pi*1.63e-11}

class EmissionMap:
  def __init__(self,ds,track):
    self.ds = ds
    self.track = track
    self.rb,self.rb_center,self.rb_width = get_refine_box(self.ds,self.ds.current_redshift,self.track)
    self.redshift = ds.redshift
    self.box_width = ds.arr(rb_width,'code_length').in_units('kpc')
    return

    def make_KCWI_frb(line,fileout,index,cquan=None):
        field = 'Emission_'+line
        ## need to generate a box that matches the KCWI FOV
        ## something I need to decide is, do I put the galaxy in the center of the
        ## emission map? or do I  use the refine box center as the center?

        obj = ds.proj(('gas',field),index,data_source=boxKCWI)
        return

    def make_MUSE_frb(line,fileout,cquan=None):
      return

    def make_LLAMAS_frb(line,fileout,cquan=None):
      return

    def convert_photons_to_ergs(line):
      cquan = self.ds.quan(1.,'steradian**-1')
      return cquan*line_energy[line]
