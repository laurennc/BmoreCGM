import yt
import numpy as np
import matplotlib.pyplot as plt
import cPickle
import trident
from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
from emission_functions import *
import yt.units as u

# ok so general structure is that I'm going to have to make an FRB
## for each instrument for each line

### !!!!!!!!!!!!!!!!!!!!!!!!!NOTE!!!!!!!!!!!!!!!!!!!!!!!!! ###
### ONE OF THE DECISIONS THAT I'VE MADE IS THAT YOU'RE ALWAYS PROJECTING
### THROUGH THE REFINE BOX AND NOTHING ELSE
### NOT SURE IF THAT'S REALLY WHAT I ALWAYS WANT TO BE DOING

## SB numbers from https://arxiv.org/pdf/1803.10781.pdf = 1e-18
## Can confirm with the LU quoted in the KCWI document
## but that gives 6e-19 (6e3 LU), and 7e-21 (70 LU) respectively
## LU = 1e5 ph cm^-2 s^01 sr^-1 = 10^-17 erg cm^-2 s^-1 arcsec^-2 for LyA
## My guess is that the integration times are different so let's go with
## the lower SB limits for now
KCWI = {'binned2x2' :{'ang_res':1.*u.arcsecond,'FOV':(33*u.arcsecond,20.4*u.arcsecond),'fig_size':(2.2,1.36),
                      'bandpass':(3500.*u.angstrom,5600.*u.angstrom),'SBlim':6e-19*u.erg/(u.s*u.arcsecond**2*u.cm**2)},
        'full_slice':{'ang_res':0.5*u.arcsecond,'FOV':(33*u.arcsecond,20.4*u.arcsecond),'fig_size':(2.2,1.36),
                      'bandpass':(3500.*u.angstrom,5600.*u.angstrom),'SBlim':7e-21*u.erg/(u.s*u.arcsecond**2*u.cm**2)}
        }
## Getting these SB numbers from https://arxiv.org/pdf/1509.05143.pdf
## where they have 27 hours of exposure so they go much deeper than
## the one hour exposure quoted on the MUSE spec website
MUSE = {'wide'  :{'ang_res':0.2*u.arcsecond,'FOV':(1.*u.arcminute,1.*u.arcminute),'fig_size':(4,4),
                  'bandpass':(4650.*u.angstrom,9300.*u.angstrom),'SBlim':1e-19*u.erg/(u.s*u.arcsecond**2*u.cm**2)},
        'narrow':{'ang_res':0.025*u.arcsecond,'FOV':(7.5*u.arcsecond,7.5*u.arcsecond),'fig_size':(0.5,0.5),
                  'bandpass':(4650.*u.angstrom,9300.*u.angstrom),'SBlim':1e-19*u.erg/(u.s*u.arcsecond**2*u.cm**2)}
        }
## Specs from http://www.mit.edu/people/rsimcoe/Llamas_pocketguide_rA.pdf
## There's nothing quoted about the SB but let's assume that they're going
## to do as well as MUSE or KCWI or it wouldn't get funded
LLAMAS = {'normal' :{'ang_res':0.75*u.arcsecond,'FOV':(1.*u.arcminute,1.*u.arcminute),'fig_size':(4,4),
          'bandpass':(3600.*u.angstrom,9700.*u.angstrom)),'SBlim':1e-19*u.erg/(u.s*u.arcsecond**2*u.cm**2)}
          }

line_energies = {'CIII_977':2.03e-11*u.erg,'CIV':1.28e-11*u.erg,
                 'OVI':1.92e-11*u.erg,'SiIV':1.42e-11*u.erg,
                 'HAlpha':3.03e-12*u.erg,'LyAlpha':1.63e-11*u.erg}

line_wavelengths = {'CIII_977':977*u.angstrom,'CIV':??*u.angstrom,
                    'OVI':1032*u.angstrom,'SiIV':??*u.angstrom,
                    'HAlpha':??*u.angstrom,'LyAlpha':1216*u.angstrom}

axis_map = {'x':0,'y':1,'z':2}

class EmissionMap:
  def __init__(self,ds,track,instrument,mode):
    self.ds = ds
    self.track = track
    self.rb,self.rb_center,self.rb_width = get_refine_box(self.ds,self.ds.current_redshift,self.track)
    self.redshift = ds.redshift
    self.box_width = ds.arr(rb_width,'code_length').in_units('kpc')
    self.instrument = instrument
    self.mode = mode
    return

  def make_frb(self,line,proj_axis,pickle_out):
      field = "Emission_"+line
      if isinstance(axis,string):
          proj_axis = axis_map[axis]

      area_axes = [0,1,2]
      i = area_axes.index(proj_axis)
      del area_axes[i]

      FOVconvert = cosmo.kpc_proper_per_arcmin(self.redshift).value*u.kpc/u.arcmin
      FOV = (self.instrument[self.mode]['FOV']*FOVconvert).to('kpc')
      num_cells = (np.ceil(self.instrument[self.mode]['FOV']/self.instrument[self.mode]['ang_res'])).value

      ## Now set up projection region
      ## Area should be size of FOV and not the refine region scaled down
      FOV = ds.arr(FOV,'kpc').to('code_length')
      LE = np.copy(self.rb_center)
      LE[area_axes] -= FOV/2.
      LE[proj_axis] -= rb_width/2.
      RE = np.copy(self.rb_center)
      RE[area_axes] += FOV/2.
      RE[proj_axis] += rb_width/2.

      box = ds.r[LE[0]:RE[0],LE[1]:RE[1],LE[2]:RE[2]]
      obj = ds.proj(('gas',field),proj_axis,data_source=box)
      frb = obj.to_frb((FOV[0].value,'code_length'),(int(num_cells[0]),int(num_cells[1])),
                        height=(FOV[1].value,'code_length'),center=rb_center)
      cPickle.dump(frb[('gas',field)],open(pickle_out,'wb'),protocol=-1)
      return frb


    def plot_frb(self,line,frb,fileout):
        if check_line_in_bandpass(line):
            ## NEED TO SET PLOT SIZE SO THAT IT'S CONSISTENT
            ## BETWEEN ALL OF THE INSTRUMENTS
            plt.figure(figsize=self.instrument[self.mode]['fig_size'])
            frb = cPickle.load(open(frb,'rb'))

            ## Convert to the units most observers use
            lenergy = line_energies[line]
            frb = frb.to(1./((u.s*u.arcsecond**2*u.cm**2)))*lenergy
            frb = np.log10(frbNAT/(1.+self.redshift)**4)

            FOVconvert = cosmo.kpc_proper_per_arcmin(self.redshift).value*u.kpc/u.arcmin
            FOV = (self.instrument[self.mode]['FOV']*FOVconvert).to('kpc')
            bd1,bd2 = -1*FOV[0]/2.,FOV[0]/2.
            bd3,bd4 = -1*FOV[1]/2./FOV[1]/2.

            plt.imshow(frb,cmap=cmap,vmin=-24,vmax=-14,extent=(bd1,bd2,bd3,bd4),
                       origin='lower',interpolation=None)
            plt.colorbar()
            plt.savefig(fileout)
            plt.close()
        else:
            PLOT A BLANK SQUARE
        return

        return

    def check_line_in_bandpass(self,line):
        if isinstance(line,str):
            line = line_wavelengths[line]
        line_observed = (1+self.redshift)*line
        bandpass = self.instrument[self.mode]['bandpass']
        if ((line_observed > bandpass[0]) & (line_observed < bandpass[1])):
            return True
        else
            return False

def emission_line_energy(wavelength,unit):
    E = h*c/(wavelength*unit)
    return E.to('erg')
