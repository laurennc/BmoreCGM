import yt
import numpy as np
import datetime
from astropy.io import fits


fn = "/Volumes/sonic/Ryan/r0043/redshift0043"
ds = yt.load(fn,file_style="%s.grid.cpu%%04i")
#val, pos = ds.find_max('Density')
pos = ds.arr([ 0.39871597,  0.46913528,  0.46808243], 'code_length')
width = ds.arr(320,'kpc').in_units('code_length')
dx = ds.arr(0.5,'kpc').in_units('code_length')
num_cells = 1200.

xL,xR = pos[0]-width,pos[0]+width
yL,yR = pos[1]-width,pos[1]+width
zL,zR = pos[2]-width,pos[2]+width
box = ds.r[xL:xR:640j,yL:yR:640j,zL:zR:640j]

primhdr = fits.Header()
primhdr['SIMNAME'] = 'Joung2012'
primhdr['REDSHIFT'] = 0.2
primhdr['DATE'] = datetime.datetime.now().isoformat()


hdr_new = fits.HDUList([ fits.PrimaryHDU(header=primhdr),
                   fits.ImageHDU(data=box['density'], name='density'),
                   fits.ImageHDU(data=box['temperature'], name='temperature'), #Flip PC2
                   fits.ImageHDU(data=box['metallicity'], name='metallicity'),
                   fits.ImageHDU(data=box['H_nuclei_density'], name='Hdensity')])

hdr_new.writeto('/Users/dalek/Desktop/simcube_z02.fits', overwrite=True)
hdr_new.close()
