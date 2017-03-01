import yt
import numpy as np
import trident
import matplotlib.pyplot as plt

def make_enzo_rays():
    fn = "/Volumes/sonic/Ryan/r0058_l10/redshift0058"
    ds = yt.load(fn,file_style="%s.grid.cpu%%04i")
    
    
## want to try Ryan's sim and perhaps Brian's to see if trident
## looks as bad with AMR outputs as it does with the SPH ones
    pos = [ 0.40328598,  0.47176743,  0.46131516] ## RYAN
    
    dxs = [ds.arr(5.,'kpc').in_units('code_length'),ds.arr(50.,'kpc').in_units('code_length')  ]
    names = ['Ryan','Ryan']
    dists = ['5kpc','50kpc']
    
    i = 0
    line_list = ['H I 1216']
    ds.unique_identifier = 'Ryan_z0'
    
    while i < len(dxs):
        dx = dxs[i]
        ray_start = [pos[0],0,pos[2]+dx.d] 
        ray_end   = [pos[0],1,pos[2]+dx.d]
        rayout = 'ray_'+names[i]+'_'+dists[i]+'.h5'    
        
        ray = trident.make_simple_ray(ds,ray_start,ray_end,lines=line_list,data_filename=rayout)
        sg = trident.SpectrumGenerator(lambda_min=1180.,lambda_max=1240.,dlambda=0.01)
        sg.make_spectrum(ray,lines=line_list)
        plt.plot(sg.lambda_field,sg.flux_field,label=dists[i])
        
        i = i + 1

    plt.legend()
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.savefig('Ryan_rays.png')
    return

## I'd also like to try and test how the H_number_density varies along
## a ray for yt2 and yt3 for both Ryan's and Brian's sims


def calculate_EW_arrays(ray,line):
    return


