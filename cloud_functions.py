import yt
import numpy as np
import matplotlib.pyplot as plt

def starting_cloud_region(ds):
    pos = [1,1,1]
    rad = 0.0625*1.543e21 ## cm    
    return ds.sphere(pos,(rad,'cm'))


def find_GasMass_dans_TempRange(region,low=5.0,high=1e7):
    idx = np.where( (region['temperature'] <= high) & (region['temperature'] >= low) )[0]
    return np.sum(region['cell_mass'][idx]).in_units('Msun')


def plot_ray_properties(rays,filebase):
    fields = ['density',('gas','H_number_density'),'metallicity','temperature']
    ads = []
    for ray in rays:
        ads = np.append(ads,ray.all_data() )

    #names = ['lowres','metalcool','maxref5','metalcoolHM12','medres']
    names = ['AMR','Unigrid']

    for field in fields:
        i = 0
        if type(field) == tuple:
            fieldout = field[1]
        else:
            fieldout = field
        for ray in ads:
            strout = 'Ray '+str(i)
            idx = np.where(np.log10(ray['x'].in_units('pc')) > 2.5  )
            plt.plot(np.log10(ray['x'].in_units('pc'))[idx],np.log10(ray[field][idx])
                ,label=names[i])           
            i = i + 1
        plt.xlabel('x [pc]')
        plt.ylabel(field)
        plt.legend()
        plt.savefig(filebase+'_'+fieldout+'.png')
        plt.close()
    return 'made plots!'

def load_rays(ray_list):
    rays = []
    for line in ray_list:
        ray = yt.load(line)
        rays = np.append(rays,ray)
    
    return rays

