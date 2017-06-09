import yt
import trident
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import cPickle

def calculate_mass_density_threshold(ad,threshold,overdensity):
    idx = np.where(ad['Density'] > threshold*overdensity)[0]
    return idx, np.sum(ad['cell_mass'][idx].in_units('Msun'))

def calculate_mass_temperature_threshold(ad,threshold):
    idx = np.where(ad['Temperature'] < threshold)[0]
    return idx, np.sum(ad['cell_mass'][idx].in_units('Msun'))

def compute_cloud_crushing_time(cloud_initial_radius,initial_overdensity,wind_velocity):
    rad = cloud_initial_radius.in_cgs()
    wind_velocity = wind_velocity.in_cgs()
    return (rad*(initial_overdensity**0.5))/wind_velocity


ds = yt.load('DD0000/DD0000')
ad = ds.all_data()


cloud_initial_radius = ds.arr(0.0625,'code_length')
initial_overdensity = 3000.
wind_velocity = ds.arr(1e7, 'km/s')

## Some basic quantities that I would want to know!
#ambient_density
#cloud_density
#shock_metallicity
#ambient_metallicity
#cloud_metallicity
#shock_metallicity
#ambient_temperature
#cloud_temperature
#shock_temperature
#wind_velocity
#cloud_initial_radius
#cloud_initial_mass

#Some quantities that I would want to store
#cloud_metal_dens = []
#cloud_metal_temp = []
#cloud_size = []
cloud_mass_dens = []
cloud_mass_temp = []
times = []

idx = np.where(ad['Temperature'] <= 10.**4.0)[0]
cloud_initial_mass = np.sum(ad['cell_mass'][idx].in_units('Msun'))

DATA_DIR = "."
dataset_list = glob.glob(os.path.join(DATA_DIR,"DD????/DD????"))
i = 0

for dataset in dataset_list:
    print 'i is ',i
    ds = yt.load(dataset)
    ad = ds.all_data()
    times.append(ds.current_time.in_units('Myr'))
    idx = np.where(ad['Density'].in_units('code_mass/code_length**3') > (1./3.)*3000.)[0]
    cloud_mass_dens.append(np.sum(ad['cell_mass'][idx].in_units('Msun')))
    idx = np.where(ad['Temperature'].in_units('K') < 10.**4.)[0]
    cloud_mass_temp.append(np.sum(ad['cell_mass'][idx].in_units('Msun')))
    i = i + 1


plt.plot(times,cloud_mass_dens/cloud_initial_mass,lw=2.0,label='dens')
plt.plot(times,cloud_mass_temp/cloud_initial_mass,lw=2.0,label='temp')
plt.legend()
plt.xlabel('Time [Myr]')
plt.ylabel('Cloud Mass / Initial Cloud Mass')
plt.savefig('cloud_mass_evolution_T4.png')
plt.close()

