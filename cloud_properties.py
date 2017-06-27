import yt
import trident
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import cPickle

def calculate_mass_density_threshold(ad,threshold,overdensity,dens_ambient):
    idx = np.where(ad['Density'].in_units('code_mass/code_length**3')
                       > threshold*overdensity*dens_ambient)[0]
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
center = ds.arr([1,1,1],'code_length')

make_plots = True
make_rays  = True

## Three Rays of Interest:
rays = {'center' : ([0.,1.,1.],[4.,1.,1.]),
	'1xradius'  : ([1.0625,0.,1.],[1.0625,2.,1.]),
	'2xradius'    : ([1.125,0.,1.],[1.125,2.,1.])
	}
ions = ['H','C','Si','O'] # ADD HI BACK!


#Some quantities that I would want to store
cloud_mass_dens = []
cloud_mass_temp = []
times = []

coldens = {'H I'    : ('H','0',[]),
	   'C III'  : ('C','2',[]),
	   'C IV'   : ('C','3',[]),
	   'Si II'  : ('Si','1',[]),
	   'Si III' : ('Si','2',[]),
	   'Si IV'  : ('Si','3',[]),
	   'O VI'   : ('O','5',[])
	  }

#line_centers = { ??? }

sp = ds.sphere(center,cloud_initial_radius)
cloud_initial_mass = sp.quantities['TotalMass']()[0].in_units('Msun')

DATA_DIR = "."
dataset_list = glob.glob(os.path.join(DATA_DIR,"DD000?/DD000?"))
i = 0

for dataset in dataset_list:
    print 'i is ',i
    ds = yt.load(dataset)
    ad = ds.all_data()
    times.append(ds.current_time.in_units('Myr'))
    idx = np.where(ad['Density'].in_units('code_mass/code_length**3') > (1./3.)*3000.)[0]
    cloud_mass_dens.append(np.sum(ad['cell_mass'][idx].in_units('Msun')))
    trident.add_ion_fields(ds,ions=ions)

    if make_rays:
       if make_plots:
          sp = yt.SlicePlot(ds,'z','Density')
	  sp.set_cmap('Density','YlGnBu')
       for ray,(start,end) in rays.items():
		ray_in = trident.make_simple_ray(ds,start_position=start,end_position=end,
						lines=ions,ftype='gas',
						data_filename=ds.basename+'_'+ray+'.h5')
		ad_ray = ray_in.all_data()
		for ion,(element,number,Narr) in coldens.items():
			field_in = element+'_p'+number+'_number_density'
			if ((element=='H') & (number=='0')):
				field_in = 'H_number_density'
			Narr.append(np.log10(np.sum(ad_ray[('gas',field_in)]*ad_ray['dl'] )).item())

			
			#line_center = ??
			#sg = trident.SpectrumGenerator(??????)
			#sg.make_spectrum(ray_in,lines=[ion])

		if make_plots:
			print 'ADD RAY!!!'
			aray = ds.ray(start,end)
			sp.annotate_ray(aray,arrow=True)

    if make_rays:
	sp.save()
    i = i + 1


if make_plots:
	plt.plot(times,cloud_mass_dens/cloud_initial_mass,'o',label='dens')
	plt.legend()
	plt.xlabel('Time [Myr]')
	plt.ylabel('Cloud Mass / Initial Cloud Mass')
	plt.savefig('cloud_mass_evolution.png')
	plt.close()

	for ion,(element,number,Narr) in coldens.items():
		plt.plot(times,Narr,'o',label='ion')
	plt.legend()
	plt.xlabel('Time [Myr]')
	plt.ylabel('log(Column Density) [cm^-3]')
	plt.savefig('coldens_evolve.png')
	plt.close()


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
