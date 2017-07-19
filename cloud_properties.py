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

sp = ds.sphere(center,cloud_initial_radius)
cloud_initial_mass = sp.quantities['TotalMass']()[0].in_units('Msun')

DATA_DIR = "."
dataset_list = glob.glob(os.path.join(DATA_DIR,"DD000?/DD000?"))
i = 0


## Three Rays of Interest:
nrays = 3
len_dslist = len(dataset_list)
rays = {'center' : ([0.,1.,1.],[4.,1.,1.]),
	'1xradius'  : ([1.0625,0.,1.],[1.0625,2.,1.]),
	'2xradius'    : ([1.125,0.,1.],[1.125,2.,1.])
	}
ions = ['H','C','Si','O'] # ADD HI BACK!

#lines_exact = ['H I 1216','Si II 1260', 'Si III 1207', 'Si IV 1403', 
#                          'C II 1335',  'C III 977', 'C IV 1548',
#		'O IV 1032']


#Some quantities that I would want to store
cloud_mass_dens = []
cloud_mass_temp = []
times = []

coldens = {'H I'    : ('H','0',1216.,np.zeros((len_dslist,nrays))),
	   'C II'   : ('C','1',1335.,np.zeros((len_dslist,nrays))),
	   'C III'  : ('C','2',977.,np.zeros((len_dslist,nrays))),
	   'C IV'   : ('C','3',1548.,np.zeros((len_dslist,nrays))),
	   'Si II'  : ('Si','1',1260.,np.zeros((len_dslist,nrays))),
	   'Si III' : ('Si','2',1207.,np.zeros((len_dslist,nrays))),
	   'Si IV'  : ('Si','3',1403.,np.zeros((len_dslist,nrays))),
	   'O VI'   : ('O','5',1032.,np.zeros((len_dslist,nrays)))
	  }


for dataset in dataset_list:
    print 'i is ',i
    ds = yt.load(dataset)
    ad = ds.all_data()
    times.append(ds.current_time.in_units('Myr'))
    idx = np.where(ad['Density'].in_units('code_mass/code_length**3') > (1./3.)*3000.)[0]
    cloud_mass_dens.append(np.sum(ad['cell_mass'][idx].in_units('Msun')))
    trident.add_ion_fields(ds,ions=ions)

    if make_rays:
       j = 0
       if make_plots:
          sp = yt.SlicePlot(ds,'z','Density')
	  sp.set_unit('Density','g/cm**3')
	  #sp.set_cmap('Density','YlGnBu')
       for ray,(start,end) in rays.items():
		ray_in = trident.make_simple_ray(ds,start_position=start,end_position=end,
						lines=ions,ftype='gas',
						data_filename=ds.basename+'_'+ray+'.h5')
		ad_ray = ray_in.all_data()
		for ion,(element,number,lambda_rest,Narr) in coldens.items():
			field_in = element+'_p'+number+'_number_density'
			if ((element=='H') & (number=='0')):
				field_in = 'H_number_density'
			Nhere = np.log10(np.sum(ad_ray[('gas',field_in)]*ad_ray['dl']))
			Narr[i,j] = Nhere
		
##			MAKE SOME SPECTRA AND PLOTS!
			if make_rays:
				lambda_min = lambda_rest - 5.
                	        lambda_max = lambda_rest + 5.
                	        line_name = ion+' '+str(int(lambda_rest))
			
				sg = trident.SpectrumGenerator(lambda_min=lambda_min, \
        				lambda_max=lambda_max, dlambda=0.0001)
    				sg.make_spectrum(ray_in,lines=line_name)

				## Make the plots (they should all have the same axes...)
				plt.plot(sg.lambda_field,sg.flux_field)
				plt.xlabel('Wavelength [A]')
				plt.ylabel('Relative Flux')
				plt.ylim(0.99,1.)
				plt.title(line_name)
				plt_out = ds.basename+'_'+ray+'_'+element+str(int(number)+1)+'.png'
				plt.savefig(plt_out)
				plt.close()

		if make_plots:
			print 'Add Ray!!!!!!!!!!!!!!!!!!!!!'
			print start,end
			aray = ds.ray(start,end)
			sp.annotate_ray(aray,arrow=True)
		j = j + 1
    if make_rays:
	sp.save()
    i = i + 1

cPickle.dump(coldens,open('coldens_rays_HOLD.cpkl','wb'),protocol=-1)

if make_plots:
	plt.plot(times,cloud_mass_dens/cloud_initial_mass,'o',label='dens')
	plt.legend()
	plt.xlabel('Time [Myr]')
	plt.ylabel('Cloud Mass / Initial Cloud Mass')
	plt.savefig('cloud_mass_evolution.png')
	plt.close()

	for ion,(element,number,lambda_rest,Narr) in coldens.items():
		for i in range(nrays):
			plt.plot(times,Narr[:,i],'o',label=rays.items()[0][0])
		plt.legend()
		plt.xlabel('Time [Myr]')
		plt.ylabel('log(Column Density) [cm^-3]')
		plt.savefig(ion+'_coldens_evolve.png')
		plt.close()

cPickle.dump(coldens,open('coldens_rays.cpkl','wb'),protocol=-1)

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
