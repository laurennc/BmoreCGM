import matplotlib
matplotlib.use('Agg')

import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.analysis_modules.star_analysis.api import StarFormationRate
import os
from radial_data_nozeros import *
import trident
import cPickle

plot_kwargs = {
    'nref10_track_2'         : {'ls':'-','color':'#e7298a'}, #,'marker':'o','markersize':'0.25'},
    'nref10_track_lowfdbk_1' : {'ls':'-','color':'#d95f02'}, #,'marker':'s','markersize':'0.25'},
    'nref10_track_lowfdbk_2' : {'ls':'-','color':'#1b9e77'}, #,'marker':'*','markersize':'0.25'},
    'nref10_track_lowfdbk_3' : {'ls':'-','color':'#7570b3'}, #,'marker':'^','markersize':'0.25'},
    'nref10_track_lowfdbk_4' : {'ls':'-','color':'#e6ab02'}, #,'marker':'p','markersize':'0.25'},
    'nref10_z1_0.5_natural'           : {'ls':'--','color':'#e7298a'}, #,'marker':'o','markersize':'0.25'},
    'nref10_z1_0.5_natural_lowfdbk_1' : {'ls':'--','color':'#d95f02'}, #,'marker':'s','markersize':'0.25'},
    'nref10_z1_0.5_natural_lowfdbk_2' : {'ls':'--','color':'#1b9e77'}, #,'marker':'*','markersize':'0.25'},
    'nref10_z1_0.5_natural_lowfdbk_3' : {'ls':'--','color':'#7570b3'}, #,'marker':'^','markersize':'0.25'},
    'nref10_z1_0.5_natural_lowfdbk_4' : {'ls':'--','color':'#e6ab02'}, #,'marker':'p','markersize':'0.25'}
}

### utility functions ###
def get_halo_center(ds, center_guess):
    proper_box_size = ds.get_parameter('CosmologyComovingBoxSize') / ds.get_parameter('CosmologyHubbleConstantNow') * 1000. # in kpc
    ad = ds.sphere(center_guess, (200., 'kpc'))
    x,y,z = np.array(ad["x"]), np.array(ad["y"]), np.array(ad["z"])
    dm_density =  ad['Dark_Matter_Density']
    imax = (np.where(dm_density > 0.9999 * np.max(dm_density)))[0]
    halo_center = [x[imax[0]], y[imax[0]], z[imax[0]]]
    #print 'We have located the main halo at :', halo_center
    return halo_center

def diskVectors(ds, center):
    sphere = ds.sphere(center,(5.,'kpc'))
    angular_momentum = sphere.quantities['AngularMomentumVector']()
    x = np.cross(angular_momentum,[1.0,0.0,0.0])
    x /= np.linalg.norm(x)
    return (angular_momentum,x)

def compute_disk_masses(filenames):
    rds = np.arange(27,43)
    rds = rds[::-1]
    timesteps = np.zeros(len(rds))
    gas_masses = np.zeros((len(rds),len(filenames)))
    stellar_masses = np.zeros((len(rds),len(filenames)))

    for j in range(len(filenames)):
        center_guess = [0.48988587,0.47121728,0.50938220]
        for i in range(len(rds)):
            filein = filenames[i]+'/RD00'+str(rds[i])+'/RD00'+str(rds[i])
            ds = yt.load(filein)
            halo_center = get_halo_center(ds,center_guess)
            center_guess = halo_center
            (Lx,x) = diskVectors(ds,halo_center)
            disk = ds.disk(halo_center,Lx,(40.,'kpc'),(20.,'kpc'))
            ## let's look for cold/dense gas
            ## ISM conditions from review paper:
            ## ??
            idx = np.where(() & ())[0]
        i = i + 1
    return disk_masses

def make_frbs(filename,center,fields,ions,fdbk=False):
    ds = yt.load(filename)
    args = filename.split('/')

    trident.add_ion_fields(ds,ions=ions)

    if fdbk== True:
        print 'in fdbk'
        halo_center = get_halo_center(ds,center)

        box_center = np.copy(halo_center)
        box_center[1] = box_center[1]+ds.arr(60.,'kpc').in_units('code_length').value

        dx = ds.arr(40.,'kpc').in_units('code_length').value
        dy = ds.arr(80.,'kpc').in_units('code_length').value
        box_left  = [box_center[0]-dx, box_center[1]-dy, box_center[2]-dx]
        box_right = [box_center[0]+dx, box_center[1]+dy, box_center[2]+dx]

        refine_box = ds.r[box_left[0]:box_right[0], box_left[1]:box_right[1], box_left[2]:box_right[2]]

        width = [(160,'kpc'),(80.,'kpc')]
        resolution = (320,160)
        for field in fields:
            fileout = args[-3]+'_'+args[-2]+'_x_'+field+'.cpkl'
            obj = ds.proj(field,'x',data_source=refine_box)
            frb = obj.to_frb(width,resolution,center=box_center)
            cPickle.dump(frb[field],open(fileout,'wb'),protocol=-1)

    else:
        print 'in else'
        halo_center = get_halo_center(ds,center)
        box_center = np.copy(halo_center)
        dx = ds.arr(60.,'kpc').in_units('code_length').value
        box_left  = [box_center[0]-dx, box_center[1]-dx, box_center[2]-dx]
        box_right = [box_center[0]+dx, box_center[1]+dx, box_center[2]+dx]

        refine_box = ds.r[box_left[0]:box_right[0], box_left[1]:box_right[1],box_left[2]:box_right[2]]
        width = [(120,'kpc'),(120.,'kpc')]
        resolution = (240,240)

        for field in fields:
            fileout = args[-3]+'_'+args[-2]+'_x_'+field+'.cpkl'
            obj = ds.proj(field,'x',data_source=refine_box)
            frb = obj.to_frb(width,resolution,center=box_center)
            cPickle.dump(frb[field],open(fileout,'wb'),protocol=-1)
    return

def fdbk_coldens_profile(frb):
    whole_box = np.zeros((560,560))
    dx = 80 #40/0.5
    dy1 = 40 #20/0.5
    dy2 = 280 #140/0.5

    #frb = frb.T
    whole_box[280-80:280+80,280-40:280+280] = np.log10(frb)
    mask = np.zeros((560,560),dtype=bool)
    mask[280-80:280+80,280-40:280+280] = True

    xL = np.linspace(-140,140,560)
    xL,yL = np.meshgrid(xL,xL)
    rp = radial_data(whole_box,working_mask=mask,x=xL,y=yL)
    return rp

def get_evolultion_Lx(filenames,center_guess):
    rds = np.arange(27,43)
    rds = rds[::-1]
    timesteps = np.zeros(len(rds))
    lvectors = np.zeros((len(rds),len(filenames)))

    for j in range(len(filenames)):
        for i in range(len(rds)):
            filein = filenames[i]+'/RD00'+str(rds[i])+'/RD00'+str(rds[i])
            ds = yt.load(filein)
            halo_center = get_halo_center(ds,center_guess)
            (Lx,x) = diskVectors(ds,halo_center)
            lvectors[j,i] = Lx

            if (j==0) & (i == 0):
                timesteps[i] = ds.current_redshift

    return

### plotting functions ###
def plot_SFHS(filenames,center,radius,pltname):
    fig,ax = plt.subplots(6,1,sharex=True,sharey=True)
    fig.set_size_inches(6,14)
    fig.subplots_adjust(hspace=0.1,wspace=0.1)
    i = 0
    for i in range(len(filenames)):
        ds = yt.load(filenames[i])
        sp = ds.sphere(center,radius)
        sfr = StarFormationRate(ds, data_source=sp)
        args = filenames[i].split('/')[-3]
        kwargs = plot_kwargs[args]
        if i < 5:
            ax[0].plot(sfr.redshift,sfr.Msol_yr,**kwargs)
            ax[i+1].plot(sfr.redshift,sfr.Msol_yr,**kwargs)
        if i > 4:
            ax[i-4].plot(sfr.redshift,sfr.Msol_yr,**kwargs)
            ax[i-4].set_xlim(1,0)
            ax[i-4].set_ylim(0,25)
            #ax[i-3].annotate('')
        i = i + 1
    ax[0].set_xlim(1,0)
    ax[0].set_ylim(0,25)
    ax[5].set_xlabel('Redshift')
    ax[2].set_ylabel('SFR [Msun/yr]')
    plt.savefig(pltname)
    return

def confirm_halo_centers(filenames,center):
    for i in range(len(filenames)):
        ds = yt.load(filenames[i])
        args = filenames[i].split('/')[-3]
        halo_center = get_halo_center(ds,center)
        sl = yt.SlicePlot(ds,'x','Density',center=halo_center,width=(200,'kpc'))
        sl.annotate_text(center,'c')
        sl.save(args)
    return

def plot_coldens_radialprofiles(filenames,fields,xlen,ylen,fileout,fdbk=False):
    if fdbk == True:
        fig,ax = plt.subplots(6,len(fields),sharex=True) #,sharey=True)
        fig.set_size_inches(xlen,ylen)
        fig.subplots_adjust(hspace=0.1,wspace=0.1)
        i = 0
        for i in range(len(filenames)):
            args = filenames[i].split('/')
            kwargs = plot_kwargs[args[-3]]
            for j in range(len(fields)):
                filein = 'coldens_cpkl/'+args[-3]+'_'+args[-2]+'_x_'+fields[j]+'.cpkl'
                #print filein
                frb = cPickle.load(open(filein,'rb'))
                rp = fdbk_coldens_profile(frb)
                if i < 5:
                    if len(fields) == 1:
			ax[0].plot(rp.r,rp.mean,**kwargs)
		    	ax[i+1].plot(rp.r,rp.mean,**kwargs)
		    else:
		    	ax[0,j].plot(rp.r,rp.mean,**kwargs)
                    	ax[i+1,j].plot(rp.r,rp.mean,**kwargs)
                if i > 4:
		    if len(fields) == 1:
		    	ax[0].set_title(fields[j].split('_')[0:2])
		        ax[i-4].plot(rp.r,rp.mean,**kwargs)
			ax[5].set_xlabel('Radius [kpc]')	
		    else:
                    	ax[0,j].set_title(fields[j].split('_')[0:2])
                    	ax[i-4,j].plot(rp.r,rp.mean,**kwargs)
        		ax[5,j].set_xlabel('Radius [kpc]')
                    #ax[i-4].set_xlim(1,0)
                    #ax[i-4].set_ylim(0,25)
        plt.savefig(fileout)
    return

###################################################################################################

filenames = ['/astro/simulations/FOGGIE/halo_008508/nref10_track_2/RD0042/RD0042',
             '/astro/simulations/FOGGIE/halo_008508/nref10_track_lowfdbk_1/RD0042/RD0042',
             '/astro/simulations/FOGGIE/halo_008508/nref10_track_lowfdbk_2/RD0042/RD0042',
             '/astro/simulations/FOGGIE/halo_008508/nref10_track_lowfdbk_3/RD0042/RD0042',
             '/astro/simulations/FOGGIE/halo_008508/nref10_track_lowfdbk_4/RD0042/RD0042',
             '/astro/simulations/FOGGIE/halo_008508/nref10_z1_0.5_natural/RD0042/RD0042',
             '/astro/simulations/FOGGIE/halo_008508/nref10_z1_0.5_natural_lowfdbk_1/RD0042/RD0042',
             '/astro/simulations/FOGGIE/halo_008508/nref10_z1_0.5_natural_lowfdbk_2/RD0042/RD0042',
             '/astro/simulations/FOGGIE/halo_008508/nref10_z1_0.5_natural_lowfdbk_3/RD0042/RD0042',
             '/astro/simulations/FOGGIE/halo_008508/nref10_z1_0.5_natural_lowfdbk_4/RD0042/RD0042']

#filenames = ['/Users/dalek/data/Jason/nref10_track_lowfdbk_1/RD0042/RD0042']
#ds = yt.load(filenames[0])
#center = [0.48988587,0.47121728,0.50938220]
#halo_center = get_halo_center(ds,center)
halo_center = np.array([0.48984, 0.47133, 0.50956])

fields = ['H_p0_number_density','O_p5_number_density','C_p3_number_density','C_p2_number_density',
          'Si_p2_number_density','Si_p3_number_density']

#plot_coldens_radialprofiles(filenames,fields,fdbk=True)

#ions = ['O VI','C IV','C III','Si III','Si IV']
#for filename in filenames:
#    make_frbs(filename,halo_center,fields,ions,fdbk=True)


Cfields  = ['C_p3_number_density','C_p2_number_density']
SiFields = ['Si_p2_number_density','Si_p3_number_density']
OFields  = ['O_p5_number_density']

#plot_coldens_radialprofiles(filenames,['H_p0_number_density'],6,14,
#			    'coldensH.pdf',fdbk=True)
#plot_coldens_radialprofiles(filenames,Cfields,12,14,'coldensC.pdf',fdbk=True)
#plot_coldens_radialprofiles(filenames,SiFields,12,14,'coldensSi.pdf',fdbk=True)
#plot_coldens_radialprofiles(filenames,OFields,6,14,'coldensO.pdf',fdbk=True)
plot_coldens_radialprofiles(filenames,fields,24,14,'coldensALL.pdf',fdbk=True)

#plot_SFHS(filenames,halo_center,(300.,'kpc'),'feedback_SFHs.pdf')
#confirm_halo_centers(filenames,halo_center)



