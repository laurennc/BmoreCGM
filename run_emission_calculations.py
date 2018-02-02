import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import cPickle
from emission_functions import *
import brewer2mpl as brew
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from astropy.cosmology import WMAP9 as cosmo
from matplotlib import colors

#base = "/Users/dalek/data/Molly/natural/nref11"
base = "/Users/dalek/data/Molly/nref11n_nref10f_refine200kpc_z4to2"
fn = base+"/RD0020/RD0020"
lines = ['HAlpha','OVI','CIV','CIII_977','SiIV']
track_name = base+"/halo_track"
args = fn.split('/')

bmap = brew.get_map('PuRd','Sequential',9)
pink_cmap = bmap.get_mpl_colormap(N=1000, gamma=2.0)

fontcs ={'fontname':'Helvetica','fontsize':16}
mpl.rc('text', usetex=True)
ds = yt.load(fn)
track = Table.read(track_name, format='ascii')
track.sort('col1')
rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)
redshift = ds.current_redshift

box_width = ds.arr(rb_width,'code_length').in_units('kpc')

def create_emission_frbs():
    ds = yt.load(fn)
    track = Table.read(track_name, format='ascii')
    track.sort('col1')
    rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)
    dx = np.unique(rb['dx']).max()
    dx = ds.arr(0.1829591205,'kpc').in_units('code_length')
    print 'Natural Refinement Resolution: ',dx.in_units('kpc')

    num_cells = np.ceil(rb_width/dx)
    #res_list = [0.5,1,5,10]

    for line in lines:
        field = 'Emission_'+line
        print line
        for index in 'xyz':
                print index
            #for res in res_list:
                #print res,' kpc'
                #dx = ds.arr(res,'kpc').in_units('code_length').value
                #num_cells = np.ceil(rb_width/dx)
                fileout = args[-3]+'_'+args[-2]+'_'+field+'_'+index+'_refwidth'
                proj = yt.ProjectionPlot(ds,index,('gas',field),width=(rb_width,'code_length'),
                                         center=rb_center,data_source=rb)
                proj.set_cmap(('gas',field),pink_cmap)
                proj.save(fileout+'.pdf')

                proj = yt.ProjectionPlot(ds,index,'density',width=(rb_width,'code_length'),
                                         center=rb_center,data_source=rb)
                proj.save(args[-3]+'_'+args[-2]+'_density_'+index+'.pdf')

                obj = ds.proj(('gas',field),index,data_source=rb)
                frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
                #fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_'+field+'_'+str(res)+'kpc.cpkl'
                fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_'+field+'_forcedres.cpkl'
                cPickle.dump(frb[('gas',field)],open(fileout,'wb'),protocol=-1)
    return

def create_phys_emis_weight_frbs():
    ds = yt.load(fn)
    track = Table.read(track_name, format='ascii')
    track.sort('col1')
    rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)
    #dx = np.unique(rb['dx']).max()
    dx = ds.arr(0.1829591205,'kpc').in_units('code_length')
    #print 'Natural Refinement Resolution: ',dx.in_units('kpc')
    num_cells = np.ceil(rb_width/dx)
    #res_list = [0.5,1,5,10]
    lines = ['HAlpha','CIII_977','CIV','SiIV','OVI']

    for line in lines:
        field = 'Emission_'+line
        print line

        for index in 'xyz':
            print index

            #obj = ds.proj('H_nuclei_density',index,data_source=rb,weight_field='H_nuclei_density')
            obj = ds.proj('H_nuclei_density',index,data_source=rb,weight_field=field)
            frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
            #fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_'+'hden_forcedres.cpkl'
            #print fileout

            #fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_hden_'+field+'_'+str(res)+'kpc.cpkl'
            fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_hden_'+field+'_forcedres.cpkl'
            cPickle.dump(frb['H_nuclei_density'],open(fileout,'wb'),protocol=-1)


            obj = ds.proj('temperature',index,data_source=rb,weight_field=field)
            frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
            #fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_temp_Whden_forcedres.cpkl'
            #print fileout

            #fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_temp_'+field+'_'+str(res)+'kpc.cpkl'
            fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_temp_'+field+'_forcedres.cpkl'
            cPickle.dump(frb['temperature'],open(fileout,'wb'),protocol=-1)

    return
	#for res in res_list:
    #	print res
        #dx = ds.arr(res,'kpc').in_units('code_length').value
        #num_cells = np.ceil(rb_width/dx)

def plot_ytProjections():
    fileout = args[-3]+'_'+args[-2]+'_'+field+'_'+index+'_refwidth'
    proj = yt.ProjectionPlot(ds,index,('gas',field),width=(rb_width,'code_length'),
                             center=rb_center,data_source=rb)
    proj.set_cmap(('gas',field),pink_cmap)
    proj.save(fileout+'.pdf')
    return

def make_emission_gif_plots():
    natural_base = '_nref11_RD0020_'
    refined_base = '_nref11n_nref10f_refine200kpc_z4to2_RD0020_'
    box_size = ds.arr(rb_width,'code_length').in_units('kpc')
    box_size = np.ceil(box_size/2.)
    res_list = [0.2,0.5,1,5,10]
    fontrc ={'fontname':'Helvetica','fontsize':20}
    mpl.rc('text', usetex=True)
    make_obs = True

    if make_obs:
        cmap = colors.ListedColormap(['Gray','HotPink','DarkTurquoise','Chartreuse'])
        bounds = [-5,1,2,3,5]
        norm = colors.BoundaryNorm(bounds,cmap.N)
    else:
        cmap = 'virdis'
        norm = None

    for line in lines:
        field = 'Emission_'+line
        for index in 'xyz':
            for res in res_list:
                if res == res_list[0]:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_forcedres.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_forcedres.cpkl'
                    pixsize = round(cosmo.arcsec_per_kpc_proper(redshift).value*0.182959,2)
                else:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_'+str(res)+'kpc.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'
                    pixsize = round(cosmo.arcsec_per_kpc_proper(redshift).value*res,2)

                frbNAT = cPickle.load(open(fileinNAT,'rb'))
                frbREF = cPickle.load(open(fileinREF,'rb'))

                if make_obs == True:
                    frbNAT = np.log10(frbNAT/(1.+redshift)**4)
                    frbREF = np.log10(frbREF/(1.+redshift)**4)
                else:
                    frbNAT = np.log10(frbNAT)
                    frbREF = np.log10(frbREF)

                bsL,bsR = -1*box_size.value,box_size.value

                fig,ax = plt.subplots(1,2)
                fig.set_size_inches(10,6)

                im = ax[0].imshow(frbNAT,extent=(bsL,bsR,bsR,bsL),
                                  vmin=-2,vmax=6,interpolation=None,
                                  cmap=cmap,norm=norm,origin='lower')
                im1= ax[1].imshow(frbREF,extent=(bsL,bsR,bsR,bsL),
                                  vmin=-2,vmax=6,interpolation=None,
                                  cmap=cmap,norm=norm,origin='lower')

                axins = inset_axes(ax[1],width="5%", height="100%",loc=3,
                                   bbox_to_anchor=(1.07, 0.0, 1, 1),
                                   bbox_transform=ax[1].transAxes,borderpad=0)

                ax[0].set_title('Natural',**fontrc)
                ax[1].set_title('Forced Refine',**fontrc)
                cb = fig.colorbar(im1, cax=axins,label=r'log( photons s$^{-1}$ cm$^{-2}$ sr$^{-1}$)')
                if line == 'CIII_977':
                    lineout = 'CIII 977'
                else:
                    lineout = line
                fig.suptitle('z=2, '+lineout+', '+str(res)+'kpc'+', '+str(pixsize)+'"',**fontrc)
                plt.savefig('z2_'+index+'_'+field+'_'+str(res)+'kpc_SBdim_obscol.pdf')
                plt.close()

    return

def make_phys_gif_plots():
    natural_base = '_nref11_RD0020_'
    refined_base = '_nref11n_nref10f_refine200kpc_z4to2_RD0020_'
    box_size = ds.arr(rb_width,'code_length').in_units('kpc')
    box_size = np.ceil(box_size/2.)
    res_list = [0.2,0.5,1,5,10]
    fontrc ={'fontname':'Helvetica','fontsize':20}
    mpl.rc('text', usetex=True)
    make_imshow = True
    properties = ['hden','temp']
    lims = {'hden':(-6,-1),'temp':(4.5,6.5)}

    for line in lines:
        field = 'Emission_'+line
        for index in 'xyz':
            for res in res_list:
                for prop in properties:
                    if res == res_list[0]:
                        fileinNAT = 'frbs/frb'+index+natural_base+prop+'_'+field+'_forcedres.cpkl'
                        fileinREF = 'frbs/frb'+index+refined_base+prop+'_'+field+'_forcedres.cpkl'
                        pixsize = round(cosmo.arcsec_per_kpc_proper(redshift).value*0.182959,2)
                    else:
                        fileinNAT = 'frbs/frb'+index+natural_base+prop+'_'+field+'_'+str(res)+'kpc.cpkl'
                        fileinREF = 'frbs/frb'+index+refined_base+prop+'_'+field+'_'+str(res)+'kpc.cpkl'
                        pixsize = round(cosmo.arcsec_per_kpc_proper(redshift).value*res,2)

                    frbNAT = cPickle.load(open(fileinNAT,'rb'))
                    frbREF = cPickle.load(open(fileinREF,'rb'))
                    frbNAT = np.log10(frbNAT/(1.+redshift)**4)
                    frbREF = np.log10(frbREF/(1.+redshift)**4)

                    bsL,bsR = -1*box_size.value,box_size.value

                    fig,ax = plt.subplots(1,2)
                    fig.set_size_inches(10,6)

                    im = ax[0].imshow(frbNAT,extent=(bsL,bsR,bsR,bsL),
                                    vmin=lims[prop][0],vmax=lims[prop][1],interpolation=None,
                                    cmap='viridis',origin='lower')
                    im1= ax[1].imshow(frbREF,extent=(bsL,bsR,bsR,bsL),
                                    vmin=lims[prop][0],vmax=lims[prop][0],interpolation=None,
                                    cmap='viridis',origin='lower')

                    axins = inset_axes(ax[1],
                            width="5%",  # width = 10% of parent_bbox width
                            height="100%",  # height : 50%
                            loc=3,bbox_to_anchor=(1.07, 0.0, 1, 1),bbox_transform=ax[1].transAxes,borderpad=0)

                    ax[0].set_title('Natural',**fontrc)
                    ax[1].set_title('Forced Refine',**fontrc)
                    cb = fig.colorbar(im1, cax=axins,label=r'log( photons s$^{-1}$ cm$^{-2}$ sr$^{-1}$)')
                    if line == 'CIII_977':
                        lineout = 'CIII 977'
                    else:
                        lineout = line
                        fig.suptitle('z=2, '+lineout+', '+str(res)+'kpc'+', '+str(pixsize)+'"',**fontrc)
                        plt.savefig('z2_'+index+'_'+prop+'_'+field+'_'+str(res)+'kpc.pdf')
                        plt.close()
    return

def make_radius_array(rb_width,frbarr):
	box_size = rb_width/2.
	num_cells = frbarr.shape[0]
	xL = np.linspace(-1*box_size,box_size,num_cells)
	xL2,yL = np.meshgrid(xL,xL)
	r = abs(xL2+1j*yL)
	#dr = np.abs([xL2[0,0]-xL2[0,1]])
	#radial = np.arange(box_size/dr)*dr + dr/2
	#nrad = len(radial)
	return r #,xL,dr,nrad

def plot_scatter_points_obscolors(ax,frbarr,r):
    colors = ['Chartreuse','DarkTurquoise','HotPink','Grey']
    lims = [10,3,2,1,-10]
    rnow = r.flatten()
    frbnow = frbarr.flatten()
    i = 0
    while i < len(lims)-1:
        idr = np.where((frbnow < lims[i]) & (frbnow > lims[i+1]))[0]
        ax.plot(rnow[idr],frbnow[idr],'.',colors=colors[i],alpha=0.37)
        i = i + 1
    return

def plot_SB_profiles(box_width):
    natural_base = '_nref11_RD0020_'
    refined_base = '_nref11n_nref10f_refine200kpc_z4to2_RD0020_'
    box_size = ds.arr(rb_width,'code_length').in_units('kpc')
    res_list = [0.2,1,5,10]
    fontrc={'fontname':'Helvetica','fontsize':20}
    mpl.rc('text',usetex=True)
    #rb_width = ds.arr(rb_width,'code_length').in_units('kpc')
    lines = ['CIII_977','OVI']

    for line in lines:
        field = 'Emission_'+line
        for index in 'xyz':
            fig,axes = plt.subplots(2,4)
            fig.set_size_inches(12,6)
            iax = 0
            for res in res_list:
                if res == res_list[0]:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_forcedres.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_forcedres.cpkl'
                else:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_'+str(res)+'kpc.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'

                frbNAT = cPickle.load(open(fileinNAT,'rb'))
                frbNAT = np.log10(frbNAT/(1+redshift)**4)
                frbREF = cPickle.load(open(fileinREF,'rb'))
                frbREF = np.log10(frbREF/(1+redshift)**4)

                ax = axes[0,iax]
                r = make_radius_array(box_width,frbREF)
                plot_scatter_points_obscolors(ax,frbREF,r)
                ax.set_ylim(-10,6)
                ax = axes[1,iax]
                r = make_radius_array(box_width,frbNAT)
                plot_scatter_points_obscolors(ax,frbNAT,r)
                ax.set_ylim(-10,6)
                iax = iax + 1

            plt.savefig('z2_'+index+'_'+field+'_SBradialplot.png')
    return


plot_SB_profiles(box_width)
