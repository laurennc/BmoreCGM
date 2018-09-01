import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import cPickle
import matplotlib as mpl
import trident
#import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
#from consistency import *
from yt.data_objects.particle_filters import add_particle_filter

from matplotlib.colorbar import Colorbar
from emission_functions import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.cosmology import WMAP9 as cosmo

#import holoviews as hv
#import pandas as pd
#import datashader as dshade
#from holoviews.operation.datashader import datashade, aggregate
#from holoviews import Store
#hv.extension('matplotlib')

def Stars(pfilter, data):
      return data[("all", "particle_type")] == 2
add_particle_filter("stars", function=Stars, filtered_type='all',
                    requires=["particle_type"])

#base = "/Users/dalek/data/Molly/natural/nref11"
base = "/Users/dalek/data/Molly/nref11n_nref10f_refine200kpc_z4to2"
fn = base+"/RD0016/RD0016"
lines = ['OVI','CIV','CIII_977','SiIV','HAlpha']
lines2 = ['O VI','C IV','C III _977','Si IV','H Alpha']
track_name = base+"/halo_track"
#track_name = base+"/halo_track_z2_z1"
args = fn.split('/')

detect_color_key = {b'nope':'#808080',
                    b'poss':'#FF69B4',
                    b'prob':'#00CED1',
                    b'defi': '#32CD32'} #'#7FFF00'}

line_color_key = {'OVI':'#2554C7','CIII_977':'#7D0552','CIV':'#438D80','SiIV':'#F75D59'}#,'HAlpha':''

detect_limits = {'nope':(-10,1),
                 'poss':(1,2),
                 'prob':(2,3),
                 'defi':(3,10)}

line_energies = {'CIII_977':4.*np.pi*2.03e-11,
                 'CIV':4.*np.pi*1.28e-11,
                 'OVI':4.*np.pi*1.92e-11,
                 'SiIV':4.*np.pi*1.42e-11,
                 'HAlpha':4.*np.pi*3.03e-12,
                 'LyAlpha':4.*np.pi*1.63e-11}

units_key = {'cell_mass':'Msun','Density':'g/cm**3','cell_volume':'kpc**3'}

roman_numerals =  {'1':'I','2':'II','3':'III','4':'IV','5':'V','6':'VI','7':'VII','8':'VIII'}
roman_numerals2 = {'I':'1','II':'2','III':'3','IV':'4','V':'5','VI':'6','VII':'7','VIII':'8'}

fontcs ={'fontname':'Helvetica','fontsize':16}
#mpl.rc('text', usetex=True)

### Fancy seaborn colormaps!! ###
#sns.set_style("white", {'axes.grid' : False})
#colors1 = plt.cm.Greys(np.linspace(0., 0.8, 192))
#sns_cmap2 = sns.blend_palette(('Pink','DeepPink','#1E90FF','DarkTurquoise','#50A638'), n_colors=40,as_cmap=True)
#colors2 = sns_cmap2(np.linspace(0.,1,64))
#colors = np.vstack((colors1, colors2))
#mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)


ds = yt.load(fn)
ds.add_particle_filter('stars')
track = Table.read(track_name, format='ascii')
track.sort('col1')
rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)
redshift = ds.current_redshift

box_width = ds.arr(rb_width,'code_length').in_units('kpc')

def add_grid(ax):
        ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
        ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
        ax.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75')
        ax.grid(which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.75',alpha=0.2)
        ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
        ax.grid(which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.75',alpha=0.2)
        return

def create_emission_frbs():
    dx = np.unique(rb['dx'])[0] #[1]
    #dx = ds.arr(0.1829591205,'kpc').in_units('code_length')
    dxs_list = [0.5,1,5,10]
    dxs_list = [ds.quan(q,'kpc').in_units('code_length') for q in dxs_list]
    res_list = np.append(dx,dxs_list)
    res_list = ds.arr(res_list,'code_length')

    for line in lines:
        field = 'Emission_'+line
        print line
        for index in 'xyz':
            print index
            for res in res_list:
                reskpc = round(res.in_units('kpc'),2)
                print reskpc,' kpc'
                num_cells = int(np.ceil(rb_width/res).value)
                #print 'num_cells: ',num_cells

                if res == res_list[0]:
                    fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_'+field+'_forcedres.cpkl'
                else:
                    fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_'+field+'_'+str(reskpc)+'kpc.cpkl'

                obj = ds.proj(('gas',field),index,data_source=rb)
                frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
                cPickle.dump(frb[('gas',field)],open(fileout,'wb'),protocol=-1)
    return

def create_coldens_frbs():
    dx = np.unique(rb['dx'])[1]
    #dx = ds.arr(0.1829591205,'kpc').in_units('code_length')
    dxs_list = [0.5,1,5,10]
    dxs_list = [ds.quan(q,'kpc').in_units('code_lenght') for q in dxs_list]
    res_list = np.append(dx,dxs_list)
    res_list = ds.arr(res_list,'code_length')

    trident.add_ion_fields(ds, ions=['Si II', 'Si III', 'Si IV',
                                     'C II', 'C III', 'C IV', 'O VI', 'Mg II'])
    lines = ['H_p0_number_density','Si_p1_number_density','Si_p2_number_density',
             'Si_p3_number_density','C_p1_number_density','C_p2_number_density',
             'C_p3_number_density','O_p5_number_density','Mg_p1_number_density']

    for line in lines:
        field = line
        print line
        for index in 'xyz':
            print index
            for res in res_list:
                reskpc = round(res.in_units('kpc'),2)
                print res,' kpc'
                num_cells = np.ceil(rb_width/res)

                if res == res_list[0]:
                    fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_'+field+'_forcedres.cpkl'
                else:
                    fileout = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_'+field+'_'+str(reskpc)+'kpc.cpkl'

                obj = ds.proj(('gas',field),index,data_source=rb)
                frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
                cPickle.dump(frb[('gas',field)],open(fileout,'wb'),protocol=-1)
    return

def create_phys_emis_weight_frbs():
    dx = np.unique(rb['dx'])[1]
    dxs_list = [0.5,1.0,5.0,10.0]
    dxs_list = [ds.quan(q,'kpc').in_units('code_length') for q in dxs_list]
    res_list = np.append(dx,dxs_list)
    res_list = ds.arr(res_list,'code_length')

    for line in lines2:
        field = 'Emission_'+line.replace(' ','')
        line_elements = line.split(' ')
        if (line_elements[0] != 'H'):
            tri_line = line_elements[0]+' '+line_elements[1]
            ion_num = str(int(roman_numerals2[line_elements[1]])-1)
            ion_field = line_elements[0]+'_p'+ion_num+'_ion_fraction'
            trident.add_ion_fields(ds,ions=[tri_line])
        else:
            trident.add_ion_fields(ds,ions=['H I'])
            ion_field = 'H_p0_ion_fraction'
        print line
        for index in 'xyz':
            print index
            for res in res_list:
                reskpc = round(res.in_units('kpc'),2)
                print reskpc,' kpc'
                num_cells = int(np.ceil(rb_width/res).value)

                if res == res_list[0]:
                    filehden = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_hden_'+field+'_forcedres.cpkl'
                    filetemp = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_temp_'+field+'_forcedres.cpkl'
                    fileion  = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_ionfrac_'+field+'_forcedres.cpkl'
                    filevel  = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_vel_'+field+'_forcedres.cpkl'
                else:
                    filehden = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_hden_'+field+'_'+str(reskpc)+'kpc.cpkl'
                    filetemp = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_temp_'+field+'_'+str(reskpc)+'kpc.cpkl'
                    fileion  = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_ionfrac_'+field+'_'+str(reskpc)+'kpc.cpkl'
                    filevel  = 'frb'+index+'_'+args[-3]+'_'+args[-2]+'_vel_'+field+'_'+str(reskpc)+'kpc.cpkl'

                obj = ds.proj('H_nuclei_density',index,data_source=rb,weight_field=field)
                frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
                cPickle.dump(frb['H_nuclei_density'],open(filehden,'wb'),protocol=-1)

                obj = ds.proj('temperature',index,data_source=rb,weight_field=field)
                frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
                cPickle.dump(frb['temperature'],open(filetemp,'wb'),protocol=-1)

                obj = ds.proj(('gas',ion_field),index,data_source=rb,weight_field=field)
                frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
                cPickle.dump(frb[('gas',ion_field)],open(fileion,'wb'),protocol=-1)

                bv = rb.quantities["BulkVelocity"]()
                rb.set_field_parameter("bulk_velocity",bv)
                obj = ds.proj(('gas','relative_velocity_'+index),index,data_source=rb,weight_field=field)
                frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
                cPickle.dump(frb[('gas','relative_velocity_'+index)].to('km/s').value,open(filevel,'wb'),protocol=-1)
    return

def plot_ytProjections():
    fileout = args[-3]+'_'+args[-2]+'_'+field+'_'+index+'_refwidth'
    proj = yt.ProjectionPlot(ds,index,('gas',field),width=(rb_width,'code_length'),
                             center=rb_center,data_source=rb)
    proj.set_cmap(('gas',field),pink_cmap)
    proj.save(fileout+'.pdf')
    return

def make_emission_gif_plots():
    natural_base = '_nref11_RD0016_'
    refined_base = '_nref11n_nref10f_refine200kpc_z4to2_RD0016_'
    nref11f_base = '_nref11f_refine200kpc_RD0016_'
    box_size = ds.arr(rb_width,'code_length').in_units('kpc')
    box_size = np.ceil(box_size/2.)
    res_list = [0.2,0.5,1.0,5.0,10.0]
    fontrc ={'fontname':'Helvetica','fontsize':20}
    mpl.rc('text', usetex=True)
    make_obs = True

    if make_obs:
        cmap = mymap
        norm = None
    else:
        cmap = 'viridis'
        norm = None

    for line in lines:
        field = 'Emission_'+line
        for index in 'xyz':
            for res in res_list:
                if res == res_list[0]:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_forcedres.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_forcedres.cpkl'
                    fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_forcedres.cpkl'
                    pixsize = round(cosmo.arcsec_per_kpc_proper(redshift).value*0.182959,2)
                else:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_'+str(res)+'kpc.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'
                    fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_'+str(res)+'kpc.cpkl'
                    pixsize = round(cosmo.arcsec_per_kpc_proper(redshift).value*res,2)

                frbNAT = cPickle.load(open(fileinNAT,'rb'))
                frbREF = cPickle.load(open(fileinREF,'rb'))
                frbN11 = cPickle.load(open(fileinN11,'rb'))

                if make_obs == True:
                    frbNAT = np.log10(frbNAT/(1.+redshift)**4)
                    frbREF = np.log10(frbREF/(1.+redshift)**4)
                    frbN11 = np.log10(frbN11/(1.+redshift)**4)
                else:
                    frbNAT = np.log10(frbNAT)
                    frbREF = np.log10(frbREF)
                    frbN11 = np.log10(frbN11)

                bsL,bsR = -1*box_size.value,box_size.value

                fig,ax = plt.subplots(1,3)
                fig.set_size_inches(14,6)


                ax[0].imshow(frbNAT,cmap=cmap,vmin=-5,vmax=3,
                             extent=(bsL,bsR,bsR,bsL),origin='lower',
                             interpolation=None)

                ax[1].imshow(frbREF,cmap=cmap,vmin=-5,vmax=3,
                             extent=(bsL,bsR,bsR,bsL),origin='lower',
                             interpolation=None)

                im2 = ax[2].imshow(frbN11,cmap=cmap,vmin=-5,vmax=3,
                             extent=(bsL,bsR,bsR,bsL),origin='lower',
                             interpolation=None)

                axins = inset_axes(ax[2],width="5%", height="100%",loc=3,
                                   bbox_to_anchor=(1.07, 0.0, 1, 1),
                                   bbox_transform=ax[2].transAxes,borderpad=0)

                ax[0].set_title('Natural',**fontrc)
                ax[1].set_title('Forced Refine 10',**fontrc)
                ax[2].set_title('Forced Refine 11',**fontrc)
                cb = fig.colorbar(im2, cax=axins,label=r'log( photons s$^{-1}$ cm$^{-2}$ sr$^{-1}$)')
                if line == 'CIII_977':
                    lineout = 'CIII 977'
                else:
                    lineout = line
                fig.suptitle('z=3, '+lineout+', '+str(res)+'kpc'+', '+str(pixsize)+'"',**fontrc)
                plt.savefig('z3_'+index+'_'+field+'_'+str(res)+'kpc_SBdim_obscol.pdf')
                plt.close()

    return

def make_small_emission_gif_plots():
    natural_base = '_nref11_RD0016_'
    refined_base = '_nref11n_nref10f_refine200kpc_z4to2_RD0016_'
    nref11f_base = '_nref11f_refine200kpc_RD0016_'
    box_size = ds.arr(rb_width,'code_length').in_units('kpc')
    box_size = np.ceil(box_size/2.)
    res_list = [0.2,0.5,1.0,5.0,10.0]
    fontrc ={'fontname':'Helvetica','fontsize':20}
    mpl.rc('text', usetex=True)
    make_obs = True

    if make_obs:
        cmap = mymap
        norm = None
    else:
        cmap = 'viridis'
        norm = None

    for line in lines:
        field = 'Emission_'+line
        for index in 'xyz':
            for res in res_list:
                if res == res_list[0]:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_forcedres.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_forcedres.cpkl'
                    fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_forcedres.cpkl'
                    pixsize = round(cosmo.arcsec_per_kpc_proper(redshift).value*0.182959,2)
                else:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_'+str(res)+'kpc.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'
                    fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_'+str(res)+'kpc.cpkl'
                    pixsize = round(cosmo.arcsec_per_kpc_proper(redshift).value*res,2)

                frbNAT = cPickle.load(open(fileinNAT,'rb'))
                frbREF = cPickle.load(open(fileinREF,'rb'))
                frbN11 = cPickle.load(open(fileinN11,'rb'))

                if make_obs == True:
                    frbNAT = np.log10(frbNAT/(1.+redshift)**4)
                    frbREF = np.log10(frbREF/(1.+redshift)**4)
                    frbN11 = np.log10(frbN11/(1.+redshift)**4)
                else:
                    frbNAT = np.log10(frbNAT)
                    frbREF = np.log10(frbREF)
                    frbN11 = np.log10(frbN11)

                bsL,bsR = -20.,20.
                num_pix = int(np.ceil(bsR/res))
                box_center = int(np.ceil(frbNAT.shape[0]/2.))
                iL,iR = box_center-num_pix,box_center+num_pix

                icA,icB = np.unravel_index(frbREF.argmax(),frbREF.shape)
                icA1,icA2 = icA-num_pix,icA+num_pix
                icB1,icB2 = icB-num_pix,icB+num_pix

                icL,icR = np.unravel_index(frbN11.argmax(),frbN11.shape)
                icL1,icL2 = icL-num_pix,icL+num_pix
                icR1,icR2 = icR-num_pix,icR+num_pix

                fig,ax = plt.subplots(1,3)
                fig.set_size_inches(14,6)


                ax[0].imshow(frbNAT[iL:iR+1,iL:iR+1],cmap=cmap,vmin=-5,vmax=3,
                             extent=(bsL,bsR,bsR,bsL),origin='lower',
                             interpolation=None)

                ax[1].imshow(frbREF[icA1:icA2+1,icB1:icB2+1],cmap=cmap,vmin=-5,vmax=3,
                             extent=(bsL,bsR,bsR,bsL),origin='lower',
                             interpolation=None)

                im2 = ax[2].imshow(frbN11[icL1:icL2+1,icR1:icR2+1],cmap=cmap,vmin=-5,vmax=3,
                             extent=(bsL,bsR,bsR,bsL),origin='lower',
                             interpolation=None)

                axins = inset_axes(ax[2],width="5%", height="100%",loc=3,
                                   bbox_to_anchor=(1.07, 0.0, 1, 1),
                                   bbox_transform=ax[2].transAxes,borderpad=0)

                ax[0].set_title('Natural',**fontrc)
                ax[1].set_title('Forced Refine 10',**fontrc)
                ax[2].set_title('Forced Refine 11',**fontrc)
                cb = fig.colorbar(im2, cax=axins,label=r'log( photons s$^{-1}$ cm$^{-2}$ sr$^{-1}$)')
                if line == 'CIII_977':
                    lineout = 'CIII 977'
                else:
                    lineout = line
                fig.suptitle('z=3, '+lineout+', '+str(res)+'kpc'+', '+str(pixsize)+'"',**fontrc)
                plt.savefig('z3_'+index+'_'+field+'_'+str(res)+'kpc_SBdim_obscol_ZOOM20kpc.pdf')
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
	dr = np.abs([xL2[0,0]-xL2[0,1]])
	radial = np.arange(box_size/dr)*dr + dr/2
	nrad = len(radial)
	return r,xL,dr,nrad,radial

def plot_scatter_points_obscolors(ax,frbarr,r,case=1):
    #colors = ['Chartreuse','DarkTurquoise','HotPink','Grey']
    colors = ['#32CD32','#00CED1','#FF69B4','#808080']
    lims = [10,3,2,1,-10]
    rnow = r.flatten()
    frbnow = frbarr.flatten()
    i = 0
    while i < len(lims)-1:
        idr = np.where((frbnow < lims[i]) & (frbnow > lims[i+1]))[0]
        if case == 1:
            if i == 0:
                ax.plot(rnow[idr],frbnow[idr],'.',color=colors[i],alpha=0.7,markersize=2)
            if i == 1:
                ax.plot(rnow[idr],frbnow[idr],'.',color=colors[i],alpha=0.25,markersize=1.75)
            else:
                ax.plot(rnow[idr],frbnow[idr],'.',color=colors[i],alpha=0.17,markersize=1.5)
        elif case == 2:
            ax.plot(rnow[idr],frbnow[idr],'.',color=colors[i],alpha=0.35,markersize=5)
        elif case == 3:
            ax.plot(rnow[idr],frbnow[idr],'.',color=colors[i])
        elif case == 4:
            ax.plot(rnow[idr],frbnow[idr],'.',color=colors[i],markersize=1.5,alpha=0.07)
        else:
            print 'CASE CAN ONLY BE GIVEN VALUES 1,2,3'
        i = i + 1
    return

def plot_SB_profiles(box_width):
    natural_base = '_nref11_RD0016_'
    refined_base = '_nref11n_nref10f_refine200kpc_z4to2_RD0016_'
    #nref11f_base = '_nref11f_refine200kpc_RD0016_'
    box_size = ds.arr(rb_width,'code_length').in_units('kpc')
    res_list = [0.2,1.0,5.0,10.0]
    fontrc={'fontname':'Helvetica','fontsize':20}
    mpl.rc('text',usetex=True)
    #rb_width = ds.arr(rb_width,'code_length').in_units('kpc')
    lines = ['CIII_977','OVI']

    case_dict = {'0.2':1,'1.0':2,'5.0':3,'10.0':3}

    for line in lines:
        field = 'Emission_'+line
        for index in 'xyz':
            fig,axes = plt.subplots(2,4,sharex=True,sharey=True)
            fig.set_size_inches(12,6)
            iax = 0
            for res in res_list:
                if res == res_list[0]:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_forcedres.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_forcedres.cpkl'
                    #fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_forcedres.cpkl'
                else:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_'+str(res)+'kpc.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'
                    #fileinN11 = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'

                frbNAT = cPickle.load(open(fileinNAT,'rb'))
                frbNAT = np.log10(frbNAT/(1+redshift)**4)
                frbREF = cPickle.load(open(fileinREF,'rb'))
                frbREF = np.log10(frbREF/(1+redshift)**4)
                #frbN11 = cPickle.load(open(fileinN11,'rb'))
                #frbN11 = np.log10(frbN11/(1+redshift)**4)

                ax = axes[1,iax]
                r,xL,dr,nrad,radial  = make_radius_array(box_width,frbREF)
                plot_scatter_points_obscolors(ax,frbREF,r,case=case_dict[str(res)])
                ax.set_ylim(-5,6)
                ax = axes[0,iax]
                r,xL,dr,nrad,radial  = make_radius_array(box_width,frbNAT)
                plot_scatter_points_obscolors(ax,frbNAT,r,case=case_dict[str(res)])
                ax.set_ylim(-5,6)
                #ax = axes[2,iax]
                #r,xL,dr,nrad,radial = make_radius_array(box_width,frbN11)
                #plot_scatter_points_obscolors(ax,frbN11,r,4) #case=case_dict[str(res)])
                #ax.set_ylim(-5,6)
                iax = iax + 1

            plt.tight_layout()
            plt.savefig('z3_'+index+'_'+field+'_SBradialplot.png')
    return


def plot_SB_profiles_all_lines(box_width):
    natural_base = '_nref11_RD0016_'
    refined_base = '_nref11n_nref10f_refine200kpc_z4to2_RD0016_'
    nref11f_base = '_nref11f_refine200kpc_RD0016_'
    box_size = ds.arr(rb_width,'code_length').in_units('kpc')
    res_list = [0.2,1.0,5.0,10.0]
    fontrc={'fontname':'Helvetica','fontsize':20}
    mpl.rc('text',usetex=True)
    #rb_width = ds.arr(rb_width,'code_length').in_units('kpc')
    lines = ['SiIV','CIII_977','CIV','OVI']

    case_dict = {'0.2':1,'1.0':2,'5.0':3,'10.0':3}

    for index in 'xyz':
        #fig,axes = plt.subplots(2,4,sharex=True,sharey=True)
        #fig.set_size_inches(12,6)
        fig,axes = plt.subplots(3,4,sharex=True,sharey=True)
        fig.set_size_inches(12,9)
        iax = 0
        for line in lines:
            field = 'Emission_'+line

            fileinNAT = 'frbs/frb'+index+natural_base+field+'_forcedres.cpkl'
            fileinREF = 'frbs/frb'+index+refined_base+field+'_forcedres.cpkl'
            fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_forcedres.cpkl'

            #fileinNAT = 'frbs/frb'+index+natural_base+field+'_'+str(res)+'kpc.cpkl'
            #fileinREF = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'
            #fileinN11 = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'

            frbNAT = cPickle.load(open(fileinNAT,'rb'))
            frbNAT = np.log10(frbNAT/(1+redshift)**4)
            frbREF = cPickle.load(open(fileinREF,'rb'))
            frbREF = np.log10(frbREF/(1+redshift)**4)
            frbN11 = cPickle.load(open(fileinN11,'rb'))
            frbN11 = np.log10(frbN11/(1+redshift)**4)

            ax = axes[1,iax]
            r,xL,dr,nrad,radial  = make_radius_array(box_width,frbREF)
            plot_scatter_points_obscolors(ax,frbREF,r,case=1)
            ax.set_ylim(-5,6)
            ax = axes[0,iax]
            r,xL,dr,nrad,radial  = make_radius_array(box_width,frbNAT)
            plot_scatter_points_obscolors(ax,frbNAT,r,case=1)
            ax.set_ylim(-5,6)
            ax = axes[2,iax]
            r,xL,dr,nrad,radial = make_radius_array(box_width,frbN11)
            plot_scatter_points_obscolors(ax,frbN11,r,4) #case=case_dict[str(res)])
            ax.set_ylim(-5,6)
            iax = iax + 1

        plt.tight_layout()
        plt.savefig('z3_'+index+'_LINES_SBradialplot_ALLSIMS.png')
    return


def make_emis_pdframe(frb,r):
    emis = frb.flatten()
    dist = r.flatten()
    detect_prob = np.chararray(np.size(dist), 4)
    detect_prob[emis >= 3] = 'definite'
    detect_prob[((emis < 3) & (emis >= 2))] = 'probable'
    detect_prob[((emis < 2) & (emis >= 1))] = 'possible'
    detect_prob[emis < 1] = 'nope'
    df = pd.DataFrame({'dist':dist,'emis':emis,'detect_prob':detect_prob})
    df.detec_prob = df.detect_prob.astype('category')
    return df

def hvPoints_by_detect_prob(df,prob,hvplot=None):
    idx = ((df['emis'] >= detect_limits[prob][0]) & (df['emis'] < detect_limits[prob][1]))
    ptshere = np.vstack((df['dist'][idx], df['emis'][idx])).T
    pltout = hv.Points(ptshere,kdims=['dist','emis']).opts(style=dict(color=detect_color_key[prob]))
    if hvplot == None:
        retplot = pltout
    else:
        retplot = hvplot * pltout
    return retplot

def holoviews_SB_profiles(box_width):
    renderer = Store.renderers['matplotlib'].instance(fig='pdf', holomap='gif')
    natural_base = '_nref11_RD0016_'
    refined_base = '_nref11n_nref10f_refine200kpc_z4to2_RD0016_'
    nref11f_base = '_nref11f_refine200kpc_RD0016_'
    box_size = ds.arr(rb_width,'code_length').in_units('kpc')
    res_list = [0.2,1.0,5.0,10.0]
    fontrc={'fontname':'Helvetica','fontsize':20}
    mpl.rc('text',usetex=True)
    #rb_width = ds.arr(rb_width,'code_length').in_units('kpc')
    lines = ['CIII_977','OVI','CIV','HAlpha','SiIV']

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
                    fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_forcedres.cpkl'
                else:
                    fileinNAT = 'frbs/frb'+index+natural_base+field+'_'+str(res)+'kpc.cpkl'
                    fileinREF = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'
                    fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_'+str(res)+'kpc.cpkl'

                frbNAT = cPickle.load(open(fileinNAT,'rb'))
                frbNAT = np.log10(frbNAT/(1+redshift)**4)
                frbREF = cPickle.load(open(fileinREF,'rb'))
                frbREF = np.log10(frbREF/(1+redshift)**4)
                frbN11 = cPickle.load(open(fileinN11,'rb'))
                frbN11 = np.log10(frbN11/(1+redshift)**4)

                r,xL,dr,nrad,radial  = make_radius_array(box_width,frbNAT)
                r2,xL,dr,nrad,radial = make_radius_array(box_width,frbN11)
                dfREF = make_emis_pdframe(frbREF,r)
                dfNAT = make_emis_pdframe(frbNAT,r)
                dfN11 = make_emis_pdframe(frbN11,r2)

                if (res < 2):
                    llim = dfREF['emis'].min()-0.1
                    pREF = hv.Points(dfREF,kdims=['dist','emis'])
                    shadeREF = datashade(pREF,color_key=detect_color_key,
                                    aggregator=dshade.count_cat('detect_prob'), y_range=(llim,6),
                                    dynamic=False).opts(plot=dict(aspect='square'))

                    pNAT =  hv.Points(dfNAT,kdims=['dist','emis'])
                    shadeNAT = datashade(pNAT,color_key=detect_color_key,
                                    aggregator=dshade.count_cat('detect_prob'), y_range=(llim,6),
                                    dynamic=False).opts(plot=dict(aspect='square'))

                    pN11 = hv.Points(dfN11,kdims=['dist','emis'])
                    shadeN11 = datashade(pN11,color_key=detect_color_key,
                                    aggregator=dshade.count_cat('detect_prob'),y_range=(llim,6),
                                    dynamic=False).opts(plot=dict(aspect='square'))

                    if res == res_list[0]:
                        colsREF,colsNAT,colsN11 = shadeREF,shadeNAT,shadeN11
                    else:
                        colsREF = colsREF + shadeREF
                        colsNAT = colsNAT + shadeNAT
                        colsN11 = colsN11 + shadeN11

                else:
                    i = 0
                    while i < len(detect_limits.keys()):
                        if i == 0:
                            pltREF = hvPoints_by_detect_prob(dfREF,detect_limits.keys()[i])
                            pltNAT = hvPoints_by_detect_prob(dfNAT,detect_limits.keys()[i])
                            pltN11 = hvPoints_by_detect_prob(dfN11,detect_limits.keys()[i])
                        else:
                            pltREF = hvPoints_by_detect_prob(dfREF,detect_limits.keys()[i],hvplot=pltREF)
                            pltNAT = hvPoints_by_detect_prob(dfNAT,detect_limits.keys()[i],hvplot=pltNAT)
                            pltN11 = hvPoints_by_detect_prob(dfN11,detect_limits.keys()[i],hvplot=pltN11)
                        i = i + 1

                    colsREF = colsREF + pltREF
                    colsNAT = colsNAT + pltNAT
                    colsN11 = colsN11 + pltN11

            pltout = (colsREF + colsNAT + colsN11).cols(len(res_list))
            fileout= 'SBprofile_z3_'+index+'_'+field
            renderer.save(pltout, fileout)
    return

def holoviews_radial_profiles(weight_by=None):
    dens = np.log10(rb['H_nuclei_density'])
    temp = np.log10(rb['Temperature'])
    Zgas = np.log10(rb['metallicity'])
    cell_mass = rb['cell_mass'].in_units('Msun')
    cell_volume = rb['cell_volume'].in_units('kpc**3')
    x = rb['x']
    y = rb['y']
    z = rb['z']

    halo_center = ds.arr(rb_center,'code_length')
    dist = np.sqrt((halo_center[0]-rb['x'])**2.+(halo_center[1]-rb['y'])**2.+(halo_center[2]-rb['z'])**2.).in_units('kpc')

    df = pd.DataFrame({'temp':temp, 'dens':dens, 'Zgas':Zgas,'cell_volume':cell_volume,
                        'x':x,'y':y,'z':z,'dist':dist,'cell_mass':cell_mass})

    temp_dist = hv.Scatter(df,kdims=['dist'],vdims=['temp'],label="Temperature ")
    dens_dist = hv.Scatter(df,kdims=['dist'],vdims=['dens'],label='Hydrogen Number Density')
    metal_dist = hv.Scatter(df,kdims=['dist'],vdims=['Zgas'],label='Metallicity')

    if weight_by == None:
        dist_plots = (datashade(temp_dist,cmap=cm.Reds, dynamic=False,x_range=(0,60),y_range=(2,8.4)).opts(plot=dict(aspect='square'))
                    + datashade(dens_dist,cmap=cm.Blues, dynamic=False,x_range=(0,60),y_range=(-6.5,2)).opts(plot=dict(aspect='square'))
                    + datashade(metal_dist,cmap=cm.BuGn, dynamic=False,x_range=(0,60),y_range=(-8.5,1.4)).opts(plot=dict(aspect='square')))
        fileout= 'basic_profile_'+args[-3]+'_'+args[-1]

    if weight_by == 'cell_mass':
        temp_shade = aggregate(hv.Scatter(df,['dist','temp']),y_range=(2,8.4),aggregator=dshade.sum('cell_mass'))
        temp_shade = temp_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.Reds))
        dens_shade = aggregate(hv.Scatter(df,['dist','dens']),y_range=(-7,2.5),aggregator=dshade.sum('cell_mass'))
        dens_shade = dens_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.Blues))
        metal_shade = aggregate(hv.Scatter(df,['dist','Zgas']),y_range=(-7,2.5),aggregator=dshade.sum('cell_mass'))
        metal_shade = metal_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.BuGn))

        dist_plots = (temp_shade + dens_shade + metal_shade)
        fileout = 'basic_profile_cell_mass_'+args[-3]+'_'+args[-1]

    if weight_by == 'cell_volume':
        temp_shade = aggregate(hv.Scatter(df,['dist','temp']),y_range=(2,8.4),aggregator=dshade.sum('cell_volume'))
        temp_shade = temp_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.Reds))
        dens_shade = aggregate(hv.Scatter(df,['dist','dens']),y_range=(-7,2.5),aggregator=dshade.sum('cell_volume'))
        dens_shade = dens_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.Blues))
        metal_shade = aggregate(hv.Scatter(df,['dist','Zgas']),y_range=(-7,2.5),aggregator=dshade.sum('cell_volume'))
        metal_shade = metal_shade.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cm.BuGn))

        dist_plots = (temp_shade + dens_shade + metal_shade)
        fileout = 'basic_profile_cell_vol_'+args[-3]+'_'+args[-1]

    renderer = Store.renderers['matplotlib'].instance(fig='pdf', holomap='gif')
    renderer.save(dist_plots, fileout)
    return

def holoviews_general_plot(xfield,yfield,ranges,fout,cmap,weight_by=None):
    xx = np.log10(rb[xfield])
    yy = np.log10(rb[yfield])

    if weight_by != None:
        weight = rb[weight_by]

    df = pd.DataFrame({xfield:xx, yfield:yy, weight_by:weight})

    if weight_by == None:
        scatterplot = hv.Scatter(df,kdims=[xfield],vdims=[yfield])
        full_plot = (datashade(scatterplot,cmap=cm.Plasma,dynamic=False,x_range=(ranges[0],ranges[1]),
                     y_range=(ranges[2],ranges[3]))).opts(plot=dict(aspect='square'))
    else:
        full_plot = aggregate(hv.Scatter(df,[xfield,yfield]),y_range=(ranges[2],ranges[3]),aggregator=dshade.sum(weight_by))
        full_plot = full_plot.opts(plot=dict(colorbar=True,aspect='square',logz=True),style=dict(cmap=cmap))

    renderer = Store.renderers['matplotlib'].instance(fig='pdf', holomap='gif')
    renderer.save(full_plot, fout)
    return

def covering_fraction_by_radius(rb_width,SB_cutoff=1.):
    natural_base = '_nref11_RD0016_'
    refined_base = '_nref11n_nref10f_refine200kpc_z4to2_RD0016_'
    nref11f_base = '_nref11f_refine200kpc_z4to2_RD0016_'

    rb_width = ds.quan(rb_width,'code_length').in_units('kpc')
    rb_width = float(rb_width.value)

    fig,ax = plt.subplots(1,3,sharey=True)
    fig.set_size_inches(12,4)

    lines = ['OVI','CIV','CIII_977','SiIV']
    res_list = ['Native',1.0,10.0]
    i,index,res = 1,'x',0.5

    #while i < len(res_list):
        #res = res_list[i]
    for line in lines:
            field = 'Emission_'+line
            if i == 0:
                fileinNAT = 'frbs/frb'+index+natural_base+field+'_forcedres.cpkl'
                fileinREF = 'frbs/frb'+index+refined_base+field+'_forcedres.cpkl'
                fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_forcedres.cpkl'
            else:
                fileinNAT = 'frbs/frb'+index+natural_base+field+'_'+str(res)+'kpc.cpkl'
                fileinREF = 'frbs/frb'+index+refined_base+field+'_'+str(res)+'kpc.cpkl'
                fileinN11 = 'frbs/frb'+index+nref11f_base+field+'_'+str(res)+'kpc.cpkl'

            frbNAT = cPickle.load(open(fileinNAT,'rb'))
            frbREF = cPickle.load(open(fileinREF,'rb'))
            frbN11 = cPickle.load(open(fileinN11,'rb'))

            frbNAT = np.log10(frbNAT/(1.+redshift)**4)
            frbREF = np.log10(frbREF/(1.+redshift)**4)
            frbN11 = np.log10(frbN11/(1.+redshift)**4)

            r,xL,dr,nrad,radial = make_radius_array(rb_width,frbREF)
            r2,xL2,dr2,nrad2,radial2 = make_radius_array(rb_width,frbN11)
            cfNAT,cfREF,cfN11 = np.zeros(nrad),np.zeros(nrad),np.zeros(nrad2)

            irad = 0
            while irad < nrad:
                minrad = irad*dr
                maxrad = minrad + dr

                idx = np.where((r>=minrad) & (r<maxrad) & (frbNAT > SB_cutoff))[0]
                idp = np.where((r>=minrad) & (r<maxrad))[0]
                cfNAT[irad] = float(len(idx))/float(len(idp))

                idx = np.where((r>=minrad) & (r<maxrad) & (frbREF > SB_cutoff))[0]
                cfREF[irad] = float(len(idx))/float(len(idp))
                irad = irad + 1

            irad = 0
            while irad < nrad2:
                minrad = irad*dr2
                maxrad = minrad + dr2
                idx = np.where((r2>=minrad) & (r2<maxrad) & (frbN11 > SB_cutoff))[0]
                idp = np.where((r2>=minrad) & (r2<maxrad))[0]
                cfN11[irad] = float(len(idx))/float(len(idp))
                irad = irad + 1

            ax[0].plot(radial,cfREF,lw=1.5,color=line_color_key[line],label=line)
            ax[1].plot(radial,cfNAT,lw=1.5,ls='--',color=line_color_key[line])
            ax[2].plot(radial2,cfN11,lw=1.5,ls='-.',color=line_color_key[line])
            #ax[i].set_title(str(res_list[i]))

            #i = i + 1
    #print radial
    ax[1].set_xlabel('Radius [kpc]')
    ax[0].set_ylabel('Covering Fraction')
    ax[0].legend()
    plt.savefig('testing_covering_frac_0.5kpc.pdf')

    return

def cumulative_distribution_function(ax,line,res,rb_width):
    field = 'Emission_'+line
    bases = ['_nref11_RD0016_','_nref11n_nref10f_refine200kpc_z4to2_RD0016_','_nref11f_refine200kpc_RD0016_']
    #colors = ['3 colors here!']

    rb_width = ds.quan(rb_width,'code_length').in_units('kpc')
    rb_width = float(rb_width.value)

    bins = np.linspace(-1.5,3.5,11)

    ls = iter(['--', '-', '-'])
    lw = iter([1.5,3.0,1.5])
    indices = 'xyz'

    hist = np.zeros(len(bins))
    tot_pix = 0
    pix_hist = np.zeros(len(bins))


    for base in bases:
        for index in indices:
            if res == 'forced':
                filein = 'frbs/frb'+index+base+field+'_forcedres.cpkl'
            else:
                filein = 'frbs/frb'+index+base+field+'_'+str(res)+'kpc.cpkl'

            frb = cPickle.load(open(filein,'rb'))
            frb = np.log10(frb/(1.+redshift)**4)
            frb = frb.flatten()
            tot_pix += len(frb)

            i = 0
            while i < len(bins):
                if i == 0:
                    num_pix = np.where(frb < bins[i])[0]
                    pix_hist[i] += len(num_pix)
                elif i == len(bins)-1:
                    num_pix = np.where(frb >= bins[i])[0]
                    pix_hist[i] += len(num_pix)
                else:
                    num_pix = np.where((frb < bins[i+1]) & (frb >= bins[i]))[0]
                    #hist[i] = float(len(num_pix))/len(frb)
                    pix_hist[i] += len(num_pix)
                i = i + 1

        for i in range(len(hist)):
            hist[i] = float(pix_hist[i]/tot_pix)

        label_out = base.split('_')[1]
        ls_here,lw_here = ls.next(),lw.next()
        ax.plot(bins,np.log10(hist),label=label_out,color=line_color_key[line],linestyle=ls_here,lw=lw_here)
        #ax.plot(bins,hist,color=line_color_key[line],marker='s',linestyle=ls_here,lw=1.5)

    #ax.set_xlabel('SB')
    #ax.set_ylabel('Pix Fraction')
    ax.set_xlim(-1,3)
    ax.set_ylim(-5.5,-0.5)
    ax.legend()
    #plt.savefig('SB_cdf_RD0016_'+field+'_forcedres.pdf')
    #plt.close()

    return

def cdf_plot_loop():
    fig,ax = plt.subplots(2,2,sharex=True,sharey=True)
    fig.set_size_inches(6,6)
    ax = ax.flatten()
    subplts = iter(ax)
    lines = ['OVI','CIV','CIII_977','SiIV']
    for line in lines:
        cumulative_distribution_function(subplts.next(),line,'forced',rb_width)

    plt.tight_layout()
    #plt.savefig('SB_'+index+'_cdf_RD0016_all_forced_LOG.pdf')
    plt.savefig('SB_all_coverfrac_RD0016_forced_LOG.pdf')
    plt.close()
    return

def hden_temp_hist(index,base,RD,line,resolution,redshift):
    hden_file = '/Users/dalek/repos/BmoreCGM/frbs/frb'+index+base+RD+'_hden_Emission_'+line+'_'+resolution+'.cpkl'
    temp_file = '/Users/dalek/repos/BmoreCGM/frbs/frb'+index+base+RD+'_temp_Emission_'+line+'_'+resolution+'.cpkl'
    emis_file = '/Users/dalek/repos/BmoreCGM/frbs/frb'+index+base+RD+'_Emission_'+line+'_'+resolution+'.cpkl'

    hden = cPickle.load(open(hden_file,'rb'))
    temp = cPickle.load(open(temp_file,'rb'))
    emis = cPickle.load(open(emis_file,'rb'))
    emis = emis/(1.+redshift)**4.0

    hden,temp,emis = np.log10(hden.flat),np.log10(temp.flat),np.log10(emis.flat)

    hist, xedges, yedges  = np.histogram2d(hden,temp,range=[[-6,1],[3.5, 6.5]],bins=120)#,emis=coldens)

    xbins = np.digitize(hden,xedges[1:-1])
    ybins = np.digitize(temp,yedges[1:-1])

    totvals = np.zeros((len(xedges[1:]),len(yedges[1:])))
    numvals = np.zeros((len(xedges[1:]),len(yedges[1:])))
    maxvals = np.zeros((len(xedges[1:]),len(yedges[1:])))

    for i in range(len(emis)):
        totvals[xbins[i],ybins[i]] += 10.**emis[i]
        numvals[xbins[i],ybins[i]] += 1
        if 10.**emis[i] > maxvals[xbins[i],ybins[i]]:
            maxvals[xbins[i],ybins[i]] = 10.**emis[i]

    avgvals = totvals/numvals
    idx = np.where(totvals == 0)
    avgvals[idx] = 0

    ## I'm going to want to plot the average values on a plot but the issue is
    ## that the nans throw off the averages. But there are so many arrays to
    ## parse between it's easiest just to do it at the end
    averages = np.zeros(4)
    idA = np.where(emis > 1.)[0]
    averages[2],averages[3] = np.average(hden[idA]),np.average(temp[idA])

    id1 = np.argwhere(np.isnan(hden))
    id2 = np.argwhere(np.isnan(temp))

    hden = np.delete(hden, id1)
    temp = np.delete(temp, id2)
    averages[0],averages[1] = np.average(hden),np.average(temp)



    return hist,avgvals,maxvals,averages

def add_phase_histograms(ax,field,index,base,RD,resolution,redshift,line):
    hden_file = '/Users/dalek/repos/BmoreCGM/frbs/frb'+index+base+RD+'_hden_Emission_'+line+'_'+resolution+'.cpkl'
    temp_file = '/Users/dalek/repos/BmoreCGM/frbs/frb'+index+base+RD+'_temp_Emission_'+line+'_'+resolution+'.cpkl'
    emis_file = '/Users/dalek/repos/BmoreCGM/frbs/frb'+index+base+RD+'_Emission_'+line+'_'+resolution+'.cpkl'

    hden = cPickle.load(open(hden_file,'rb'))
    temp = cPickle.load(open(temp_file,'rb'))
    emis = cPickle.load(open(emis_file,'rb'))
    emis = emis/(1.+redshift)**4.0

    hden,temp,emis = np.log10(hden.flat),np.log10(temp.flat),np.log10(emis.flat)

    id1 = np.argwhere(np.isnan(hden))
    id2 = np.argwhere(np.isnan(temp))
    hden = np.delete(hden, id1)
    temp = np.delete(temp, id2)
    emis = np.delete(emis, id1)

    idN = np.where(emis < 1.)[0]
    idP = np.where((emis >= 1.) & (emis < 2.))[0]
    idB = np.where((emis >= 2.) & (emis < 3.))[0]
    idG = np.where((emis >= 3.))[0]

    if field == 'temp':
        bins = np.arange(3.5,6.5,0.1)
        ax.hist(temp[idN],bins=bins, orientation='horizontal', color='k', edgecolor='w',alpha=0.5,normed=True)
        ax.hist(temp[idP],bins=bins, orientation='horizontal', color='DeepPink', edgecolor='w',alpha=0.5,normed=True)
        ax.hist(temp[idB],bins=bins, orientation='horizontal', color='DarkTurquoise', edgecolor='w',alpha=0.5,normed=True)
        ax.hist(temp[idG],bins=bins, orientation='horizontal', color='#50A638', edgecolor='w',alpha=0.5,normed=True)

        ax.set_yticks(np.linspace(3.5,6.5,7))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_ylim([3.5,6.5])
        ax.set_axis_off()

    elif field == 'hden':
        bins = np.linspace(-6,1,50)
        ax.hist(hden[idN],bins=bins, orientation='vertical', color='k', edgecolor='w',alpha=0.5,normed=True)
        ax.hist(hden[idP],bins=bins, orientation='vertical', color='DeepPink', edgecolor='w',alpha=0.5,normed=True)
        ax.hist(hden[idB],bins=bins, orientation='vertical', color='DarkTurquoise', edgecolor='w',alpha=0.5,normed=True)
        ax.hist(hden[idG],bins=bins, orientation='vertical', color='#50A638', edgecolor='w',alpha=0.5,normed=True)

        ax.set_xticks([]) # Ensures we have the same ticks as the scatter plot !
        ax.set_yticklabels([])
        ax.set_xlim([-6,1])
        ax.set_axis_off()

    else:
        print 'NEED TO GIVE AN APPROPRIATE FIELD'

    return

def make_weighted_phase_diagrams(index,base,RD,resolution,redshift):
    ds = yt.load(fn)
    #lines = ['CIII_977']
    for line in lines:
        hist,avgvals,maxvals,averages = hden_temp_hist(index,base,RD,line,resolution,ds.current_redshift)

        fig = plt.figure(1, figsize=(9,8))
        gs = gridspec.GridSpec(3,3, height_ratios=[0.2,1,0.05], width_ratios=[1,0.2,0.2])
        gs.update(left=0.075, right=0.93, bottom=0.09, top=0.95, wspace=0.00, hspace=0.00)

        # --------------------------------------------------------
        ax1 = plt.subplot(gs[1,0])
        # --------------------------------------------------------
        plt1 = ax1.imshow(np.log10(avgvals.T),extent=[-6,1,3.5,6.5],
                          interpolation='nearest',origin='lower',cmap=mymap,vmin=-5,vmax=3)

        add_grid(ax1)
        ax1.set_xlim([-6,1])
        ax1.set_ylim([3.5,6.5])
        ax1.set_xlabel(r' ') # Force this empty !
        ax1.set_xticks(np.linspace(-6,1,8)) # Force this for consistency with hists!
        ax1.set_xticklabels(np.linspace(-6,1,8))
        ax1.set_ylabel(r'Temperature [log(K)]')
        ax1.set_xlabel(r'Hydrogen Number Density [log(cm^{-3})]')
        x0,x1 = ax1.get_xlim()
        y0,y1 = ax1.get_ylim()
        ax1.set_aspect((x1-x0)/(y1-y0))

        # --------------------------------------------------------
        #cbax = plt.subplot(gs[2,0]) # Place it where it should be.
        # --------------------------------------------------------
        #cb = Colorbar(ax = cbax, mappable = plt1, orientation = 'horizontal', ticklocation = 'bottom')
        #cb.set_label(r'Colorbar !', labelpad=5)

        # --------------------------------------------------------
        ax1v = plt.subplot(gs[1,1])
        # --------------------------------------------------------
        add_phase_histograms(ax1v,'temp',index,base,RD,resolution,redshift,line)

        # --------------------------------------------------------
        ax1h = plt.subplot(gs[0,0])
        # --------------------------------------------------------
        add_phase_histograms(ax1h,'hden',index,base,RD,resolution,redshift,line)

        plt.savefig('phase'+base+RD+'_weighted_Emission_'+line+'_'+resolution+'.pdf')
        plt.close()
    return

def hden_temp_ionfrac_hist(index,base,RD,line,resolution,redshift):
    hden_file = '/Users/dalek/repos/BmoreCGM/frbs/frb'+index+base+RD+'_hden_Emission_'+line+'_'+resolution+'.cpkl'
    temp_file = '/Users/dalek/repos/BmoreCGM/frbs/frb'+index+base+RD+'_temp_Emission_'+line+'_'+resolution+'.cpkl'
    ionf_file = '/Users/dalek/repos/BmoreCGM/frbs/frb'+index+base+RD+'_ionfrac_Emission_'+line+'_'+resolution+'.cpkl'

    hden = cPickle.load(open(hden_file,'rb'))
    temp = cPickle.load(open(temp_file,'rb'))
    ionf = cPickle.load(open(ionf_file,'rb'))

    hden,temp,ionf = np.log10(hden.flat),np.log10(temp.flat),ionf.flat

    hist, xedges, yedges  = np.histogram2d(hden,temp,range=[[-6,1],[3.5, 6.5]],bins=120)#,emis=coldens)

    xbins = np.digitize(hden,xedges[1:-1])
    ybins = np.digitize(temp,yedges[1:-1])

    totvals = np.zeros((len(xedges[1:]),len(yedges[1:])))
    numvals = np.zeros((len(xedges[1:]),len(yedges[1:])))
    maxvals = np.zeros((len(xedges[1:]),len(yedges[1:])))

    for i in range(len(ionf)):
        totvals[xbins[i],ybins[i]] += ionf[i]
        numvals[xbins[i],ybins[i]] += 1
        if ionf[i] > maxvals[xbins[i],ybins[i]]:
            maxvals[xbins[i],ybins[i]] = ionf[i]

    avgvals = totvals/numvals
    idx = np.where(totvals == 0)
    avgvals[idx] = 0

    ## I'm going to want to plot the average values on a plot but the issue is
    ## that the nans throw off the averages. But there are so many arrays to
    ## parse between it's easiest just to do it at the end
    #averages = np.zeros(4)
    #idA = np.where(emis > 1.)[0]
    #averages[2],averages[3] = np.average(hden[idA]),np.average(temp[idA])

    #id1 = np.argwhere(np.isnan(hden))
    #id2 = np.argwhere(np.isnan(temp))

    #hden = np.delete(hden, id1)
    #temp = np.delete(temp, id2)
    #averages[0],averages[1] = np.average(hden),np.average(temp)

    return hist,avgvals,maxvals #,averages

ion_lims = {'OVI':(0.,0.2),'CIV':(0.,0.25)}
def make_ionfrac_weighted_phase_diagrams(index,base,RD,resolution,redshift):
    ds = yt.load(fn)
    #lines = ['CIII_977']
    for line in lines:
        hist,avgvals,maxvals = hden_temp_ionfrac_hist(index,base,RD,line,resolution,ds.current_redshift)

        fig = plt.figure(1, figsize=(9,8))
        gs = gridspec.GridSpec(3,3, height_ratios=[0.2,1,0.05], width_ratios=[1,0.2,0.2])
        gs.update(left=0.075, right=0.93, bottom=0.09, top=0.95, wspace=0.00, hspace=0.00)

        # --------------------------------------------------------
        ax1 = plt.subplot(gs[1,0])
        # --------------------------------------------------------
        if line in ion_lims.keys():
            plt1 = ax1.imshow(avgvals.T,extent=[-6,1,3.5,6.5],vmin=ion_lims[line][0],vmax=ion_lims[line][1],
                          interpolation='nearest',origin='lower',cmap='Purples')#mymap,vmin=-5,vmax=3)
        else:
            plt1 = ax1.imshow(avgvals.T,extent=[-6,1,3.5,6.5],
                          interpolation='nearest',origin='lower',cmap='Purples')
        add_grid(ax1)
        ax1.set_xlim([-6,1])
        ax1.set_ylim([3.5,6.5])
        ax1.set_xlabel(r' ') # Force this empty !
        ax1.set_xticks(np.linspace(-6,1,8)) # Force this for consistency with hists!
        ax1.set_xticklabels(np.linspace(-6,1,8))
        ax1.set_ylabel(r'Temperature [log(K)]')
        ax1.set_xlabel(r'Hydrogen Number Density [log(cm^{-3})]')
        x0,x1 = ax1.get_xlim()
        y0,y1 = ax1.get_ylim()
        ax1.set_aspect((x1-x0)/(y1-y0))
        plt.colorbar(plt1)
        # --------------------------------------------------------
        #cbax = plt.subplot(gs[2,0]) # Place it where it should be.
        # --------------------------------------------------------
        #cb = Colorbar(ax = cbax, mappable = plt1, orientation = 'horizontal', ticklocation = 'bottom')
        #cb.set_label(r'Colorbar !', labelpad=5)

        # --------------------------------------------------------
        #ax1v = plt.subplot(gs[1,1])
        # --------------------------------------------------------
        #add_phase_histograms(ax1v,'temp',index,base,RD,resolution,redshift,line)

        # --------------------------------------------------------
        #ax1h = plt.subplot(gs[0,0])
        # --------------------------------------------------------
        #add_phase_histograms(ax1h,'hden',index,base,RD,resolution,redshift,line)

        plt.savefig('phase'+base+RD+'_weighted_Emission_'+line+'_ionfrac_'+resolution+'.pdf')
        plt.close()
    return



#make_weighted_phase_diagrams('x','_nref11n_nref10f_refine200kpc_z4to2_','RD0020','forcedres',ds.current_redshift)
#make_weighted_phase_diagrams('x','_nref11_','RD0016','forcedres',ds.current_redshift)
#make_ionfrac_weighted_phase_diagrams('x','_nref11_','RD0020','forcedres',ds.current_redshift)
#make_ionfrac_weighted_phase_diagrams('x','_nref11n_nref10f_refine200kpc_z4to2_','RD0020','forcedres',ds.current_redshift)
#create_emission_frbs()
#make_small_emission_gif_plots()
#create_phys_emis_weight_frbs()
#make_emission_gif_plots()
cdf_plot_loop()
#plot_SB_profiles_all_lines(box_width)
