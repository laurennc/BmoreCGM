import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import cPickle
from emission_functions import *
from astropy.cosmology import WMAP9 as cosmo
#import seaborn as sns
#sns.set_style("whitegrid", {'axes.grid' : False})

base = "/Users/dalek/data/Molly/nref11n_nref10f_refine200kpc_z4to2"
fn = base+"/RD0027/RD0027"
track_name = base+"/halo_track_z2_z1"
args = fn.split('/')
lines = ['LyAlpha','CIII_977','OVI']

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
    ds = yt.load(fn)
    track = Table.read(track_name, format='ascii')
    track.sort('col1')
    rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)
    lines = ['LyAlpha','CIII_977','OVI']

    ang_res = 1.0 ## arcsecond
    dx = (ang_res / cosmo.arcsec_per_kpc_proper(ds.current_redshift)).value
    dx = ds.arr(dx,'kpc').in_units('code_length')

    for line in lines:
        field = 'Emission_'+line
        print line
        for index in 'xyz':
                num_cells = np.ceil(rb_width/dx)
                fileout = args[-3]+'_'+args[-2]+'_'+field+'_'+index+'_luvoir'
                obj = ds.proj(('gas',field),index,data_source=rb)
                frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
                fileout = 'luvoir/frb'+index+'_'+args[-3]+'_'+args[-2]+'_'+field+'_1arcsec_luvoir.cpkl'
                cPickle.dump(frb[('gas',field)],open(fileout,'wb'),protocol=-1)
    return

def create_phys_emis_weight_frbs():
    ds = yt.load(fn)
    track = Table.read(track_name, format='ascii')
    track.sort('col1')
    rb,rb_center,rb_width = get_refine_box(ds,ds.current_redshift,track)
    lines = ['CIII_977','LyAlpha','OVI']

    ang_res = 1 ## arcsecond
    dx = (ang_res / cosmo.arcsec_per_kpc_proper(ds.current_redshift)).value
    dx = ds.arr(dx,'kpc').in_units('code_length')

    for line in lines:
        field = 'Emission_'+line
        print line

        for index in 'xyz':
            print index
            num_cells = np.ceil(rb_width/dx)
            obj = ds.proj('H_nuclei_density',index,data_source=rb,weight_field=field)
            frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
            fileout = 'luvoir/frb'+index+'_'+args[-3]+'_'+args[-2]+'_hden_'+field+'_1arcsec_luvoir.cpkl'
            cPickle.dump(frb['H_nuclei_density'],open(fileout,'wb'),protocol=-1)

            obj = ds.proj('temperature',index,data_source=rb,weight_field=field)
            frb = obj.to_frb((rb_width,'code_length'),(num_cells,num_cells),center=rb_center)
            fileout = 'luvoir/frb'+index+'_'+args[-3]+'_'+args[-2]+'_temp_'+field+'_1arcsec_luvoir.cpkl'
            cPickle.dump(frb['temperature'],open(fileout,'wb'),protocol=-1)

    return


def make_luvoir_plot(RD,SBlim):
    lines = ['LyAlpha','CIII_977','OVI']
    ds = yt.load(fn)
    for line in lines:
        field = 'Emission_'+line
        filein = 'luvoir/frbs/frbx_nref11n_nref10f_refine200kpc_z4to2_'+RD+'_'+field+'_1arcsec_luvoir.cpkl'
        frb = cPickle.load(open(filein,'rb'))
        frb = np.log10(frb/(1.+ds.current_redshift)**4)

        test = np.ma.masked_where((frb < SBlim),frb)
        plt.imshow(frb,cmap='Greys',vmin=-5,vmax=3)
        plt.colorbar()
        plt.imshow(test,cmap='GnBu',vmin=1,vmax=3)
        plt.colorbar()

        plt.savefig('z1_x_'+field+'_1arcsec_luvoir.pdf')
        plt.close()
    return


def hden_temp_hist(RD,line,redshift):
    hden_file = '/Users/dalek/repos/BmoreCGM/luvoir/frbs/frbx_nref11n_nref10f_refine200kpc_z4to2_'+RD+'_hden_Emission_'+line+'_1arcsec_luvoir.cpkl'
    temp_file = '/Users/dalek/repos/BmoreCGM/luvoir/frbs/frbx_nref11n_nref10f_refine200kpc_z4to2_'+RD+'_temp_Emission_'+line+'_1arcsec_luvoir.cpkl'
    emis_file = '/Users/dalek/repos/BmoreCGM/luvoir/frbs/frbx_nref11n_nref10f_refine200kpc_z4to2_'+RD+'_Emission_'+line+'_1arcsec_luvoir.cpkl'

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

    return hist,avgvals,maxvals

def make_weighted_phase_diagrams(RD):
    ds = yt.load(fn)
    for line in lines:
        hist,avgvals,maxvals = hden_temp_hist(RD,line,ds.current_redshift)

        fig,ax = plt.subplots(1,1)
        im1 = ax.imshow(np.log10(hist.T),extent=[-6,1,3.5,6.5],
                            interpolation='nearest',origin='lower',cmap='Blues')
        add_grid(ax)
        ax.set_xlim(-6,1.5)
        ax.set_ylim(3.5,6.5)
        ax.set_xlabel('log(Hydrogen Number Density) [cm^-3]')
        ax.set_ylabel('log(Temperature) [K]')
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        plt.colorbar(im1)
        plt.savefig('z1_x_phase_weighted_Emission_'+line+'_1arcsec_luvoir.pdf')
        plt.close()

        fig,ax = plt.subplots(1,1)
        #im2 = ax.imshow(np.log10(avgvals.T),extent=[-6,1,3.5,6.5],
        #                    interpolation='nearest',origin='lower',
        #                    vmin=-3,vmax=6,cmap='jet')

        avgvals = np.log10(avgvals)
        test = np.ma.masked_where((avgvals < 1),avgvals)
        test2 = np.ma.masked_where((avgvals >= 1),avgvals)

        im2a = ax.imshow(test2.T,cmap='Greys',vmin=-5,vmax=3,origin='lower',
                         extent=[-6,1,3.5,6.5],interpolation='nearest')
        im2b = ax.imshow(test.T,cmap='GnBu',vmin=1,vmax=6,origin='lower',
                         extent=[-6,1,3.5,6.5],interpolation='nearest')
        add_grid(ax)
        ax.set_xlim(-6,1.5)
        ax.set_ylim(3.5,6.5)
        ax.set_xlabel('log(Hydrogen Number Density) [cm^-3]')
        ax.set_ylabel('log(Temperature) [K]')
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        plt.colorbar(im2b)
        plt.savefig('z1_x_phase_weighted_Emission_'+line+'_avgEmis_1arcsec_luvoir.pdf')
        plt.close()

        fig,ax = plt.subplots(1,1)
        im3 = ax.imshow(np.log10(maxvals.T),extent=[-6,1,3.5,6.5],
                            interpolation='nearest',origin='lower',
                            vmin=-3,vmax=6,cmap='jet')
        add_grid(ax)
        ax.set_xlim(-6,1.5)
        ax.set_ylim(3.5,6.5)
        ax.set_xlabel('log(Hydrogen Number Density) [cm^-3]')
        ax.set_ylabel('log(Temperature) [K]')
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        plt.colorbar(im3)
        plt.savefig('z1_x_phase_weighted_Emission_'+line+'_maxEmis_1arcsec_luvoir.pdf')
        plt.close()

    return

#create_emission_frbs()
#create_phys_emis_weight_frbs()
#make_luvoir_plot('RD0027',1.)
make_weighted_phase_diagrams('RD0027')
