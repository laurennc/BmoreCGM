import matplotlib
matplotlib.use('Agg')

import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from emission_functions import *

import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec


def get_refine_box(ds, zsnap, track):
    ## find closest output, modulo not updating before printout
    diff = track['col1'] - zsnap
    this_loc = track[np.where(diff == np.min(diff[np.where(diff > 1.e-6)]))]
    print "using this loc:", this_loc
    x_left = this_loc['col2'][0]
    y_left = this_loc['col3'][0]
    z_left = this_loc['col4'][0]
    x_right = this_loc['col5'][0]
    y_right = this_loc['col6'][0]
    z_right = this_loc['col7'][0]
    refine_box_center = [0.5*(x_left+x_right), 0.5*(y_left+y_right), 0.5*(z_left+z_right)]
    refine_box = ds.r[x_left:x_right, y_left:y_right, z_left:z_right]
    refine_width = np.abs(x_right - x_left)
    return refine_box, refine_box_center, refine_width


sns.set_style("white", {'axes.grid' : False})
colors1 = plt.cm.Greys(np.linspace(0., 0.8, 192))
sns_cmap2 = sns.blend_palette(('Pink','DeepPink','#1E90FF','DarkTurquoise','#50A638'), n_colors=40,as_cmap=True)
colors2 = sns_cmap2(np.linspace(0.,1,64))
colors = np.vstack((colors1, colors2))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

RDs = np.arange(19,140,1)

track_name = baseREF+"/halo_track"
track = Table.read(track_name, format='ascii')
track.sort('col1')

for RD in RDs:
    if RD < 100:
        RD = 'RD00'+str(RD)
    else:
        RD = 'RD0'+str(RD)

    print '###############################'
    print '########## '+RD+' #############'
    print '###############################'

    #baseREF = "/lou/s2m/mpeeples/halo_008508/nref10_refine200kpc_z4to2"
    fnREF = baseREF+"/"+RD+"/"+RD
    #baseNAT = "/lou/s1z/jtumlins/Halos/halo_008508/natural/nref10"
    fnNAT = baseNAT+"/"+RD+"/"+RD
    fileout = 'EMIS_'+RD+'.pdf'

    dsREF = yt.load(fnREF)
    dsNAT = yt.load(fnNAT)
    rbR, rb_centerR, rb_widthR = get_refine_box(dsREF, dsREF.current_redshift, track)
    rbN, rb_centerN, rb_widthN = get_refine_box(dsNAT, dsNAT.current_redshift, track)

    spN = dsNAT.sphere(rb_centerN,(30.,'kpc'))
    idx = np.where(spN['Density'] == spN['Density'].max())[0][0]
    posNAT = dsNAT.arr([float(spN['x'][idx].value),float(spN['y'][idx].value),float(spN['z'][idx].value)],'code_length')
    spR = dsNAT.sphere(rb_centerR,(30.,'kpc'))
    idx = np.where(spR['Density'] == spR['Density'].max())[0][0]
    posREF = dsREF.arr([float(spR['x'][idx].value),float(spR['y'][idx].value),float(spR['z'][idx].value)],'code_length')


    fig,ax = plt.subplots(2,3)#,sharex=True,sharey=True)
    fig.set_size_inches(20,13)


    ## ok what I should do is find the most dens position, make that the center, then do
    ## +/- 30 kpc
    LR,RR = posREF-dsREF.quan(30.,'kpc').to('code_length'),posREF+dsREF.quan(30.,'kpc').to('code_length')
    LN,RN = posNAT-dsNAT.quan(30.,'kpc').to('code_length'),posNAT+dsNAT.quan(30.,'kpc').to('code_length')

    rbR = dsREF.r[LR[0]:RR[0],LR[1]:RR[1],LR[2]:RR[2]]
    rbN = dsNAT.r[LN[0]:RN[0],LN[1]:RN[1],LN[2]:RN[2]]

    ### SiIV Ref ###
    emisSiR = yt.ProjectionPlot(dsREF,'z',('gas','Emission_SiIV'),data_source=rbR,center=posREF,width=(60,'kpc'))
    emisSiR.set_zlim(('gas','Emission_SiIV'),10**-5,10**3)
    frbSiR = np.log10(emisSiR.frb[('gas','Emission_SiIV')]/(1.+dsREF.current_redshift)**4)
    cPickle.dump(frbSiR,open('frbz_nref11n_nref10f_200kpc_'+RD+'_Emission_SiIV_ytres.cpkl','wb'),protcol=-1)

    im = ax[0,0].imshow(frbSiR,vmin=-5.,vmax=3.,interpolation='none',cmap=mymap,origin='lower')
    im.axes.set_xticks([])
    im.axes.set_yticks([])
    ax[0,0].text(0.95, 0.95, 'SiIV',
            verticalalignment='top', horizontalalignment='right',
            transform=ax[0,0].transAxes,
            color='black', fontsize=30)

    ### CIV Ref ###
    emisCR = yt.ProjectionPlot(dsREF,'z',('gas','Emission_CIV'),data_source=rbR,center=posREF,width=(60,'kpc'))
    emisCR.set_zlim(('gas','Emission_CIV'),10**-5,10**3)
    frbCR = np.log10(emisCR.frb[('gas','Emission_CIV')]/(1.+dsREF.current_redshift)**4)
    cPickle.dump(frbCR,open('frbz_nref11n_nref10f_200kpc_'+RD+'_Emission_CIV_ytres.cpkl','wb'),protcol=-1)

    im2 = ax[0,1].imshow(frbCR,vmin=-5.,vmax=3.,interpolation='none',cmap=mymap,origin='lower')
    im2.axes.set_xticks([])
    im2.axes.set_yticks([])
    ax[0,1].text(0.95, 0.95, 'CIV',
            verticalalignment='top', horizontalalignment='right',
            transform=ax[0,1].transAxes,
            color='black', fontsize=30)

    ### OVI Ref ###
    emisOR = yt.ProjectionPlot(dsREF,'z',('gas','Emission_OVI'),data_source=rbR,center=posREF,width=(60,'kpc'))
    emisOR.set_zlim(('gas','Emission_OVI'),10**-5,10**3)
    frbOR = np.log10(emisOR.frb[('gas','Emission_OVI')]/(1.+dsREF.current_redshift)**4)
    cPickle.dump(frbOR,open('frbz_nref11n_nref10f_200kpc_'+RD+'_Emission_OVI_ytres.cpkl','wb'),protcol=-1)

    im3 = ax[0,2].imshow(frbOR,vmin=-5.,vmax=3.,interpolation='none',cmap=mymap,origin='lower')
    im3.axes.set_xticks([])
    im3.axes.set_yticks([])
    ax[0,2].text(0.95, 0.95, 'OVI',
            verticalalignment='top', horizontalalignment='right',
            transform=ax[0,2].transAxes,
            color='black', fontsize=30)

    ## SiIV Nat ###
    emisSiN = yt.ProjectionPlot(dsNAT,'z',('gas','Emission_SiIV'),data_source=rbN,center=posNAT,width=(60,'kpc'))
    emisSiN.set_zlim(('gas','Emission_SiIV'),10**-5,10**3)
    frbSiN = np.log10(emisSiN.frb[('gas','Emission_SiIV')]/(1.+dsNAT.current_redshift)**4)
    cPickle.dump(frbSiN,open('frbz_nref11n_'+RD+'_Emission_SiIV_ytres.cpkl','wb'),protcol=-1)

    im = ax[1,0].imshow(frbSiN,vmin=-5.,vmax=3.,interpolation='none',cmap=mymap,origin='lower')
    im.axes.set_xticks([])
    im.axes.set_yticks([])


    ## CIV Nat ###
    emisCN = yt.ProjectionPlot(dsNAT,'z',('gas','Emission_CIV'),data_source=rbN,center=posNAT,width=(60,'kpc'))
    emisCN.set_zlim(('gas','Emission_CIV'),10**-5,10**3)
    frbCN = np.log10(emisCN.frb[('gas','Emission_CIV')]/(1.+dsNAT.current_redshift)**4)
    cPickle.dump(frbCN,open('frbz_nref11n_'+RD+'_Emission_CIV_ytres.cpkl','wb'),protcol=-1)

    im2 = ax[1,1].imshow(frbCN,vmin=-5.,vmax=3.,interpolation='none',cmap=mymap,origin='lower')
    im2.axes.set_xticks([])
    im2.axes.set_yticks([])

    ### OVI Nat ###
    emisON = yt.ProjectionPlot(dsNAT,'z',('gas','Emission_OVI'),data_source=rbN,center=posNAT,width=(60,'kpc'))
    emisON.set_zlim(('gas','Emission_OVI'),10**-5,10**3)
    frbON = np.log10(emisON.frb[('gas','Emission_OVI')]/(1.+dsNAT.current_redshift)**4)
    cPickle.dump(frbOR,open('frbz_nref11n_'+RD+'_Emission_OVI_ytres.cpkl','wb'),protcol=-1)

    im3 = ax[1,2].imshow(frbON,vmin=-5.,vmax=3.,interpolation='none',cmap=mymap,origin='lower')
    im3.axes.set_xticks([])
    im3.axes.set_yticks([])

    plt.tight_layout()
    fig.subplots_adjust(bottom=0.065)


    cbaxes2 = fig.add_axes([0.363, 0.035, 0.273, 0.02]) #0.015
    cb2 = plt.colorbar(im2, cax = cbaxes2,orientation='horizontal')
    axcb = cb2.ax
    text = axcb.yaxis.label

    plt.savefig(fileout)
    plt.close()
