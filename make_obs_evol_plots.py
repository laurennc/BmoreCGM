import matplotlib
matplotlib.use('Agg')

import yt
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
import trident

def get_refine_box(ds, zsnap, track):
    ## find closest output, modulo not updating before printout
    #diff = track['col1'] - zsnap
    #this_loc = track[np.where(diff == np.min(diff[np.where(diff > 1.e-6)]))]
    #print "using this loc:", this_loc
    x_left = np.interp(zsnap, track['col1'], track['col2'])
    y_left = np.interp(zsnap, track['col1'], track['col3'])
    z_left = np.interp(zsnap, track['col1'], track['col4'])
    x_right = np.interp(zsnap, track['col1'], track['col5'])
    y_right = np.interp(zsnap, track['col1'], track['col6'])
    z_right = np.interp(zsnap, track['col1'], track['col7'])
    refine_box_center = [0.5*(x_left+x_right), 0.5*(y_left+y_right), 0.5*(z_left+z_right)]
    refine_box = ds.r[x_left:x_right, y_left:y_right, z_left:z_right]
    refine_width = np.abs(x_right - x_left)
    return refine_box, refine_box_center, refine_width

fontrc ={'fontname':'Osaka','fontsize':30}

baseREF = "/lou/s2m/mpeeples/halo_008508/nref11n/nref11n_nref10f_refine200kpc_z4to2"
baseNAT = "/lou/s2m/mpeeples/halo_008508/nref11n/natural"
track_name = "/u/lcorlies/halo_track"
track = Table.read(track_name, format='ascii')
track.sort('col1')

#RDs = np.arange(19,355,1)
RDs = [20,21]

h1_color_map = sns.blend_palette(("white","#ababab","#565656","black","#4575b4","#984ea3","#d73027","darkorange","#ffe34d"), as_cmap=True)
h1_proj_min = 11
h1_proj_max = 21

o6_color_map = sns.blend_palette(("white","black","#4daf4a","#4575b4","#984ea3","#d73027","darkorange"), as_cmap=True)
o6_min = 11
o6_max = 15

si3_color_map = "magma"
si3_min = 11
si3_max = 16

for RD in RDs:
    if RD < 100:
        RD = 'DD00'+str(RD)
    else:
        RD = 'DD0'+str(RD)

    print '###############################'
    print '########## '+RD+' #############'
    print '###############################'

    fnREF = baseREF+"/"+RD+"/"+RD
    fnNAT = baseNAT+"/"+RD+"/"+RD
    fileout = '/u/lcorlies/OBSplots/OBS_'+RD+'.pdf'

    dsREF = yt.load(fnREF)
    dsNAT = yt.load(fnNAT)
    trident.add_ion_fields(dsREF,ions=['O VI','Si III'])
    trident.add_ion_fields(dsNAT,ions=['O VI','Si III'])
    rbR, rb_centerR, rb_widthR = get_refine_box(dsREF, dsREF.current_redshift, track)
    rbN, rb_centerN, rb_widthN = get_refine_box(dsNAT, dsNAT.current_redshift, track)
    rb_widthR,rb_widthN = rb_widthR*0.98,rb_widthN*0.98

    fig,ax = plt.subplots(2,3)#,sharex=True,sharey=True)
    fig.set_size_inches(20,12)

    ### HI Ref ###
    print 'HI REF'
    #val,pos = dsREF.find_max('dark_matter_density')
    HIR = yt.ProjectionPlot(dsREF,'x','H_number_density',center=rb_centerR,data_source=rbR,
                              width=(rb_widthR,'code_length'))
    HIR.set_zlim('H_number_density',10**h1_proj_min,10**h1_proj_max)
    frbHR = np.log10(HIR.frb['H_number_density'])

    im = ax[0,0].imshow(frbHR,vmin=h1_proj_min,vmax=h1_proj_map,interpolation='none',cmap=h1_color_map,origin='lower')
    im.axes.set_xticks([])
    im.axes.set_yticks([])
    ax[0,0].text(0.95, 0.95, 'HI',
                verticalalignment='top', horizontalalignment='right',
                transform=ax[0,0].transAxes,color='white', fontsize=30)

    ### SiIII Ref ###
    print 'SiIII REF'
    SiIIIR = yt.ProjectionPlot(dsREF,'x','Si_p2_number_density',center=rb_centerR,data_source=rbR,
                              width=(rb_widthR,'code_length'))
    SiIIIR.set_zlim('Si_p2_number_density',10**si3_min,10**si3_max)
    frbSR = np.log10(SiIIIR.frb['Si_p2_number_density'])

    im2 = ax[0,1].imshow(frbSR,vmin=si3_min,vmax=si3_max,interpolation='none',cmap=si3_color_map,origin='lower')
    im2.axes.set_xticks([])
    im2.axes.set_yticks([])
    ax[0,1].text(0.95, 0.95, 'SiIII',
                verticalalignment='top', horizontalalignment='right',
                transform=ax[0,1].transAxes,color='black', fontsize=30)

    ### OVI Ref ###
    print 'OVI REF'
    OVIR = yt.ProjectionPlot(dsREF,'x','O_p5_number_density',center=rb_centerR,data_source=rbR,
                               width=(rb_widthR,'code_length'))
    OVIR.set_zlim('O_p5_number_density',10**o6_min,10**o6_max)
    frbOR = np.log10(OVIR.frb['O_p5_number_density'])

    im3 = ax[0,2].imshow(frbOR,vmin=o6_min,vmax=o6_max,interpolation='none',cmap=o6_color_map,origin='lower')
    im3.axes.set_xticks([])
    im3.axes.set_yticks([])
    ax[0,2].text(0.95, 0.95, 'OVI',
                verticalalignment='top', horizontalalignment='right',
                transform=ax[0,2].transAxes,color='white', fontsize=30)

    ## HI  Nat ###
    #val,pos = dsNAT.find_max('dark_matter_density')
    print 'HI NAT'
    HIN = yt.ProjectionPlot(dsNAT,'x','H_number_density',center=rb_centerN,data_source=rbN,
                              width=(rb_widthN,'code_length'))
    HIN.set_zlim('H_number_density',10**h1_proj_min,10**h1_proj_max)
    frbHN = np.log10(HIN.frb['H_number_density'])

    im = ax[1,0].imshow(frbHN,vmin=h1_proj_min,vmax=h1_proj_max,interpolation='none',cmap=h1_color_map,origin='lower')
    im.axes.set_xticks([])
    im.axes.set_yticks([])


    ## SiIII Nat ###
    print 'SiIII NAT'
    SiIIIN = yt.ProjectionPlot(dsNAT,'x','Si_p2_number_density',center=rb_centerN,data_source=rbN,
                              width=(rb_widthN,'code_length'))
    SiIIIN.set_zlim('Si_p2_number_density',10**si3_min,10**si3_max)
    frbSN = np.log10(SiIIIN.frb['Si_p2_number_density'])

    im2 = ax[1,1].imshow(frbSN,vmin=si3_min,vmax=si3_max,interpolation='none',cmap=si3_color_map,origin='lower')
    im2.axes.set_xticks([])
    im2.axes.set_yticks([])

    ### OVI Nat ###
    print 'OVI NAT'
    OVIN = yt.ProjectionPlot(dsNAT,'x','O_p5_number_density',center=rb_centerN,data_source=rbN,
                               width=(rb_widthN,'code_length'))
    OVIN.set_zlim('O_p5_number_density',10**o6_min,10**o6_max)
    frbON = np.log10(OVIN.frb['O_p5_number_density'])
    im3 = ax[1,2].imshow(frbON,vmin=o6_min,vmax=o6_max,interpolation='none',cmap=o6_color_map,origin='lower')
    im3.axes.set_xticks([])
    im3.axes.set_yticks([])

    ## Reformat ##
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.065)

    ## Colorbars ##
    print 'COLORBARS'
    cbaxes = fig.add_axes([0.0328, 0.03, 0.273, 0.02])
    cb = plt.colorbar(im, cax = cbaxes,orientation='horizontal')
    axcb = cb.ax
    text = axcb.yaxis.label
    font = mpl.font_manager.FontProperties(family='Osaka', size=20)
    text.set_font_properties(font)

    cbaxes2 = fig.add_axes([0.363, 0.03, 0.273, 0.02])
    cb2 = plt.colorbar(im2, cax = cbaxes2,orientation='horizontal')
    axcb = cb2.ax
    text = axcb.yaxis.label
    font = mpl.font_manager.FontProperties(family='Osaka', size=20)
    text.set_font_properties(font)

    cbaxes3 = fig.add_axes([0.689, 0.03, 0.273, 0.02])
    cb3 = plt.colorbar(im3, cax = cbaxes3,orientation='horizontal')
    axcb = cb3.ax
    text = axcb.yaxis.label
    font = mpl.font_manager.FontProperties(family='Osaka', size=20)
    text.set_font_properties(font)

    plt.savefig(fileout)
    plt.close()
