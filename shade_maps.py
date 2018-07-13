
""" a module for datashader renders of phase diagrams"""
from functools import partial
import datashader as dshader
from datashader.utils import export_image
import datashader.transfer_functions as tf
import holoviews as hv
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import yt
import trident
import numpy as np
from astropy.table import Table
from get_refine_box import get_refine_box as grb
from consistency import ion_frac_color_key, phase_color_key, metal_color_key, axes_label_dict, logfields
from holoviews.operation.datashader import datashade, aggregate
from holoviews import Store
hv.extension('matplotlib')

def prep_dataset(fname, trackfile, ion_list=['H I'], region='trackbox'):
    """prepares the dataset for rendering by extracting box or sphere"""
    data_set = yt.load(fname)

    trident.add_ion_fields(data_set, ions=ion_list)
    for ion in ion_list:
        print("Added ion "+ion+" into the dataset.")

    track = Table.read(trackfile, format='ascii')
    track.sort('col1')
    refine_box, refine_box_center, refine_width = \
            grb(data_set, data_set.current_redshift, track)
    print('Refine box corners: ', refine_box)
    print('           center : ', refine_box_center)

    if region == 'trackbox':
        all_data = refine_box
    elif region == 'sphere':
        sph = data_set.sphere(center=refine_box_center, radius=(500, 'kpc'))
        all_data = sph
    else:
        print("your region is invalid!")

    return all_data, refine_box, refine_width

def categorize_by_temp(temp):
    """ define the temp category strings"""
    phase = np.chararray(np.size(temp), 4)
    phase[temp < 9.] = 'hot'
    phase[temp < 6.] = 'warm'
    phase[temp < 5.] = 'cool'
    phase[temp < 4.] = 'cold'
    return phase

def categorize_by_fraction(f_ion):
    """ define the ionization category strings"""
    frac = np.chararray(np.size(f_ion), 4)
    frac[f_ion > -10.] = 'all'
    frac[f_ion > 0.01] = 'low' # yellow
    frac[f_ion > 0.1] = 'med'  # orange
    frac[f_ion > 0.2] = 'high' # red
    return frac

def categorize_by_metallicity(metallicity):
    """ define the metallicity category strings"""
    metal_label = np.chararray(np.size(metallicity), 5)
    metal_label[metallicity < 10.] = 'high'
    metal_label[metallicity < 0.005] = 'solar'
    metal_label[metallicity < 0.000001] = 'low'
    metal_label[metallicity < 0.0000001] = 'poor'
    return metal_label

def scale_lvec(lvec):
    lvec[lvec > 0.] = (np.log10(lvec[lvec > 0.]) - 25.) / 8.
    lvec[lvec < 0.] = (-1. * np.log10(-1.*lvec[lvec < 0.]) + 25.) / 8.
    return lvec

def prep_dataframe(all_data, refine_box, refine_width, field1, field2):
    """ add fields to the dataset, create dataframe for rendering
        The enzo fields x, y, z, temperature, density, cell_vol, cell_mass,
        and metallicity will always be included, others will be included
        if they are requested as fields. """

    # obtain fields that we'll use no matter what the input fields.
    density = all_data['density']
    cell_mass = all_data['cell_volume'].in_units('kpc**3') * density.in_units('Msun / kpc**3')
    cell_size = np.array(all_data["cell_volume"])**(1./3.)

    x = all_data['x'].ndarray_view() + cell_size * (np.random.rand(np.size(cell_size)) * 2. - 1.)
    y = all_data['y'].ndarray_view() + cell_size * (np.random.rand(np.size(cell_size)) * 2. - 1.)
    z = all_data['z'].ndarray_view() + cell_size * (np.random.rand(np.size(cell_size)) * 2. - 1.)
    x = (x - refine_box.center[0].ndarray_view()) / (refine_width/2.)
    y = (y - refine_box.center[1].ndarray_view()) / (refine_width/2.)
    z = (z - refine_box.center[2].ndarray_view()) / (refine_width/2.)

    density = np.log10(density)
    temperature = np.log10(all_data['temperature'])
    mass = np.log10(cell_mass)
    phase = categorize_by_temp(temperature)
    metal = categorize_by_metallicity(all_data['metallicity'])

    # build data_frame with mandatory fields
    data_frame = pd.DataFrame({'x':x, 'y':y, 'z':z, 'temperature':temperature, \
                               'density':density, 'cell_mass': mass, \
                               'phase':phase, 'metal':metal})
    data_frame.phase = data_frame.phase.astype('category')
    data_frame.metal = data_frame.metal.astype('category')

    # now add the optional fields
    print("you have requested fields ", field1, field2)

    # add those two fields
    #logfields = ('density', 'temperature', 'entropy', 'O_p5_ion_fraction',
    #            'C_p3_ion_fraction', 'Si_p3_ion_fraction', 'O_p5_number_density',
    #            'C_p3_number_density', 'Si_p3_number_density')
    logfields = ('density','temperature','entropy','O_p5_number_density',
                 'C_p3_number_density','Si_p3_number_density')
    if field1 not in data_frame.columns:
        print("Did not find "+field1+" in the dataframe, will add it.")
        if field1 in logfields:
            print("Field 1, "+field1+" is a log field.")
            data_frame[field1] = np.log10(all_data[field1])
        else:
            data_frame[field1] = all_data[field1]
    if field2 not in data_frame.columns:
        print("Did not find "+field2+" in the dataframe, will add it.")
        if field2 in logfields:
            print("Field 2, "+field2+" is a log field.")
            data_frame[field2] = np.log10(all_data[field2])
        else:
            data_frame[field2] = all_data[field2]

    if ((field1[0:4] == 'O_p5') or (field2[0:4] == 'O_p5')):
        data_frame['o6frac'] = categorize_by_fraction(np.array(data_frame['O_p5_ion_fraction']))
        data_frame.o6frac = data_frame.o6frac.astype('category')
    if ((field1[0:4] == 'C_p3') or (field2[0:4] == 'C_p3')):
        data_frame['c4frac'] = categorize_by_fraction(np.array(data_frame['C_p3_ion_fraction']))
        data_frame.c4frac = data_frame.c4frac.astype('category')
    if ((field1[0:5] == 'Si_p3') or (field2[0:4] == 'Si_p3')):
        data_frame['si4frac'] = categorize_by_fraction(np.array(data_frame['Si_p3_ion_fraction']))
        data_frame.si4frac = data_frame.si4frac.astype('category')

    return data_frame


def wrap_axes(filename, field1, field2, ranges):
    """intended to be run after render_image, take the image and wraps it in
        axes using matplotlib and so offering full customization."""

    img = mpimg.imread(filename+'.png')
    #img2 = np.flip(img,0)
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1, 0.1, 0.88, 0.88])
    #ax.imshow(img2)
    ax.imshow(img)

    xtext = ax.set_xlabel(axes_label_dict[field1], fontname='Arial', fontsize=20)
    ax.set_xticks(np.arange((ranges[0][1] - ranges[0][0]) + 1.) * 1000. / (ranges[0][1] - ranges[0][0]))
    ax.set_xticklabels([ str(int(s)) for s in np.arange((ranges[0][1] - ranges[0][0]) + 1.) +  ranges[0][0] ], fontname='Arial', fontsize=20)

    ytext = ax.set_ylabel(axes_label_dict[field2], fontname='Arial', fontsize=20)
    ax.set_yticks(np.arange((ranges[1][1] - ranges[1][0]) + 1.) * 1000. / (ranges[1][1] - ranges[1][0]))
    print [ str(int(s)) for s in np.arange((ranges[1][1] - ranges[1][0]) + 1.) +  ranges[1][0] ]
    ax.set_yticklabels([ str(int(s)) for s in np.arange((ranges[1][0] - ranges[1][1]) + 1.) +  ranges[1][0] ], fontname='Arial', fontsize=20)

    plt.savefig(filename)


def render_image(frame, field1, field2, count_cat, x_range, y_range, filename):
    """ renders density and temperature 'Phase' with linear aggregation"""

    export = partial(export_image, background='white', export_path="./")
    cvs = dshader.Canvas(plot_width=1000, plot_height=1000,
                         x_range=x_range, y_range=y_range)
    agg = cvs.points(frame, field1, field2, dshader.count_cat(count_cat))

    if 'frac' in count_cat:
        color_key = ion_frac_color_key
    elif 'phase' in count_cat:
        color_key = phase_color_key
    elif 'metal' in count_cat:
        color_key = metal_color_key

    img = tf.shade(agg, color_key=color_key, how='eq_hist')
    export(img, filename)
    return img



def drive(fname, trackfile, ion_list=['H I', 'C IV', 'Si IV', 'O VI']):
    """this function drives datashaded phase plots"""

    all_data, refine_box, refine_width = \
        prep_dataset(fname, trackfile, ion_list=ion_list, region='sphere')

    data_frame = prep_dataframe(all_data, refine_box, refine_width, field1, field2)

    for ion in ['o6', 'c4', 'si4']:
        render_image(data_frame, 'density', 'temperature', ion+'frac',
                     (-31, -20), (2,8), 'RD0020_phase_'+ion)
        render_image(data_frame, 'x', 'y', ion+'frac',
                     (-3,3), (-3,3), 'RD0020_proj_'+ion)

    render_image(data_frame, 'temperature', 'logf_o6', 'phase', (2, 8), (-5, 0), 'RD0020_ionfrac')
    render_image(data_frame, 'density', 'temperature', 'phase', (-31, -20), (2, 8), 'RD0020_phase')
    render_image(data_frame, 'x', 'y', 'phase', (-3,3), (-3,3), 'RD0020_proj')
    render_image(data_frame, 'x', 'mass', 'phase', (-3.1, 3.1), (-1, 8), 'RD0020_mass')
    render_image(data_frame, 'x', 'lz', 'phase', (-1.1, 1.1), (-1.1, 1.1), 'RD0020_lz')
