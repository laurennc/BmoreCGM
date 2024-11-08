"""
useful stuff for cloud analysis JT090618
useful stuff added LC010219
"""

from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
import numpy as np
import trident
import cPickle
import glob
import os
import yt
import copy
import re

def reduce_ion_vector(vx, ion):
    """ this function takes in two vectors for velocity and ionization
        fraction and chunks the ionization fraction into a uniform velocity
        grid. JT 082018"""
    v = np.arange(3001) - 1500
    ion_hist = v * 0.
    index = np.clip(np.around(vx) + 1500, 0, 2999)
    for i in np.arange(np.size(ion)):
        ion_hist[int(index[i])] = ion_hist[int(
            index[i])] + ion[int(i)]

    return v, ion_hist


def get_fion_threshold(ion_to_use, coldens_fraction):
    cut = 0.999
    total = np.sum(ion_to_use)
    ratio = 0.001
    while ratio < coldens_fraction:
        part = np.sum(
            ion_to_use[ion_to_use > cut * np.max(ion_to_use)])
        ratio = part / total
        cut = cut - 0.001

    threshold = cut * np.max(ion_to_use)
    number_of_cells_above_threshold = np.size(
        np.where(ion_to_use > threshold))

    return threshold, number_of_cells_above_threshold


def get_sizes(ray_df, species, x, axis_to_use, ion_to_use, coldens_threshold):

    threshold, number_of_cells = get_fion_threshold(
        ion_to_use, coldens_threshold)

    cell_mass = np.array(ray_df['cell_mass'])
    dx = np.array(ray_df['dx'])
    axis_velocity = np.array(ray_df[axis_to_use+'-velocity'])

    ion_density = copy.deepcopy(ion_to_use)

    # insert a cloud flag vector that IDs the cloud
    cloud_flag = np.zeros(np.size(dx), dtype=np.int8)

    indexsizes = []
    kpcsizes = []
    column_densities = []
    masses = []
    centers = []
    velocities = []
    indices = []
    xs = []
    for m in np.arange(100):  # can find up to 100 peaks
        i = np.squeeze(np.where(np.array(ion_to_use) > threshold))

        if np.size(i) >= 1:
            startindex = np.min(i)
            f = ion_to_use[startindex]
            index = startindex
            ion_to_use[startindex] = 0.0
            sum_mass = cell_mass[startindex]
            sum_coldens = ion_density[startindex] * dx[index]
            count = 0
            velsum = 0.
            while (f > threshold) and (index < np.size(x)-1):
                count += 1
                cloud_flag[index] = m+1 # place cloud number in flag vector
                if (count > 10000):
                    os.sys.exit('stuck in the size finding loop')
                index += 1
                if index == np.size(x):  # this means we're at the edge
                    index = np.size(x)-1
                    f = 0.0
                else:
                    f = ion_to_use[index]
                    ion_to_use[index] = 0.0
                    sum_mass = sum_mass + cell_mass[index]
                    velsum = velsum + \
                        cell_mass[index] * axis_velocity[index]
                    sum_coldens = sum_coldens + \
                        ion_density[index] * dx[index]

            x_coord = x[startindex:index]
            ion_d = ion_density[startindex:index]
            ion_center = np.sum(x_coord * ion_d) / np.sum(ion_d)

            indexsizes.append(index - startindex)
            kpcsizes.append(x[startindex]-x[index])
            column_densities.append(sum_coldens)
            masses.append(sum_mass)
            # should end up with mass-weighted velocity along LOS
            velocities.append(velsum / sum_mass)
            centers.append(ion_center)
            indices.append(index)
            xs.append(x[index])

    size_dict = {'coldens_threshold': coldens_threshold}
    size_dict[species+'_xs'] = xs
    size_dict[species+'_indices'] = indices
    size_dict[species+'_kpcsizes'] = kpcsizes
    size_dict[species+'_indexsizes'] = indexsizes
    size_dict[species+'_coldens'] = column_densities
    size_dict[species+'_n_cells'] = number_of_cells
    size_dict[species+'_cell_masses'] = masses
    size_dict[species+'_centers'] = centers
    size_dict[species+'_velocities'] = velocities

    ray_df[species+'_cloud_flag'] = cloud_flag

    return size_dict

def get_trident_ray(ds, ray_start, ray_end, line_list, **kwargs):
    '''
    input: simulation dataset, the ray start and end points, and the line list;
    returns: the trident ray with the physical information we want to keep track of
    '''
    #out_tri_name = kwargs.get('out_tri_name', "temp.h5")
    ## first, figure out the fields
    field_list = ['metallicity', 'H_p0_number_density']

    ## now, to make sure we have the relevant info for the relevant lines

    ## now, make the ray
    triray = trident.make_simple_ray(ds, start_position=ray_start.copy(),
                              end_position=ray_end.copy(),
                             # data_filename=out_tri_name,
                              lines=line_list,
                              ftype='gas',
                              fields=['metallicity', 'H_p0_number_density'])
    return triray

def parse_for_ray_parameters(filename):
    hdu = fits.open(filename)
    if 'RAYSTART' not in hdu[0].header.keys():
        return False,False
    else:
        start,end = hdu[0].header['raystart'], hdu[0].header['rayend']
        if 'unitary' in start and end:
            if ',' in start and end:
                start = re.split(' unitary,| unitary',start)[:-1]
                end = re.split(' unitary,| unitary',end)[:-1]
            else:
                start,end = start.split()[1:-2],end.split()[1:-2]
        elif ',' in start and end:
            start,end = re.split(',',start),re.split(',',end)
        else:
            start,end = start.split(),end.split()

        start = [float(i) for i in start]
        end = [float(i) for i in end]

    return np.array(start),np.array(end)

def examine_column_density_thresholds(directory,enzo_filename):
    spectra_list = glob.glob(os.path.join(directory,'*.fits'))
    ions = ['SiII','SiIII','SiIV','CII','CIII','CIV','OVI','MgII','NeVIII']
    line_list = ['Si II','Si III','Si IV','C II','C III','C IV','O VI','Mg II','Ne VIII']
    field_list = ['Si_p1','Si_p2','Si_p3','C_p1','C_p2','C_p3','O_p5','Mg_p1','Ne_p7']
    t = Table(names=ions,dtype=[float,float,float,float,float,float,float,float,float])

    ds = yt.load(enzo_filename)
    trident.add_ion_fields(ds, ions=['H I','Si II', 'Si III', 'Si IV', 'C II','C III', 'C IV', 'O VI', 'Mg II','Ne VIII'])

    useable_spectra = []
    for spectrum in spectra_list:
        ## grab the ray data that I need
        print spectrum
        start,end = parse_for_ray_parameters(spectrum)
        if isinstance(start,type(False)):
            continue
        lray = trident.make_simple_ray(ds,start_position=start,end_position=end,lines=line_list)
        ldata = lray.all_data()
        thresh_here = []
        useable_spectra.append(spectrum)

        for ion in field_list:
            ion_to_use = ldata[ion+'_number_density']*ldata['dl']
            threshold,num_cells = get_fion_threshold(ion_to_use,0.8)
            thresh_here.append(threshold)
        t.add_row(thresh_here)

    spec_col = Column(data=useable_spectra, name='name')
    t.add_column(spec_col)
    cPickle.dump(t,open('threshold_table.cpkl','wb'),protocol=-1)
    return t
