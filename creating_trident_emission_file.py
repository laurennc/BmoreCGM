import h5py
import numpy as np
import matplotlib.pyplot as plt

nhden = 29
nred  = 25
ntemp = 161
nline = 15

path_to_Cloudy_data = '/Users/dalek/data/cloudy_data/hm_2011/'

redshift_bins = ['0.0000e+00','1.2202e-01','2.5893e-01','4.1254e-01','5.8489e-01',
    '7.7828e-01','9.9526e-01','1.2387e+00','1.5119e+00','1.8184e+00',
    '2.1623e+00','2.5481e+00','2.9811e+00','3.4668e+00','4.0119e+00',
    '4.6234e+00','5.3096e+00','6.0795e+00','6.9433e+00','7.9125e+00',
    '9.0000e+00','1.0220e+01','1.1589e+01','1.3125e+01','1.4849e+01']

zbins = np.array(redshift_bins).astype(float)
tbins = np.linspace(1,9,161)
hbins = np.linspace(-10,4,29)

roman_numerals = {'1':'I','2':'II','3':'III','4':'IV','5':'V','6':'VI','7':'VII','8':'VIII'}

#lines = ['HI_1216','HI_6563','CIV_1548','CIV_1551','OVI_1032','OVI_1038',
#        'CIII_977','CIII_1907','CIII_1910','SiII_1814','SiIII_1207',
#        'SiIV_1394','SiIV_1403','MgII_2796','MgII_2803']

## a more flexible way to do this would be to read the lines from the header
## of the .dat files and transform them to a user-friendly name
## break by underscore, throw out all underscores, add surrounding words

def set_up_line_structure(path_to_Cloudy_data):
  fin = open(path_to_Cloudy_data+'z_'+redshift_bins[0]+'/hm_2011_run1.dat')
  flines = fin.readlines()
  words = flines[11].split()
  lines = []
  i = 1
  while i < len(words):
    pieces = words[i].split('_')
    pieces = list(filter(None,pieces))
    lout = pieces[0]+roman_numerals[pieces[1]]+'_'+pieces[2]
    lines = np.append(lines,lout)
    i = i + 1
  return list(lines)

def create_dataset_list(lines):
  datasets = {}
  for line in lines:
    datasets[line] = np.zeros((nhden,nred,ntemp))
  return datasets

def remove_zeros(datasets):
    keys = datasets.keys()
    for key in keys:
        ## this value seems to be the double float floor that Cloudy used for the calculations
        datasets[key][datasets[key] == 0.] = -313.65269999999998
    return datasets

def add_attrs(datasets):
    for key in datasets.keys():
        datasets[key].attrs.create('Parameter1',hbins)
        datasets[key].attrs.create('Parameter2',zbins)
        datasets[key].attrs.create('Temperature',tbins)
    return datasets

lines = set_up_line_structure(path_to_Cloudy_data)
datasets = create_dataset_list(lines)
nlines = len(lines)

for i in range(nred):
  for j in range(nhden):
    fin  = path_to_Cloudy_data+'z_'+redshift_bins[i]+'/hm_2011_run'+str(j+1)+'.dat'
    try:
        data = np.genfromtxt(fin)
        for k in range(nlines):
            datasets[lines[k]][j,i,:] = data[:,k]
            #datasets[lines[k]][j,i,:] = -313.65269999999998
    except ValueError:
        holder = 0
        #datasets[lines[k]][j,i,:] = -313.65269999999998
        #print fin

## I need to remove zeros from the datasets
datasets = remove_zeros(datasets)
## I need to finally build the hdf5 structure

hfT = h5py.File('emissivity_tables.h5', 'w')
for line in lines:
    hfT.create_dataset(line,data=datasets[line])

add_attrs(hfT)

hfT.close()
