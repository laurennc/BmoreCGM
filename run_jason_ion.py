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
import seaborn as sns

sns.set_style("white", {'axes.grid' : False})

RDname = '/Users/dalek/data/Molly/natural/nref11/RD0020/RD0020'
trackfile = '/Users/dalek/data/Molly/natural/nref11/halo_track'
ion_list = ['O VI','C IV']

ad,rb,rb_width = prep_dataset(RDname,trackfile,ion_list,region='trackbox')

df = prep_dataframe(ad,rb,rb_width,'O_p5_ion_fraction','C_p3_ion_fraction')

img = render_image(df,'density','temperature','o6frac',(-30.,-21.),(2,8.5),'ovi_frac_test')
wrap_axes('ovi_frac_test','density','temperature',[(-30.,-21.),(2,8.5)])
