#! /usr/bin/env python

'''
AUTHOR: Lauren Corlies & Molly Peeples
DATE: 12/18/16
NAME: output_spectra_html.py
DESCRIPTION: makes spectra plots and html files for blobs
'''

import yt
import trident

from astropy.table import Table, vstack

import matplotlib.pyplot as plt
import numpy as np

import os
import sys
import glob

import argparse

#-----------------------------------------------------------------------------------------------------

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description = \
                                     "makes some plots")

    parser.add_argument('--clobber', dest='clobber', action='store_true')
    parser.add_argument('--no-clobber', dest='clobber', action='store_false', help="default is no clobbering")

    parser.set_defaults(clobber=False)

    args = parser.parse_args()
    return args



#-----------------------------------------------------------------------------------------------------
def make_trident_plots(ds):
    line_list = ['H','C','N','O','Mg','S','Si','Ne']

    # Go through the cloud perpendicular and parallel to the wind
    ray_starts = [[1,0,1], [0,1,1]]
    ray_ends   = [[1,2,1], [4,1,1]]
    i = 0   
 
    while i < len(ray_starts):
        ray_start,ray_end = ray_starts[i], ray_ends[i]
        ray = trident.make_simple_ray(ds,start_position=ray_start,end_position=ray_end,
                                      data_filename="ray.h5",lines=line_list,ftype='gas')

        sg = trident.SpectrumGenerator('COS-G130M')
        sg.make_spectrum(ray,lines=line_list)
        sg.plot_spectrum('Spec'+str(i)+'_'+ds.basename+'.png')
        
        os.remove('ray.h5')
        i = i + 1




#-----------------------------------------------------------------------------------------------------

def make_quicklook_page(ds):
    # put some stuff on an html page
    outfilename = ds.basename + "_quicklook.html"
    outfile = open(outfilename, "w")
    info = """<html>
    <head><title>"""+ds.basename+"""</title>
    <h1 style="text-align:center;font-size:350%">"""+ds.basename+"""</h1></head>
    <body>
    <p style="text-align:center;font-size=250%">t = """+str(ds.current_time)+"""<br>
    <hr />
    """
    outfile.write(info)

    ## now add the figures
    print "figure names should be dynamic or passed, but not currently"
    spec1_string = 'Spec1_'+ds.basename+'.png'
    spec2_string = 'Spec2_'+ds.basename+'.png'

    addfig = r"""<img src='"""+spec1_string+r"""' style="width:35%">"""
    outfile.write(addfig)
    addfig = r"""<img src='"""+spec2_string+r"""' style="width:35%">"""
    outfile.write(addfig)


    outfile.write("</body></html>")
    outfile.close()




#-----------------------------------------------------------------------------------------------------

def look_at_output(clobber):
    ## is there an output defined? if not, do all!
    print " ---> WARNING OUTPUT NAMES ARE HARDCODED, ASSUMES DD#### FORMAT <---- "
    DATA_DIR = "."
    PLOTS_DIR = "./plots"
    dataset_list = glob.glob(os.path.join(DATA_DIR, 'DD????/DD????'))

    print " ---> I'm not checking the plots have been made I'm just gonna make them <---- "
    run_catalog = Table(names=("rootname","time","Perpendicular","Parallel"),dtype=("S200","f8","S200","S200"))
    for output in dataset_list:
        ds = yt.load(output)
        spec1_string = 'Spec1_'+ds.basename+'.png'
        spec2_string = 'Spec2_'+ds.basename+'.png'
        if clobber or not (os.path.exists(density_string) and os.path.exists(density_string)):
            make_trident_plots(ds)
        if clobber:
            make_quicklook_page(ds)
        ds_link = '<a href=\"'+ ds.basename + '_quicklook.html\">' + ds.basename + '</a>'
        spec1_html_string = '<a href=\"'+spec1_string+'\"><img height=\"200px\" src=\"'+density_string+'\"></a>'
        spec2_html_string = '<a href=\"'+spec2_string+'\"><img height=\"200px\" src=\"'+temperature_string+'\"></a>'
        run_catalog.add_row([ds_link, ds.current_time, spec1_html_string, spec2_html_string])

    ## make a table with the thumbnails and basic information
    run_catalog.write("run_catalog.temp",format="jsviewer")
    os.system('sed "s/&lt;/</g" run_catalog.temp | sed "s/&gt;/>/g" > run_catalog.html')
    os.system('rm run_catalog.temp')

    ## now generate QL pages


    return "yay simulations!"

#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    ## which outputs?
    args = parse_args()
    clobber = args.clobber

    message = look_at_output(clobber)
    sys.exit(message)

