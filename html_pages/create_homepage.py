import numpy as np
import yt
from astropy.table import Table
import os
from astropy.units.quantity import Quantity

def add_table_entry():
	return

def create_page():
	t = Table.read('test_list.dat',format='ascii',names=('testname','boxsize','unit','amr','dens_cloud','density_ratio','machnum','cloud_type'))

	run_catalog = Table(names=('Run Name','Box Size','AMR On?','Cloud Density [g/cm^3]','Density Ratio','Mach Number','Cloud Type' ),dtype=('S200','S200','S    200','float64','float64','float64','int64'))

	for i in range(len(t['testname'])):	
		repo_link = '<a href=\\' + t['testname'][i]+'\\run_display.html > '+ t['testname'][i] + ' </a>'	
		box_unit = str(Quantity(t['boxsize'][i],t['unit'][i]))
		if t['amr'][i] == 1:
			amr_out = 'Yes'
		else: 
			amr_out = 'No'	

		run_catalog.add_row([repo_link,box_unit,amr_out,t['dens_cloud'][i],t['density_ratio'][i],t['machnum'][i],t['cloud_type'][i]  ])


	os.system('cp main_page_head.html main_page_display.html')
	run_catalog.write("output.temp",format="jsviewer")
	os.system('sed "s/&lt;/</g" output.temp | sed "s/&gt;/>/g" >> main_page_display.html')

	os.system('cat main_page_tail.html >> main_page_display.html')
	os.system('rm output.temp')	

	return 'main page built!'


create_page()




