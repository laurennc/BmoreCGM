import numpy as np
from astropy.table import Table
import json

def parse_MPC_comet_json(filename):
    fopen = open(filename)
    data = json.load(fopen)
    names = ['name','day_perihelion','epoch_day','epoch_month','epoch_year',
             'G','H','month_perihelion','node','orbit_type','peri',
             'perihelion_dist','provis_packed_desig','ref','year_perihelion','e','i']
    dtype = [str,int,int,int,int,float,float,int,float,str,float,float,str,str,int,float,float]
    t = Table(names=names,dtype=dtype)
    count_bad = 0
    print data[0].keys()
    for row in range(len(data)):
        if data[row].keys() == data[0].keys():
            name                = data[row]['Designation_and_name']
            day_perihelion      = data[row]['Day_of_perihelion']
            epoch_day           = data[row]['Epoch_day']
            epoch_month         = data[row]['Epoch_month']
            epoch_year          = data[row]['Epoch_year']
            G                   = data[row]['G']
            H                   = data[row]['H']
            month_perihelion    = data[row]['Month_of_perihelion']
            node                = data[row]['Node']
            orbit_type          = data[row]['Orbit_type']
            peri                = data[row]['Peri']
            perihelion_dist     = data[row]['Perihelion_dist']
            provis_packed_desig = data[row]['Provisional_packed_desig']
            ref                 = data[row]['Ref']
            year_perihelion     = data[row]['Year_of_perihelion']
            e                   = data[row]['e']
            i                   = data[row]['i']
            all_vals = [name,day_perihelion,epoch_day,epoch_month,epoch_year,G,H,
                        month_perihelion,node,orbit_type,peri,perihelion_dist,
                        provis_packed_desig,ref,year_perihelion,e,i]
            t.add_row(all_vals)
        else:
            count_bad = count_bad + 1
    print "Number of bad rows found: ",count_bad
    return t

def parse_MPC_table(filename):
    fopen = open(filename)
    names = ['name','abs_mag','slope','epoch','anomaly','perihelion',
			 'long_ascend_node','inclination','eccentricity','mean_daily_motion',
			 'semimajor_axis','uncertainty_param','reference','num_observ',
			 'num_opposition','year_first_obs','year_last_obs','arc_length',
			 'days','rms_residual']

    dtype = [str,float,float,str,float,float,float,float,float,float,
	         float,str,str,int,int,int,int,str,str,float]

    t = Table(names=names, dtype=dtype)

    content = fopen.readlines()
    i = 0
    for line in content:
		#print i
		name              =  line[0:8]
		abs_mag           =  line[8:14]
		if abs_mag.isspace():
 			abs_mag = -9999
		abs_mag           =  float(abs_mag)
		slope             =  line[14:20]
		if slope.isspace():
			slope = -9999
		slope             =  float(slope)
		epoch             =  line[20:26]
		anomaly           =  float(line[26:36])
		perihelion        =  float(line[37:47])
		long_ascend_node  =  float(line[48:58])
		inclination       =  float(line[59:69])
		eccentricity      =  float(line[70:80])
		mean_daily_motion =  float(line[80:92])
		semimajor_axis    =  float(line[92:104])
		uncertainty_param =  line[105]
		reference         =  line[107:117]
		num_observ        =  int(line[117:123])
		num_opposition    =  int(line[123:127])
		if num_opposition > 1:
			year_first_obs    =  int(line[127:131])
			year_last_obs     =  int(line[132:137])
			arc_length = ''
			days = ''
		else:
			arc_length        = line[127:142]
			days              = line[132:136]
			year_first_obs = 0
			year_last_obs = 0

		rms_residual      =  float(line[137:142])

		data_vals = [name,abs_mag,slope,epoch,anomaly,perihelion,
		             long_ascend_node,inclination,eccentricity,mean_daily_motion,
		             semimajor_axis,uncertainty_param,reference,num_observ,
	       	             num_opposition,year_first_obs,year_last_obs,arc_length,
			     days,rms_residual]

		t.add_row(data_vals)

		i = i + 1
    return t