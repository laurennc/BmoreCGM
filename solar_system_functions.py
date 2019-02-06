import numpy as np
from astropy.table import Table
from astropy.table import Column
import json

def parse_MPC_comet_json(filename):
    fopen = open(filename)
    data = json.load(fopen)
    names = ['day_perihelion','epoch_day','epoch_month','epoch_year',
             'G','H','month_perihelion','node','orbit_type','peri',
             'perihelion_dist','provis_packed_desig','ref','year_perihelion','e','i']
    dtype = [int,int,int,int,float,float,int,float,str,float,float,str,str,int,float,float]
    t = Table(names=names,dtype=dtype)
    name_array = []
    count_bad = 0
    #print data[0].keys()
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
            all_vals = [day_perihelion,epoch_day,epoch_month,epoch_year,G,H,
                        month_perihelion,node,orbit_type,peri,perihelion_dist,
                        provis_packed_desig,ref,year_perihelion,e,i]
            name_array.append(name)
            t.add_row(all_vals)
        else:
            count_bad = count_bad + 1
    #print "Number of bad rows found: ",count_bad
    name_array = np.array(name_array)
    name_col = Column(name='name',data=name_array)
    t.add_column(name_col)
    return t

##### NOTE!!  #####
#Making executive decision to use objects with 30 keys since they are by far
#the most common object types. If we decide we want to make an even bigger
#sample, then I can work on generalizing the code.

## The one excpetion most worth considering is the designation fields
## I think it may be some of the biggest/best observed ones that don't have
## a "Principle desig" that's more numeric than their "Name"

def parse_MPC_json_table(filename):
    fopen = open(filename)
    data = json.load(fopen)

    names = ['rms','Peri','Tp','Epoch','Orbit_type','Ref','Node',#'Other_desigs','Name',
             'G','Perturbers_2','H','M','Num_opps','Perturbers','Orbital_period',
             'U','Num_obs','Arc_years','semimajor_axis','eccentricity','Last_obs',
             'inclination','Perihelion_dist','Number','n','Semilatus_rectum',
             'Hex_flags','Computer','Synodic_period','Aphelion_dist','Principal_desig']
    dtype = [float,float,float,float,str,str,float,float,str,float,float,int,str,float,
             float,int,str,float,float,str,float,float,int,float,float,str,str,float,float,str]

    t = Table(names=names,dtype=dtype)

    name_array = []
    count_bad = 0
    for row in range(len(data)):
        #print row
        if data[row].keys() == data[0].keys():
            rms              = data[row]['rms']
            peri             = data[row]['Peri']
            t_p              = data[row]['Tp']
            epoch            = data[row]['Epoch']
            orbit_type       = data[row]['Orbit_type']
            ref              = data[row]['Ref']
            node             = data[row]['Node']
            #other_desigs     = data[row]['Other_desigs']
            name             = data[row]['Name']
            gval             = data[row]['G']
            perturb2         = data[row]['Perturbers_2']
            hval             = data[row]['H']
            mval             = data[row]['M']
            num_opps         = data[row]['Num_opps']
            perturb           = data[row]['Perturbers']
            orbit_period     = data[row]['Orbital_period']
            uval             = data[row]['U']
            num_obs          = data[row]['Num_obs']
            arc_years        = data[row]['Arc_years']
            semimajor        = data[row]['a']
            eccentricity     = data[row]['e']
            last_obs         = data[row]['Last_obs']
            inclination      = data[row]['i']
            perihel_dist     = data[row]['Perihelion_dist']
            number           = data[row]['Number'][1:-1]
            nval             = data[row]['n']
            semilatus_rectum = data[row]['Semilatus_rectum']
            hex_flags        = data[row]['Hex_flags']
            computer         = data[row]['Computer']
            synodic_period   = data[row]['Synodic_period']
            aphel_dist       = data[row]['Aphelion_dist']
            principal_desig  = data[row]['Principal_desig']

            all_vals = [rms,peri,t_p,epoch,orbit_type,ref,node,#other_desigs,name
                        gval,perturb2,hval,mval,num_opps,perturb,orbit_period,
                        uval,num_obs,arc_years,semimajor,eccentricity,last_obs,
                        inclination,perihel_dist,number,nval,semilatus_rectum,
                        hex_flags,computer,synodic_period,aphel_dist,principal_desig]
            #print len(all_vals)
            name_array.append(name)
            t.add_row(all_vals)
        else:
            count_bad = count_bad + 1

    #print "Number of bad rows found: ",count_bad
    name_array = np.array(name_array)
    name_col = Column(name='name',data=name_array)
    t.add_column(name_col)
    return t

def parse_MPC_ascii_table(filename):
    fopen = open(filename)
    names = ['abs_mag','slope','epoch','anomaly','perihelion',
			 'long_ascend_node','inclination','eccentricity','mean_daily_motion',
			 'semimajor_axis','uncertainty_param','reference','num_observ',
			 'num_opposition','year_first_obs','year_last_obs','arc_length',
			 'days','rms_residual']

    dtype = [float,float,str,float,float,float,float,float,float,
	         float,str,str,int,int,int,int,str,str,float]

    t = Table(names=names, dtype=dtype)
    name_array = []

    content = fopen.readlines()
    i = 0
    for line in content:
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

        data_vals = [abs_mag,slope,epoch,anomaly,perihelion,
		             long_ascend_node,inclination,eccentricity,mean_daily_motion,
		             semimajor_axis,uncertainty_param,reference,num_observ,
		             num_opposition,year_first_obs,year_last_obs,arc_length,
		             days,rms_residual]

        name_array.append(name)
        t.add_row(data_vals)

        i = i + 1
        
    name_array = np.array(name_array)
    name_col = Column(name='name',data=name_array)
    t.add_column(name_col)
    return t

## This is using data from this paper: https://arxiv.org/abs/1601.02087
def parse_asteroid_color_table(filename):
    f = open(filename)
    names = ['number','name1','name2','class','r','r_err','i','i_err',
             'z','z_err','semimajor_axis','eccentricity','inclination','H']
    dtypes = [int,object,object,object,float,float,
              float,float,float,float,float,float,float,float]
    t = Table(names=names, dtype=dtypes)
    for line in f.readlines():
        vals = line.split()
        if vals[0] == '#':
            continue
        t.add_row(vals)

    return t

def calculate_semimajor_axis(perihelion,eccentricity):
    aphelion = perihelion*(1.+eccentricity)/(1.-eccentricity)
    return 0.5*(aphelion+perihelion)
