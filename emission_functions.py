import numpy as np
from scipy import interpolate
import yt
import yt.units as u
#import builtins

patt = "/Users/dalek/data/cloudy_data/emissivity_tables/bertone_z2/bertone_z2_run%i.dat"

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


def scale_by_metallicity(values,assumed_Z,wanted_Z):
        wanted_ratio = (10.**(wanted_Z))/(10.**(assumed_Z))
        return values*wanted_ratio

def parse_cloudy_lines_file(filein,fileout):
    f = open(filein,'r')
    fh = open(filein.split('.lines')[0]+'.dat','r')
    fo = open(fileout,'w')

    foL = fh.readlines()
    hden = foL[5]

    linenum = 0
    header = True
    fo.write('# \n# Cooling Map File \n# Loop values: \n# \n')
    fo.write(hden)
    fo.write('\n# Data Columns: \n # log10 Te [K] \n')
    fo.write('# log10 Emissivities / n_H^2 [erg s^-1 cm^3] \n# \n')
    header_words = '#Te'

    for line in f:
        print linenum
        if linenum == 0:
            words = line.split()
            temp = str(np.log10(float(words[-1][:-1])))
            linenum = linenum + 1
            continue
        if linenum == 1:
            if header == True:
                print line
                words = line.split('depth')[1]
                header_words = header_words + words+' \n'
                fo.write(header_words)
                header = False
            linenum = linenum + 1
            continue
        if linenum == 2:
            nums = line.split()[1:]
            #print nums
            #vals = vals + str(np.log10(float(x))) for x in nums
            line_out = temp +'\t'+ '\t'.join(nums) + '\n'
            #print line_out
            fo.write(line_out)
            linenum = 0
            continue
    fo.write('\n')
    fo.close()
    return

def make_Cloudy_table(table_index):
        hden_n_bins, hden_min, hden_max = 23, -9, 2
        T_n_bins, T_min, T_max = 71, 2, 8
        #patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/grid_galquas/emis/z0/g1q1/g1q1_run%i.dat"
        #patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/test_hden_step/g1q1_run%i.dat"

        hden=np.linspace(hden_min,hden_max,hden_n_bins)
        T=np.linspace(T_min,T_max, T_n_bins)
        table = np.zeros((hden_n_bins,T_n_bins))
        for i in range(hden_n_bins):
                #print i,patt%(i+1)
                table[i,:]=[float(l.split()[table_index]) for l in open(patt%(i+1)) if l[0] != "#"]
        return hden,T,table



emission_units = 's**-1 * cm**-3 * steradian**-1'
ytEmU = u.s**-1 * u.cm**-3 * u.steradian**-1

hden_pts,T_pts,table_HA = make_Cloudy_table(2)
hden_pts,T_pts = np.meshgrid(hden_pts,T_pts)
pts = np.array((hden_pts.ravel(),T_pts.ravel())).T

sr_HA = table_HA.T.ravel()
bl_HA = interpolate.LinearNDInterpolator(pts,sr_HA)

def _Emission_HAlpha(field,data):
    H_N = np.log10(np.array(data['H_nuclei_density']))
    Temperature = np.log10(np.array(data['Temperature']))
    dia1 = bl_HA(H_N,Temperature)
    idx = np.isnan(dia1)
    dia1[idx] = -200.
    emission_line=(10.**dia1)*((10.**H_N)**2.0)
    emission_line = emission_line/(4.*np.pi*3.03e-12)
    return emission_line*ytEmU

yt.add_field(('gas','Emission_HAlpha'),units=emission_units,function=_Emission_HAlpha)

hden_pts,T_pts,table_LA = make_Cloudy_table(1)
sr_LA = table_LA.T.ravel()
bl_LA = interpolate.LinearNDInterpolator(pts,sr_LA)

def _Emission_LyAlpha(field,data):
    H_N=np.log10(np.array(data["H_nuclei_density"]))
    Temperature=np.log10(np.array(data["Temperature"]))
    dia1 = bl_LA(H_N,Temperature)
    idx = np.isnan(dia1)
    dia1[idx] = -200.
    emission_line = (10**dia1)*((10.0**H_N)**2.0)
    emission_line = emission_line/(4.*np.pi*1.63e-11)
    return emissoin_line*ytEmU

yt.add_field(('gas','Emission_LyAlpha'),units=emission_units,function=_Emission_LyAlpha)

hden1, T1, table_SiIV_1 = make_Cloudy_table(13)
hden1, T1, table_SiIV_2 = make_Cloudy_table(14)
sr_SiIV_1 = table_SiIV_1.T.ravel()
sr_SiIV_2 = table_SiIV_2.T.ravel()
bl_SiIV_1 = interpolate.LinearNDInterpolator(pts,sr_SiIV_1)
bl_SiIV_2 = interpolate.LinearNDInterpolator(pts,sr_SiIV_2)

def _Emission_SiIV(field,data):
        H_N=np.log10(np.array(data["H_nuclei_density"]))
        Temperature=np.log10(np.array(data["Temperature"]))
        dia1 = bl_SiIV_1(H_N,Temperature)
        dia2 = bl_SiIV_2(H_N,Temperature)
        idx = np.isnan(dia1)
        dia1[idx] = -200.
        dia2[idx] = -200.
        emission_line=((10.0**dia1)+(10**dia2))*((10.0**H_N)**2.0)
        emission_line = emission_line/(4.*np.pi*1.42e-11)
        emission_line = scale_by_metallicity(emission_line,0.0,np.log10(np.array(data['metallicity'])))
	return emission_line*ytEmU

yt.add_field(("gas","Emission_SiIV"),units=emission_units,function=_Emission_SiIV)

hden1,T1,table_CIII_977 = make_Cloudy_table(7)
sr_CIII_977 = table_CIII_977.T.ravel()
bl_CIII_977 = interpolate.LinearNDInterpolator(pts,sr_CIII_977)

def _Emission_CIII_977(field,data):
        H_N=np.log10(np.array(data["H_nuclei_density"]))
        Temperature=np.log10(np.array(data["Temperature"]))
        dia1 = bl_CIII_977(H_N,Temperature)
        idx = np.isnan(dia1)
        dia1[idx] = -200.
        emission_line=(10.0**dia1)*((10.0**H_N)**2.0)
        emission_line = emission_line/(4.*np.pi*2.03e-11)
        emission_line = scale_by_metallicity(emission_line,0.0,np.log10(np.array(data['metallicity'])))
	return emission_line*ytEmU

yt.add_field(("gas","Emission_CIII_977"),units=emission_units,function=_Emission_CIII_977)

hden1, T1, table_CIV_1 = make_Cloudy_table(3)
hden1, T1, table_CIV_2 = make_Cloudy_table(4)
sr_CIV_1 = table_CIV_1.T.ravel()
sr_CIV_2 = table_CIV_2.T.ravel()
bl_CIV_1 = interpolate.LinearNDInterpolator(pts,sr_CIV_1)
bl_CIV_2 = interpolate.LinearNDInterpolator(pts,sr_CIV_2)

def _Emission_CIV(field,data):
        H_N=np.log10(np.array(data["H_nuclei_density"]))
        Temperature=np.log10(np.array(data["Temperature"]))
        dia1 = bl_CIV_1(H_N,Temperature)
        dia2 = bl_CIV_2(H_N,Temperature)
        idx = np.isnan(dia1)
        dia1[idx] = -200.
        dia2[idx] = -200.
        emission_line=((10.0**dia1)+(10**dia2))*((10.0**H_N)**2.0)
        emission_line = emission_line/(4.*np.pi*1.28e-11)
        emission_line = scale_by_metallicity(emission_line,0.0,np.log10(np.array(data['metallicity'])))
	return emission_line*ytEmU

yt.add_field(("gas","Emission_CIV"),units=emission_units,function=_Emission_CIV)

hden1, T1, table_OVI_1 = make_Cloudy_table(5)
hden1, T1, table_OVI_2 = make_Cloudy_table(6)
sr_OVI_1 = table_OVI_1.T.ravel()
sr_OVI_2 = table_OVI_2.T.ravel()
bl_OVI_1 = interpolate.LinearNDInterpolator(pts,sr_OVI_1)
bl_OVI_2 = interpolate.LinearNDInterpolator(pts,sr_OVI_2)

def _Emission_OVI(field,data):
        H_N=np.log10(np.array(data["H_nuclei_density"]))
        Temperature=np.log10(np.array(data["Temperature"]))
        dia1 = bl_OVI_1(H_N,Temperature)
        dia2 = bl_OVI_2(H_N,Temperature)
        idx = np.isnan(dia1)
        dia1[idx] = -200.
        dia2[idx] = -200.
        emission_line=((10.0**dia1)+(10**dia2))*((10.0**H_N)**2.0)
        emission_line = emission_line/(4.*np.pi*1.92e-11)
        emission_line = scale_by_metallicity(emission_line,0.0,np.log10(np.array(data['metallicity'])))
	return emission_line*ytEmU

yt.add_field(("gas","Emission_OVI"),units=emission_units,function=_Emission_OVI)
