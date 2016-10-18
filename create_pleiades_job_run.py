import os
import numpy as np

def update_enzo_params(f_in,dens_c,mach,dens_amb,rad_c,temp_c,cloud_type):
	fin = open(f_in)
	lines = fin.readlines()
	fout = open('temp.out','w')
	for line in lines:
		line = line.replace("CloudCrusherCloudDensity     = 1000.0","CloudCrusherCloudDensity     = "+str(dens_c))
		line = line.replace("CloudCrusherShockMachNumber = 8.8","CloudCrusherShockMachNumber = "+str(mach))  
		line = line.replace("DensityUnits = 1.67e-24","DensityUnits = "+str(dens_amb))  
		line = line.replace("CloudCrusherCloudRadius      = 0.0625","CloudCrusherCloudRadius      = "+str(rad_c))
		line = line.replace("CloudCrusherCloudTemperature = 1000.0","CloudCrusherCloudTemperature = "+str(temp_c))  
		line = line.replace("CloudCrusherCloudType        = 2","CloudCrusherCloudType        = "+str(cloud_type))
		fout.write(line+'\n')
	os.system('mv temp.out '+fin)
	return

def update_runscript(f_in,dir_name):
	fin = open(f_in)
	lines = fin.readlines()
	fout = open('temp.out','w')
	for line in lines:
		line = line.replace("#PBS -N NAME_HERE","#PBS -N "+dir_name)
		fout.write(line+'\n')
	return


## need to be able to type python ..py a b c d e f
## and then it should execute

dens_c = [100.,1000.,10000.]
mach = [1.,5.,10.,15.]
dens_amb = [1e-28,1e-27,1e-26,1e-25]
rad_c = [0.01]
temp_c = [1000.,10000.]
cloud_type = [0,1,2]

#run_dc1e2_m1_da1e29_rc?_tc_1e3_ct0

a,b,c,d,e,f = 0,0,0,0,0,0
fout = open('test.out','w')

HOME = "/Volumes/sonic/CloudCrusher/TestPleiadesCode"
PLOTS = "/Volumes/sonic/CloudCrusher/TestPleiadesCode/QuickLook"

#### NEED TO ADD IN UPDATED TO A SCRIPT THAT WILL SUBMIT MULTIPLE JOBS
frun = open('jobs_to_run.exe','w')
#### NEED TO ADD TO A FILE THAT'S KEEPING TRACK OF THE RUN PARAMETERS
flist = open('run_list.dat','w')
### NEED TO ADD SCRIPT TO AUTOMATE THE FILES THAT I WANT AND WHERE THEY SHOULD GO
### NEED TO ADD COMMAND TO CREATE THE QUICK LOOK PAGES
#### 
### I SHOULD STORE ALL OF THESE THINGS IN A NEW DIRECTORY SO I CAN EASILY TRANSFER THAT HOME INSTEAD OF THE DATA
### 

while a < len(dens_c):
	while b < len(mach):
		while c < len(dens_amb):
			while d < len(rad_c):
				while e < len(temp_c):
					while f < len(cloud_type):
						fout = "run_cd"+str(int(dens_c[a]))+"_m"+str(int(mach[b]))+"_da"+str(dens_amb[c])+"_rc"+str(rad_c[d])+"_tc"+str(int(temp_c[e]))+"_ct"+str(int(cloud_type[f]))
						os.system('mkdir ' + fout)
						os.system('cp CloudCrusher.enzo '+fout+'/')
						os.system('cp RunScript.sh '+fout+'/')
						os.system('cp enzo.exe '+fout+'/')
						os.system('cp cool_rates.in '+fout+'/')
						os.system('cp output_quick_look.py '+fout+'/')
						os.system('mkdir '+PLOTS+'/'+fout)
						flist.write(fout+' '+str(dens_c[a])+' '+str(mach[b])+' '+str(dens_amb[c])+' '+str(temp_c[e])+' '+str(cloud_type[f]    )+'\n')                                         
						os.chdir(fout)
						update_enzo_params('CloudCrusher.enzo',dens_c[a],mach[b],dens_amb[c],rad_c[d],temp_c[e],cloud_type[f])
						update_runscript('RunScript.sh',fout)
						frun.write('cd '+'fout'+'/ \n')
						frun.write('qsub RunScript.sh \n')
						ps.chdir(HOME)
						f = f + 1
					f = 0
					e = e + 1
				e = 0
				d = d + 1
			d = 0
			c = c + 1
		c = 0
		b = b + 1
	b = 0
	a = a + 1

frun.close()
flist.close()


