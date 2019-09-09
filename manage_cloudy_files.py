import numpy as np
import os, sys
import shutil

def update_cloudy_file(filein,temperature):
    f = open(filein)
    lines = f.readlines()
    f_out_name = filein + '.TEMP'
    f_out = open(f_out_name,'w')
    
    for line in lines:
        words = line.split()
        #print words, len(words)
        if len(words) == 1:
           new_line = line
        else:
	    if words[1] == 'temperature':
	       new_line = 'constant temperature '+str(temperature)+' K linear \n'
            else:
               new_line = line
        f_out.write(new_line)
    
    f_out.close()
    f.close()
    shutil.move(f_out_name,filein)
    return

def build_cloudy_lines_file(filein,fileout,temp_val):
    f_in = open(filein)
    f_out = open(fileout,'a+')

    lines = f_in.readlines() 
    vals = lines[1].split()
    line_out2 = "   ".join(vals[1:])
    line_out = str(temp_val)+"   "+line_out2 +' \n'
    f_out.write(line_out)
    f_in.close()
    f_out.close()
    return
 

def run_cloudy_loop(filein):
    temp_vals = np.linspace(3.,8.,51)
    fileout =  'hden-3.dat'    

    for temp in temp_vals:
        temp_update = 10**temp
        update_cloudy_file(filein,temp_update)
        print 'Temperature = '+str(temp)
        os.system(r"/Users/lcorlies/codes/c17.01/source/cloudy.exe -r hden-3")
        build_cloudy_lines_file('hden-3.lines',fileout,temp)
        
    return

run_cloudy_loop('hden-3.in')



