import numpy as np

def prepare_shiftc_file(start,stop,ref=True):
    DDs = np.arange(start,stop,1)
    fout = open('untar.exe','w')
    for DD in DDs:
        if DD < 1000:
            DDstr = 'DD0'+str(DD)
        else:
            DDstr = 'DD'+str(DD)
        if ref == True:
            dirout = '/lou/s2m/lcorlies/nref11n/nref10_refine200kpc_z2to1/'
        else:
            dirout = '/lou/s2m/lcorlies/nref11n/natural'
        lineout = 'shiftc --extract-tar '+DDstr+'.tar '+dirout+' \n'
        fout.write(lineout)
    fout.close()
    return

prepare_shiftc_file(400,500,ref=False)

#shiftc --extract-tar DD0354.tar
