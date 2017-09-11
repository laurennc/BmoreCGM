import yt
import numpy as np
import matplotlib.pyplot as plt
import trident

filenames = ['/Users/dalek/data/Jason/nref10/RD0042/RD0042','/Users/dalek/data/Jason/nref11_mustrefine/RD0042/RD0042']
outputnames = ['nref10','nref11_mustrefine']
ions = ['O VI','C IV','C III','Si III']
width = [(55.,'kpc'),(105.,'kpc')]
halo_center = [0.48984336853027843, 0.47133064270019531, 0.50956535339355846]
fields = ['Density','H_number_density',('gas','O_p5_number_density'),('gas','C_p2_number_density')]
axis = 'z'

make_projection = True

## for the early must refine regions
mrl = [0.48966962,0.47133064,0.5093916]
mrr = [0.49001712,0.47202564,0.5097391]
xc = (mrr[0]-mrl[0])/2.0 + mrl[0]
yc = (mrr[1]-mrl[1])/2.0 + mrl[1]
zc = (mrr[2]-mrl[2])/2.0 + mrl[2]

box_center = [xc,yc,zc]

keywords = {'Density':('Blues_r','g/cm**2',1e-4,1e-1)}
	    #'O_p5_number_density':(}


if make_projection:
	for i in range(len(filenames)):
		ds = yt.load(filenames[i])
		trident.add_ion_fields(ds, ions=ions)
		for field in fields:
			proj = yt.ProjectionPlot(ds,axis,field,center=box_center,width=width)
			if field in keywords:
			        cmap,unit,zlow,zhigh = keywords[field]  
				proj.set_cmap(field,cmap)
				proj.set_zlim(field,zlow,zhigh)
				proj.set_unit(field,unit)
			proj.save(outputnames[i])

