import numpy as np
import matplotlib.pyplot as plt


def make_thin_projections(ds,field,weight,axis,center,width,thickness,unit,resolution,fileout):
	thickness = yt.YTQuantity(thickness,unit)
	width = yt.YTQuantity(width,unit)
	thickness,width = thickness.in_units('cm'),width.in_units('cm')
	thickness = float(thickness.value)/1.10202798453e+26
	width = float(width.value)/1.10202798453e+26

	if axis=='x':
		axis = 0
	elif axis=='y':
		axis = 1
	elif axis=='z':
		axis = 2
	else:
		print 'no correct axis'

	center = np.array(center)
	LE,RE = center.copy(),center.copy()
	
	LE[axis] -= thickness/2.0
	RE[axis] += thickness/2.0

	area_axes = [0,1,2]
	i = area_axes.index(axis)
	del area_axes[i]

	LE[area_axes] -= width/2.0
	RE[area_axes] += width/2.0

	region = ds.region(center,LE,RE)
	obj = ds.proj(field,axis,data_source=region,weight_field=weight)

	frb = obj.to_frb(width,resolution,center=center)

	cPickle.dump(frb[field],open(fileout,'wb'),protocol=-1)

	return frb[field]


	return
