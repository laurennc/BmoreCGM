import numpy as np
import cPickle
import matplotlib.pyplot as plt
from yt.mods import *

def diskVectors(pf,center):
	center = np.array(center)
	sphere = pf.h.sphere(center,(5,'kpc'))
	angular_momentum = sphere.quantities['AngularMomentumVector']()
	x = np.cross(angular_momentum,[1.0,0.0,0.0]
	x /= np.linalg.norm(x))
	return center,angular_momentum,x

def getFrameAxes(pf,center,view='edge-on'):
	center,angular_momentum,x = diskVectors(pf,center)
	axis = angular_momentum
	tmp_north = np.cross(angular_momentum,[1.0,0.0,0.0])
	north_vector = tmp_north / np.linalg.norm(tmp_north)
	if view == 'edge-on':
		north_vector = angular_momentum
		axis = x
	return axis,north_vector


def getFrame(pf,fType,view,var,fileout,center,radius):
	axis,north_vector = getFrameAxes(pf,center,view=view)
	if fType == 'Projection'
		frame = lp.OffAxisProjectionPlot(pf,axis,var)
	elif fType == 'Slice':
		frame = lp.OffAxisSlicePlot(pf,axis,var)
	else:
		frame = None

	frb = frame[var].image.get_array()
	return frb


##CHECK IF THE OFF AXIS PROJECTIONS CAN HAVE A THICKNESS AND WIDTH
## ... ALSO I don't want either of these... I want a full cube...
##check out the code that I used for Mark and Lluis


	

	return

