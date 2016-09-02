import yt
import numpy as np
import matplotlib.pyplot as plt
import cPickle

## The correct file name goes here!
fn = "/Volumes/TOSHIBA EXT/r0058_l10/redshift0058"
ds = yt.load(fn)

## The halo center and virial radius
pos = [0.40328598,0.47176743,0.46131516]
rvir = 320.

## Here I define the thickness as the length over which you'll project. So here through projection 1 Mpc through the center
## Width is the length on one side of the final projection
## So for a projection along the x-direction, the thickness would be along x and the length of y,z would be width

thickness = yt.YTQuantity(1., 'Mpc')
width = yt.YTQuantity(1.,'Mpc')
axis = 'x' ## axis to project over (x,y,z)

##The thickness and the width need to be in code units to correctly set up the edges for the yt data region. Here I'm doing the conversion but I'm thinking you might just be able to do:

thickness,width = thickness.in_units('code_length'),width.in_units('code_length')

##I got it to print the conversion between cgs and code which is what's here
#thickness,width = thickness.in_units('cm'),width.in_units('cm')
#thickness = float(thickness.value)/1.10202798453e+26
#width = float(width.value)/1.10202798453e+26

#so first it's just deciding which axis you're projecting over
if axis=='x':
    axis = 0
elif axis=='y':
    axis = 1
elif axis=='z':
    axis = 2
else:
    print 'no correct axis'


## The region is defined by setting the center and then adding the correct thickness and width to the x,y,z axes
center = np.array(pos)
LE,RE = center.copy(),center.copy()

## Here it adds the thickness to the projection axis
LE[axis] -= thickness/2.0
RE[axis] += thickness/2.0

## Here it adds the width to the two remaining axes that will be shown in the projection
area_axes = [0,1,2]
i = area_axes.index(axis)
del area_axes[i]

LE[area_axes] -= width/2.0
RE[area_axes] += width/2.0


##Now that you have the center and the edges, you can create the yt data region
region = ds.region(center,LE,RE)

## you can interact with this like any other data object
## if you want to get the projected quantities out you ask yt for a projection of the region
obj2 = ds.proj('Temperature','x',data_source=region,weight_field='Density')   

##This is a different yt object. I generally prefer at this point to work with the data itself as a numpy array. You do that by asking for a "fixed resolution buffer" where you can set the resolution to whatever you want. The yt projections use [800,800] but you can choose the values to correspond to something physical given the width of the box that you set
resolution = [800,800]
frb_out = obj2.to_frb(width,resolution,center=center)

## Then to get the actual data
frb = np.array(frb_out['Temperature']) ## where here it's weighted by the density because that's how we created the projection

## it's then super easy to plot the data:
plt.imshow(np.log10(frb),vmin=4,vmax=7)

## or to save the data for later use
cPickle.dump(frb,open('fileout.cpkl','wb'),protocol=-1)

## so then it's trivial to reload
frb_again = cPickle.load(open('fileout.cpkl','rb'))





## Additional code for finding the conversion between cgs and code units
#reg = ds.unit_registry
#
#for un in reg.keys():
#    if un.startswith('code_'):
#        fmt_tup = (un, reg.lut[un][0], str(reg.lut[un][1]))
#        print ("Unit name:      {:<15}\nCGS conversion: {:<15}\nDimensions:     {:<15}\n".format(*fmt_tup)



)
