from astropy.table import Table

track = Table.read('nref11_track_box', format='ascii')
#track.sort('col1')

f = open('nref11_cool_small_region','w')

for i in range(len(track['col1'])):
  x_left =  track['col2'][i]
  y_left =  track['col3'][i]
  z_left =  track['col4'][i]
  x_right = track['col5'][i]
  y_right = track['col6'][i]
  z_right = track['col7'][i]
  dx = (x_right-x_left)*0.2
  dy = (y_right-y_left)*0.15
  dz = (z_right-z_left)*0.15
  xL2,xR2 = x_left+3*dx, x_left+4*dx
  yL2,yR2 = y_left+3*dy, y_left+4*dy
  zL2,zR2 = z_left+3*dz, z_left+4*dz
  line = str(track['col1'][i])+' '+str(xL2)+' '+' '+str(yL2)+' '+' '+str(zL2)+' '
  line = line+str(xR2)+' '+str(yR2)+' '+str(zR2)+' 11\n'
  f.write(line)

f.close()

def plotting_the_boxes():
    topleft = [x_left,y_left,z_left]
    botmright = [x_right,y_right,z_right]
    topright = [x_right,y_left,z_left]
    botmleft = [x_left,y_right,z_left]
    ray1 = ds.ray(topleft, topright)
    ray2 = ds.ray(topleft, botmleft)
    ray3 = ds.ray(topright, botmright)
    ray4 = ds.ray(botmleft, botmright)

    dx = (x_right-x_left)*0.2
    dy = (y_right-y_left)*0.15
    dz = (z_right-z_left)*0.15

    xL2,xR2 = x_left+3*dx, x_left+4*dx
    yL2,yR2 = y_left+3*dy, y_left+4*dy
    zL2,zR2 = z_left+3*dz, z_left+4*dz
    topleft2 = [xL2,yL2,zL2]
    botmright2 = [xR2,yR2,zR2]
    topright2 = [xR2,yL2,zL2]
    botmleft2 = [xL2,yR2,zL2]
    ray5 = ds.ray(topleft2, topright2)
    ray6 = ds.ray(topleft2, botmleft2)
    ray7 = ds.ray(topright2, botmright2)
    ray8 = ds.ray(botmleft2, botmright2)


    proj = yt.SlicePlot(ds,'z','Density',center=refine_box_center,width=(120.,'kpc'))
    proj.annotate_ray(ray1)
    proj.annotate_ray(ray2)
    proj.annotate_ray(ray3)
    proj.annotate_ray(ray4)
    proj.annotate_ray(ray5)
    proj.annotate_ray(ray6)
    proj.annotate_ray(ray7)
    proj.annotate_ray(ray8)
    proj.annotate_grids()
    proj.save()
    return
