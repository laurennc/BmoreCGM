from emissivity_balance import *
import yt

ds = yt.load("/Users/dalek/data/Molly/natural/nref11/RD0020/RD0020")

add_line_fields(ds,['OVI_1032','CIII_977'])

ad = ds.all_data()
print ad[('gas','OVI_1032_emissivity_photons')].min(),ad[('gas','OVI_1032_emissivity_photons')].max()





proj = yt.ProjectionPlot(ds,'x',('gas','OVI_1032_emissivity_photons'),center=rb_center,width=(100.,'kpc'))
proj.set_zlim(('gas','OVI_1032_emissivity_photons'),1e-1,1e7)
proj.save('trident_emis')

proj = yt.ProjectionPlot(ds,'x',('gas','CIII_977_emissivity_photons'),center=rb_center,width=(100.,'kpc'))
proj.set_zlim(('gas','CIII_977_emissivity_photons'),1e-1,1e7)
proj.save('trident_emis')


proj = yt.ProjectionPlot(ds,'x',('gas','Emission_OVI'),center=rb_center,width=(100.,'kpc'))
proj.set_zlim(('gas','Emission_OVI'),1e-1,1e7)
proj.save('lauren_emis')

proj = yt.ProjectionPlot(ds,'x',('gas','Emission_CIII_977'),center=rb_center,width=(100.,'kpc'))
proj.set_zlim(('gas','Emission_CIII_977'),1e-1,1e7)
proj.save('lauren_emis')
