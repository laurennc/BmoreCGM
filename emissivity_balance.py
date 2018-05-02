"""
Emissivity fields using Cloudy data.

"""

from yt.fields.field_detector import \
    FieldDetector
from yt.fields.local_fields import add_field
from yt.utilities.linear_interpolators import \
    TrilinearFieldInterpolator, \
    UnilinearFieldInterpolator
import six
from yt.utilities.physical_constants import mh
from yt.funcs import mylog
import numpy as np
import h5py
import copy
import os
import warnings
import yt.units as u
from trident.config import \
    ion_table_filepath
from trident.line_database import \
    LineDatabase, \
    uniquify
from trident.utilities import \
    _determine_dataset_sampling_type
from trident.roman import \
    from_roman
from trident.ion_balance import \
    _log_nH, _log_T, _redshift, _determine_sampling_type

### USEFUL CONSTANTS AND UNITS ###
h = 6.626*10**(-27) * u.erg* u.s
c = 3e8 * u.m / u.s
###
emission_units_photons = 's**-1 * cm**-3 * steradian**-1'
ytEmU_photons = u.ergs * u.s**-1 * u.cm**-3 * u.steradian**-1
####
emission_units_ergs = 'erg * s**-1 * cm**-3 * arcsec**-2'
ytEmU_ergs = u.ergs * u.s**-1 * u.cm**-3 * u.steradian**-1


### NEED TO ADD THIS TO THE TRIDENT CONFIGURE FILE
emis_table_filepath = "/Users/dalek/repos/BmoreCGM/emissivity_tables.h5"

# NOTE:
# I don't set a floor in emissivity that's zero because it affects all of the
# projections and makes taking logs much harder. Instead, the floor is set
# by the double precision limit in Cloudy

table_store = {}

class EmissivityTable(object):
    def __init__(self,filename=None,line=None):
        """
        Base class for adding emissivity fields

        Used to load in an HDF5 file that contains the values for the
        line emissivity as a function of density, temperature, and redshift
        (for the extragalactic UV background). Currently, the emissivity
        gets scaled linearly by metallicity.

        **Parameters**

        :filename: string, option

            Name of the HDF5 file that contains the emissitivy table data.

            Default: it uses the table specified in ~/.trident/config

        :line: string, optional

            The line for which you want to create an EmissivityTable

            Default: None
        """
        if filename is None:
            filename = emis_table_filepath
        self.filename = filename
        self.parameters = []
        self.line_emissivity = []
        self._load_hdf5_emis_table(line)

    def _load_hdf5_emis_table(self,line):
        """
        Read in the HDF5 ion balance table
        """

        input =  h5py.File(self.filename,'r')
        self.line_emissivity = input[line].value
        #### WHY DONT I JUST ACTUALLY NAME THE PARAMETERS???
        ## not -1 here because we don't have the extra column for number of ions
        for par in range(1,len(self.line_emissivity.shape)):
            name = "Parameter%d" % par
            self.parameters.append(input[line].attrs[name])
        self.parameters.append(input[line].attrs['Temperature'])
        input.close()


### THEN THERE ARE THE METHODS THAT BUILD THE HELPER FIELDS
# def _log_nH(field, data):
# def _redshift(field, data):
# def _log_T(field, data):
# def _determine_sampling_type(ds, sampling_type, particle_type):
### THESE SHOULD BE A STRAIGHT COPY FROM ORIGINAL TRIDENT

def add_line_fields(ds,lines,ftype='gas',
                    emissivity_table=None,
                    field_suffix=False,
                    force_override=False,
                    sampling_type='auto',
                    particle_type=None):
    """
    Preferred method for adding line emissivity fields to a yt dataset.

    Line emissivities need to be generated on a line-by-line basis in Cloudy.
    The list of available lines can be found by trying:

    ### PUT LINES AND THEIR ENERGIES IN A FILE SOMEWHERE!
    ### THEN REFERENCE HERE!

    For each line selected, two fields will be added for two different
    units systems (example for CIII 977):
        * photon s^-1 cm^-2 sr^-1 (ftype, 'CIII_977_emissivity_photons')
        * ergs s^-1 cm^-2 arcsec^-2 (ftype, 'CIII_977_emissivity_ergs')

    The built in yt units can handle conversions between sr/arcsec and
    other standard units. These fields are to facilitate line-specific
    energies.

    Fields are added assuming collisional ionization equilibrium and
    photoionization in the optically thin limit from a redshift-dependent
    metagalactic ionizing background using the ionization_table specified.

    **WARNING**: The "ftype" must match the field type that you're using for
    the field interpolation.  So for particle-based codes, this must be the
    ftype of the gas particles (e.g., `PartType0`, `Gas`).  Using the
    default of `gas` in this instance will interpolate on the grid-based
    fields, which will give the wrong answers for particle-based codes,
    since the ion field interpolation will take place on the already
    deposited grid-based fields.

    **Parameters**

    :ds: yt dataset object

        This is the dataset to which the line emissivity field will be added.

    :lines: list of strings

        List of strings matching possible emission lines. The list of available
        lines can be found by MORE HERE

        If set to 'all', creates **all** lines available in the database.

    :ftype: string, option

        The field type of the field to add. It is the first string in the
        field tuple e.g. "gas" in ("gas", "O_p5_ion_fraction")
        ftype must correspond to the ftype of the 'density', and 'temperature'
        fields in your dataset you wish to use to generate the ion field.
        Default: "gas"

    :emissivity_table: string, optional

        Path to an appropriately formatted HDF5 table that can be used to
        compute the ion fraction as a function of density, temperature,
        and redshift.  When set to None, it uses the table
        specified in ~/.trident/config
        Default: None

    :field_suffix: boolean, optional

        Determines whether or not to append a suffix to the field name that
        indicates what ionization table was used.  Useful when using generating
        ion_fields that already exist in a dataset.

    :force_override: boolean, optional

        Set to True if you wish to clobber existing ion fields with any
        created with this functionality.  Otherwise, existing fields will
        remain untouched.
        Default: False

    :sampling_type: string, optional

        Set to 'particle' if the field should be for particles.
        Set to 'cell' if the field should be for grids/cells.
        Set to 'auto' for this to be determined automatically.
        Default: 'auto'

    :particle_type: boolean, optional

        This is deprecated in favor of 'sampling_type'.
        Set to True if you are adding ion fields to particles, as specified
        by the 'ftype'. Set to False if you are not. Set to 'auto', if
        you the want the code to autodetermine if the field specified by the
        'ftype' is particle or not.
        Default: 'auto'

    **Example**
    To add H-alpha and CIII 977 to a dataset, you would run:

    >>> import yt
    >>> import trident
    >>> ds =  yt.load('path/to/file')
    >>> trident.add_emissivity_fields(ds,lines=['HI_6563','CIII_977'])
    """

    #line_list = []
    line_list = lines
    sampling_type = \
        _determine_sampling_type(ds,sampling_type,particle_type)

    if emissivity_table is None:
        emissivity_table = emis_table_filepath

    ## No line parsing needed because we're making the user provide
    ## lines in the correct form since they're so specific

    # make sure the line list is unique
    line_list = uniquify(line_list)

    for line in line_list:
        add_line_emissivity_fields(line,ds,ftype,emissivity_table,
                                   field_suffix=field_suffix, force_override=force_override,
                                   sampling_type=sampling_type)

def add_line_emissivity_fields(line,ds,ftype="gas",emissivity_table=None,
                               field_suffix=False,force_override=False,
                               sampling_type="auto",particle_type=None):
    """
    Add line emissivity fields to a yt dataset for the desired line.

    ..note::

        The preferred method for add line emissivity fields to a datset is
        using :class:`~trident.add_line_fields`

        This is a subroutine used within that method.

    Fields are added assuming collisional ionization equilibrium and
    photoionization in the optically thin limit from a redshift-dependent
    metagalactic ionizing background using the ionization_table specified.

    **WARNING**: The "ftype" must match the field type that you're using for
    the field interpolation.  So for particle-based codes, this must be the
    ftype of the gas particles (e.g., `PartType0`, `Gas`).  Using the
    default of `gas` in this instance will interpolate on the grid-based
    fields, which will give the wrong answers for particle-based codes,
    since the ion field interpolation will take place on the already
    deposited grid-based fields.

    **Parameters**

    :line: string
        Line for which the emissivity is desired. Ion and
        wavelength in angstroms must be specified (e.g. 'CIII_977')

    :ds: yt dataset object
        This is the dataset to which the line emissivity field will be added.

    :ftype: string, optional
        The field type of the field to add.  it is the first string in the
        field tuple e.g. "gas" in ("gas", "O_p5_ion_fraction")
        ftype must correspond to the ftype of the 'density', and 'temperature'
        fields in your dataset you wish to use to generate the ion field.
        Default: "gas"

    :emissivity_table: string, optional
        Path to an appropriately formatted HDF5 table that can be used to
        compute the ion fraction as a function of density, temperature,
        and redshift.  By default, it uses the table specified in
        ~/.trident/config

    :field_suffix: boolean, optional
        Determines whether or not to append a suffix to the field name that
        indicates what ionization table was used

    :force_override: boolean, optional

        Set to True if you wish to clobber existing ion fields with any
        created with this functionality.  Otherwise, existing fields will
        remain untouched.
        Default: False

    :sampling_type: string, optional

        Set to 'particle' if the field should be for particles.
        Set to 'cell' if the field should be for grids/cells.
        Set to 'auto' for this to be determined automatically.
        Default: 'auto'

    :particle_type: boolean, optional

        This is deprecated in favor of 'sampling_type'.
        Set to True if you are adding ion fields to particles, as specified
        by the 'ftype'.  Set to False if you are not.  Set to 'auto', if
        you want the code to autodetermine if the field specified by the
        'ftype' is particle or not.
        Default: 'auto'

    """

    sampling_type = \
      _determine_sampling_type(ds, sampling_type, particle_type)

    if emissivity_table is None:
        emissivity_table = emis_table_filepath

    if (ftype, "log_nH") not in ds.derived_field_list:
        ds.add_field((ftype, "log_nH"), function=_log_nH, units="",
                     sampling_type=sampling_type,
                     force_override=force_override)

    if (ftype, "redshift") not in ds.derived_field_list:
        ds.add_field((ftype, "redshift"), function=_redshift, units="",
                     sampling_type=sampling_type,
                     force_override=force_override)

    if (ftype, "log_T") not in ds.derived_field_list:
        ds.add_field((ftype, "log_T"), function=_log_T, units="",
                     sampling_type=sampling_type,
                     force_override=force_override)

    photon_field = line+'_emissivity_photons'
    ergs_field   = line+'_emissivity_ergs'

    if not line in table_store:
        lineTable = EmissivityTable(emissivity_table,line)
        table_store[line] = {'emissivity':copy.deepcopy(lineTable.line_emissivity),
                             'parameters':copy.deepcopy(lineTable.parameters)}
        del lineTable

    ds.add_field((ftype,photon_field), function=_line_emis_photons_field,
                 units=emission_units_photons,
                 sampling_type=sampling_type, force_override=force_override)

    ds.add_field((ftype,ergs_field), function=_line_emis_ergs_field,
                 units=emission_units_ergs,
                 sampling_type=sampling_type, force_override=force_override)

    if sampling_type == 'particle':
        new_fieldP = ds.add_smoothed_particle_field((ftype,photon_field))
        new_fieldE = ds.add_smoothed_particle_field((ftype,ergs_field))
        if ftype != "gas":
            ds.field_info.alias(("gas",alias_fieldP),new_fieldP)
            ds.derived_field_list.append(("gas",alias_fieldP))
            ds.field_info.alias(("gas",alias_fieldE),new_fieldE)
            ds.derived_field_list.append(("gas",alias_fieldE),new_fieldE)
    return

#### DON'T FORGET TO ADD THE 4PI FOR STERADIANS!!!!
def _line_emis_ergs_field(field,data):
    if isinstance(field.name,tuple):
        ftype = field.name[0]
        field_name = field.name[1]
    else:
        ftype="gas"
        field_name = field.name

    line_parts = field_name.split('_')
    table_line = line_parts[0]+('_')+line_parts[1]
    n_parameters = len(table_store[table_line]['parameters'])

    if n_parameters == 1:
        lineEmis = table_store[table_line]['emissivity']
        t_param = table_store[table_line]['parameters'][0]
        bds = t_param.astype("=f8")

        interp = UnilinearFieldInterpolator(lineEmis,bds,'log_T',truncate=True)

    elif n_parameters == 3:
        lineEmis = table_store[table_line]['emissivity']
        n_param = table_store[table_line]['parameters'][0]
        z_param = table_store[table_line]['parameters'][1]
        t_param = table_store[table_line]['parameters'][2]
        bds = [n_param.astype("=f8"),z_param.astype("=f8"),t_param.astype("=f8")]

        interp = TrilinearFieldInterpolator(lineEmis,bds,
                                            [(ftype,"log_nH"),
                                             (ftype,"redshift"),
                                             (ftype, "log_T")],
                                             truncate=True)
    else:
        raise RuntimeError("This data file format is not supported.")

    ##interpolate first then multiply by density**2 later
    ## MORE HERE ##
    ## THIS may be a place where I can save time!

    hden = np.power(10,data["log_nH"])
    emissivity = np.power(10,interp(data))*(hden**2.)*ytEmU_ergs/(4.*np.pi)

    if line_parts[0][1].isupper():
        atom = line_parts[0][0]
    else:
        atom = line_parts[0][0:2]

    #### NEED TO SCALE BY THE METALLICITY
    # try the species metallicity
    metallicity_field = "%s_metallicity" % atom
    if (ftype,metallicity_field) in data.ds.field_info:
        total_Z = data[ftype,metallicity_field]/solar_abundance[atom]
        emissivity = emissivity * total_Z
        return emissivity.to(emission_units_ergs)

    if atom == 'H' or atom == 'He':
        return emissivity.to(emission_units_ergs)
    else:
        emissivity = emissivity*data[ftype,"metallicity"]/(data.ds.quan(1.,'Zsun'))
        return emissivity.to(emission_units_ergs)



def _line_emis_photons_field(field,data):
    if isinstance(field.name,tuple):
        ftype = field.name[0]
        field_name = field.name[1]
    else:
        ftype="gas"
        field_name = field.name

    line_parts = field_name.split('_')
    table_line = line_parts[0]+('_')+line_parts[1]
    n_parameters = len(table_store[table_line]['parameters'])

    if n_parameters == 1:
        lineEmis = table_store[table_line]['emissivity']
        t_param = table_store[table_line]['parameters'][0]
        bds = t_param.astype("=f8")

        interp = UnilinearFieldInterpolator(lineEmis,bds,'log_T',truncate=True)

    elif n_parameters == 3:
        lineEmis = table_store[table_line]['emissivity']
        n_param = table_store[table_line]['parameters'][0]
        z_param = table_store[table_line]['parameters'][1]
        t_param = table_store[table_line]['parameters'][2]
        bds = [n_param.astype("=f8"),z_param.astype("=f8"),t_param.astype("=f8")]

        interp = TrilinearFieldInterpolator(lineEmis,bds,
                                            [(ftype,"log_nH"),
                                             (ftype,"redshift"),
                                             (ftype, "log_T")],
                                             truncate=True)
    else:
        raise RuntimeError("This data file format is not supported.")

    hden = np.power(10,data["log_nH"])
    emissivity = np.power(10,interp(data))
    emissivity = emissivity*(hden**2.)*ytEmU_photons/(4.*np.pi)
    ##### NEED TO DIVIDE BY THE LINE ENERGY
    line_parts = field_name.split('_')
    lEnergy = emission_line_energy(float(line_parts[1]),u.angstrom)
    emissivity = emissivity / lEnergy

    ## Figure out which atom we're dealing with
    if line_parts[0][1].isupper():
        atom = line_parts[0][0]
    else:
        atom = line_parts[0][0:2]
    #### NEED TO SCALE BY THE METALLICITY
    # try the species metallicity
    metallicity_field = "%s_metallicity" % atom
    if (ftype,metallicity_field) in data.ds.field_info:
        total_Z = data[ftype,metallicity_field]/solar_abundance[atom]
        return emissivity * total_Z

    if atom == 'H' or atom == 'He':
        return emissivity
    else:
        return emissivity*data[ftype,"metallicity"]/(data.ds.quan(1.,'Zsun'))



solar_abundance = {
    'H' : 1.00e+00, 'He': 1.00e-01, 'Li': 2.04e-09,
    'Be': 2.63e-11, 'B' : 6.17e-10, 'C' : 2.45e-04,
    'N' : 8.51e-05, 'O' : 4.90e-04, 'F' : 3.02e-08,
    'Ne': 1.00e-04, 'Na': 2.14e-06, 'Mg': 3.47e-05,
    'Al': 2.95e-06, 'Si': 3.47e-05, 'P' : 3.20e-07,
    'S' : 1.84e-05, 'Cl': 1.91e-07, 'Ar': 2.51e-06,
    'K' : 1.32e-07, 'Ca': 2.29e-06, 'Sc': 1.48e-09,
    'Ti': 1.05e-07, 'V' : 1.00e-08, 'Cr': 4.68e-07,
    'Mn': 2.88e-07, 'Fe': 2.82e-05, 'Co': 8.32e-08,
    'Ni': 1.78e-06, 'Cu': 1.62e-08, 'Zn': 3.98e-08}

def emission_line_energy(wavelength,unit):
    E = h*c/(wavelength*unit)
    return E.to('erg')
