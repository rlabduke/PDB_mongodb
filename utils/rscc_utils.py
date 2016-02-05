import sys
import iotbx.phil
from cctbx import miller
import mmtbx.utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
import iotbx.pdb
from libtbx import group_args
from cctbx.array_family import flex
from cStringIO import StringIO
from cctbx import maptbx
import utils
import pdb_utils
from iotbx import mtz

core_params_str = """\
atom_radius = None
  .type = float
  .help = Atomic radius for map CC calculation. Determined automatically if \
          if None is given
  .expert_level = 2
hydrogen_atom_radius = None
  .type = float
  .help = Atomic radius for map CC calculation for H or D.
  .expert_level = 2
resolution_factor = 1./4
  .type = float
use_hydrogens = None
  .type = bool
"""

master_params_str = """\
%s
scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations
detail = atom residue *automatic
  .type = choice(multi=False)
  .help = Level of details to show CC for
map_1
  .help = First map to use in map CC calculation
{
 type = Fc
   .type = str
   .help = Electron density map type. Example xmFobs-yDFcalc (for \
           maximum-likelihood weighted map) or xFobs-yFcalc (for simple \
           unweighted map), x and y are any real numbers.
 fill_missing_reflections = False
   .type = bool
 isotropize = False
   .type = bool
}
map_2
  .help = Second map to use in map CC calculation
{
 type = 2mFo-DFc
   .type = str
   .help = Electron density map type. Example xmFobs-yDFcalc (for \
           maximum-likelihood weighted map) or xFobs-yFcalc (for simple \
           unweighted map), x and y are any real numbers.
 fill_missing_reflections = True
   .type = bool
 isotropize = True
   .type = bool
}
pdb_file_name = None
  .type = str
  .help = PDB file name.
reflection_file_name = None
  .type = str
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
data_labels = None
  .type = str
  .help = Labels for experimental data.
high_resolution = None
  .type=float
low_resolution = None
  .type=float
pdb_id  = None
  .type = str
  .help = the pdb id for csv output
show_csv = False
   .type = bool
   .help = write csv to stdout
show_default_output = False
   .type = bool
"""%core_params_str

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def extract_data_and_flags(params, crystal_symmetry=None):
  data_and_flags = None
  if(params.reflection_file_name is not None):
    reflection_file = reflection_file_reader.any_reflection_file(
      file_name = params.reflection_file_name)
    reflection_file_server = reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry   = True,
      reflection_files = [reflection_file])
    parameters = mmtbx.utils.data_and_flags_master_params().extract()
    parameters.force_anomalous_flag_to_be_equal_to = False
    if(params.data_labels is not None):
      parameters.labels = [params.data_labels]
    #### Not Relevant
    # if(params.high_resolution is not None):
      # parameters.high_resolution = params.high_resolution
    # if(params.low_resolution is not None):
      # parameters.low_resolution = params.low_resolution
    #### Not Relevant #### 
    data_and_flags = mmtbx.utils.determine_data_and_flags(
      reflection_file_server = reflection_file_server,
      parameters             = parameters,
      data_description       = "X-ray data",
      extract_r_free_flags   = False, # XXX
      log                    = StringIO())
  return data_and_flags

def pdb_to_xrs(pdb_file_name, scattering_table):
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  xray_structure = pdb_inp.xray_structure_simple()
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq() # VERY important to do.
  mmtbx.utils.setup_scattering_dictionaries(
    scattering_table = scattering_table,
    xray_structure = xray_structure,
    d_min = None)
  return group_args(
    xray_structure = xray_structure,
    pdb_hierarchy  = pdb_hierarchy)

def get_data_labels(mtz_file) :
  mtz_object = mtz.object(file_name=mtz_file)
  labels = mtz_object.column_labels()
  if 'FOBS' in labels and 'SIGFOBS' in labels : return 'FOBS,SIGFOBS'

def get_fmodel(pdb_file,mtz_file,params,log=sys.stderr) :
  params.pdb_file_name = pdb_file
  params.reflection_file_name = mtz_file
  params.data_labels = get_data_labels(mtz_file)
  utils.broadcast(m='Data labels : %s' % params.data_labels)
  utils.broadcast(m="Input PDB file name: %s"%params.pdb_file_name, log=log)
  pdbo = pdb_to_xrs(pdb_file_name=params.pdb_file_name,
    scattering_table=params.scattering_table)
  pdbo.xray_structure.show_summary(f=log, prefix="  ")
  m = "Input reflection file name: %s"
  utils.broadcast(m=m % params.reflection_file_name, log=log)
  data_and_flags = extract_data_and_flags(params = params)
  data_and_flags.f_obs.show_comprehensive_summary(f=log, prefix="  ")
  # create fmodel
  r_free_flags = data_and_flags.f_obs.array(
    data = flex.bool(data_and_flags.f_obs.size(), False))
  fmodel = mmtbx.utils.fmodel_simple(
    xray_structures     = [pdbo.xray_structure],
    scattering_table    = params.scattering_table,
    f_obs               = data_and_flags.f_obs,
    r_free_flags        = r_free_flags)
  return fmodel

def set_radius(d_min):
  if(d_min < 1.0):                    atom_radius = 1.0
  elif(d_min >= 1.0 and d_min<2.0):   atom_radius = 1.5
  elif(d_min >= 2.0 and d_min < 4.0): atom_radius = 2.0
  else:                               atom_radius = 2.5
  return atom_radius

def get_rscc_diff(pdb_file,reflection_file,log=None) :
#fmodel, pdb_hierarchy, params=None, log=None, show_results=False):
  params = master_params().extract()
  fmodel = get_fmodel(pdb_file   = pdb_file,
                      mtz_file   = reflection_file,
                      params     = params)
  pdb_hierarchy = pdb_utils.get_pdb_hierarchy(pdb_file)
  if(params is None): params =master_params().extract()
  if(log is None): log = sys.stdout
  # compute map coefficients
  e_map_obj = fmodel.electron_density_map()
  coeffs_1 = e_map_obj.map_coefficients(
    map_type     = params.map_1.type,
    fill_missing = params.map_1.fill_missing_reflections,
    isotropize   = params.map_1.isotropize)
  coeffs_2 = e_map_obj.map_coefficients(
    map_type     = params.map_2.type,
    fill_missing = params.map_2.fill_missing_reflections,
    isotropize   = params.map_2.isotropize)
  coeffs_3 = e_map_obj.map_coefficients(
    map_type     = 'mFo-DFc',
    fill_missing = True,
    isotropize   = True)
  # compute maps
  fft_map_1 = coeffs_1.fft_map(resolution_factor = params.resolution_factor)
  fft_map_1.apply_sigma_scaling()
  map_1 = fft_map_1.real_map_unpadded()
  fft_map_2 = miller.fft_map(
    crystal_gridding     = fft_map_1,
    fourier_coefficients = coeffs_2)
  fft_map_2.apply_sigma_scaling()
  map_2 = fft_map_2.real_map_unpadded()
  fft_map_3 = miller.fft_map(
    crystal_gridding     = fft_map_1,
    fourier_coefficients = coeffs_3)
  fft_map_3.apply_sigma_scaling()
  map_3 = fft_map_3.real_map_unpadded()
  # compute cc
  utils.broadcast(m="Map correlation and map values", log=log)
  overall_cc = flex.linear_correlation(x = map_1.as_1d(),
    y = map_2.as_1d()).coefficient()
  print >> log, "  Overall map cc(%s,%s): %6.4f"%(params.map_1.type,
    params.map_2.type, overall_cc)

  detail='atom'
  atom_radius = set_radius(d_min=fmodel.f_obs().d_min())
  use_hydrogens = params.use_hydrogens
  if(use_hydrogens is None):
    if(params.scattering_table == "neutron" or fmodel.f_obs().d_min() <= 1.2):
      use_hydrogens = True
    else:
      use_hydrogens = False
  hydrogen_atom_radius = params.hydrogen_atom_radius
  if(hydrogen_atom_radius is None):
    if(params.scattering_table == "neutron"):
      hydrogen_atom_radius = atom_radius
    else:
      hydrogen_atom_radius = 1
  results = compute(
    pdb_hierarchy        = pdb_hierarchy,
    unit_cell            = fmodel.xray_structure.unit_cell(),
    map_1                = map_1,
    map_2                = map_2,
    map_3                = map_3,
    detail               = detail,
    atom_radius          = atom_radius,
    use_hydrogens        = use_hydrogens,
    hydrogen_atom_radius = hydrogen_atom_radius)
  return results

def compute(pdb_hierarchy,
            unit_cell,
            map_1,
            map_2,
            map_3,
            detail,
            atom_radius,
            use_hydrogens,
            hydrogen_atom_radius):
  results = []
  for chain in pdb_hierarchy.chains():
    for residue_group in chain.residue_groups():
      for conformer in residue_group.conformers():
        for residue in conformer.residues():
          r_id_str = "%2s %1s %3s %4s %1s"%(chain.id, conformer.altloc,
            residue.resname, residue.resseq, residue.icode)
          for atom in residue.atoms():
            a_id_str = "%s %4s"%(r_id_str, atom.name)
            if(atom.element_is_hydrogen()): rad = hydrogen_atom_radius
            else: rad = atom_radius
            if(not (atom.element_is_hydrogen() and not use_hydrogens)):
              map_value_1 = map_1.eight_point_interpolation(
                unit_cell.fractionalize(atom.xyz))
              map_value_2 = map_2.eight_point_interpolation(
                unit_cell.fractionalize(atom.xyz))
              map_value_3 = map_3.eight_point_interpolation(
                unit_cell.fractionalize(atom.xyz))
              sel = maptbx.grid_indices_around_sites(
                unit_cell  = unit_cell,
                fft_n_real = map_1.focus(),
                fft_m_real = map_1.all(),
                sites_cart = flex.vec3_double([atom.xyz]),
                site_radii = flex.double([atom_radius]))
              cc = flex.linear_correlation(x=map_1.select(sel),
                y=map_2.select(sel)).coefficient()
              result = group_args(
                chain_id    = chain.id,
                res_num     = residue.resseq,
                res_name    = residue.resname,
                atom_name   = atom.name,
                alt_loc     = conformer.altloc,
                ins_code    = residue.icode,
                atom        = atom,
                id_str      = a_id_str,
                cc          = cc,
                map_value_1 = map_value_1,
                map_value_2 = map_value_2,
                map_value_3 = map_value_3,
                b           = atom.b,
                occupancy   = atom.occ,
                n_atoms     = 1)
              results.append(result)
  return results

