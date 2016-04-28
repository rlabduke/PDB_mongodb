import sys,os
from cctbx import miller
import mmtbx.utils
from libtbx import group_args
from cctbx.array_family import flex
from cctbx import maptbx
import utils
from cStringIO import StringIO
import iotbx.pdb
from iotbx.file_reader import any_file
from iotbx import mtz

def get_data_labels(mtz_file) :
  mtz_object = mtz.object(file_name=mtz_file)
  labels = mtz_object.column_labels()
  if 'FOBS' in labels and 'SIGFOBS' in labels : return 'FOBS,SIGFOBS'
  elif 'FOBS' in labels : return 'FOBS'

def set_radius(d_min):
  if(d_min < 1.0):                    atom_radius = 1.0
  elif(d_min >= 1.0 and d_min<2.0):   atom_radius = 1.5
  elif(d_min >= 2.0 and d_min < 4.0): atom_radius = 2.0
  else:                               atom_radius = 2.5
  return atom_radius

def excise_unk_get_lines(pdb_file) :
  lines = []
  fle = open(pdb_file,'r')
  for l in fle :
    if ' UNK ' in l and l.strip().startswith('HETATM') : continue
    lines.append(l.strip())
  return lines

def get_rscc_diff(pdb_file,reflection_file,high_resolution,log=None) :
  if not log : log = sys.stderr
  print >> log, '*' * 20 + '  rscc  ' + '*' * 20
  print >> log, 'pdb file : %s' % pdb_file
  print >> log, 'hkl mtz : %s' % reflection_file
  print >> log, 'high_resolution : %.2f' % high_resolution
  pdb_in = any_file(pdb_file)
  hierarchy = pdb_in.file_object.hierarchy
  inputs = mmtbx.utils.process_command_line_args([pdb_file,reflection_file])

  # try to get data labels
  data_labels = get_data_labels(reflection_file)
  parameters = mmtbx.utils.data_and_flags_master_params().extract()
  parameters.force_anomalous_flag_to_be_equal_to = False
  if data_labels : parameters.labels = [data_labels]
  if high_resolution : parameters.high_resolution = high_resolution
  data_and_flags = mmtbx.utils.determine_data_and_flags(
    reflection_file_server = inputs.get_reflection_file_server(),
    keep_going             = True, # don't stop if free flags are not present
    parameters             = parameters,
    log                    = StringIO())
  f_obs = data_and_flags.f_obs
  r_free_flags = data_and_flags.f_obs.array(
        data = flex.bool(data_and_flags.f_obs.size(), False))
  pdb_lines = excise_unk_get_lines(inputs.pdb_file_names[0])
  xrs = iotbx.pdb.input(
    lines = pdb_lines).xray_structure_simple()
    #file_name=inputs.pdb_file_names[0]).xray_structure_simple()
  fmodel = mmtbx.utils.fmodel_simple(
    f_obs=f_obs,
    r_free_flags=r_free_flags,
    scattering_table="n_gaussian",
    xray_structures=[xrs],
    bulk_solvent_correction=True,
    skip_twin_detection=True)

  e_map_obj = fmodel.electron_density_map()
  coeffs_1 = e_map_obj.map_coefficients(
    map_type     = 'Fc',
    fill_missing = False,
    isotropize   = False)
  coeffs_2 = e_map_obj.map_coefficients(
    map_type     = '2mFo-DFc',
    fill_missing = True,
    isotropize   = True)
  coeffs_3 = e_map_obj.map_coefficients(
    map_type     = 'mFo-DFc',
    fill_missing = False,
    isotropize   = True)
  # compute maps
  fft_map_1 = coeffs_1.fft_map(resolution_factor = 1./4)
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
  print >> log, "  Overall map cc(%s,%s): %6.4f"%('Fc','2mFo-DFc',overall_cc)

  detail='atom'
  atom_radius = set_radius(d_min=fmodel.f_obs().d_min())
  print >> sys.stderr, 'radius : %.3f' % atom_radius
  results = compute(
    pdb_hierarchy        = hierarchy,
    unit_cell            = fmodel.xray_structure.unit_cell(),
    map_1                = map_1,
    map_2                = map_2,
    map_3                = map_3,
    detail               = detail,
    atom_radius          = atom_radius)
  return results

def compute(pdb_hierarchy,
            unit_cell,
            map_1,
            map_2,
            map_3,
            detail,
            atom_radius) :
  results = []
  for chain in pdb_hierarchy.chains():
    for residue_group in chain.residue_groups():
      for conformer in residue_group.conformers():
        for residue in conformer.residues():
          r_id_str = "%2s %1s %3s %4s %1s"%(chain.id, conformer.altloc,
            residue.resname, residue.resseq, residue.icode)
          for atom in residue.atoms():
            a_id_str = "%s %4s"%(r_id_str, atom.name)
            rad = atom_radius
            if not atom.element_is_hydrogen() :
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

