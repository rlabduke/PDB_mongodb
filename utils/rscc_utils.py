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

def set_radius(d_min):
  if(d_min < 1.0):                    atom_radius = 1.0
  elif(d_min >= 1.0 and d_min<2.0):   atom_radius = 1.5
  elif(d_min >= 2.0 and d_min < 4.0): atom_radius = 2.0
  else:                               atom_radius = 2.5
  return atom_radius

def get_rscc_diff(pdb_file,reflection_file,log=None) :
#fmodel, pdb_hierarchy, params=None, log=None, show_results=False):
  #params = master_params().extract()
  #fmodel = get_fmodel(pdb_file   = pdb_file,
  #                    mtz_file   = reflection_file,
  #                    params     = params)
  #pdb_hierarchy = pdb_utils.get_pdb_hierarchy(pdb_file)
  #if(params is None): params =master_params().extract()
  #if(log is None): log = sys.stdout
  # compute map coefficients
  #import mmtbx.maps.utils
  #from libtbx.utils import null_out
  #if pdb_file.endswith('.pdb') : mfn = pdb_file.replace('.pdb','_map.mtz')
  #else : mfn = pdb_file + '_map.mtz'
  #mmtbx.maps.utils.create_map_from_pdb_and_mtz(
  #  pdb_file=pdb_file,
  #  mtz_file=reflection_file,
  #  output_file=mfn,
  #  fill=False,
  #  out=null_out())
  #assert (os.path.isfile(mfn))
  #maps_in = any_file(mfn)
  #assert (len(maps_in.file_server.miller_arrays) == 3)
  # generate_water_omit_map
  # function  to get theses
  pdb_in = any_file(pdb_file)
  hierarchy = pdb_in.file_object.hierarchy
  inputs = mmtbx.utils.process_command_line_args([pdb_file,reflection_file])
  determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
    reflection_file_server = inputs.get_reflection_file_server(),
    keep_going             = True, # don't stop if free flags are not present
    log                    = StringIO())
  f_obs = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  xrs = iotbx.pdb.input(
    file_name=inputs.pdb_file_names[0]).xray_structure_simple()
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

