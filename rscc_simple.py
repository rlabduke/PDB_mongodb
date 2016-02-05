import os,sys
import argparse
from mmtbx import real_space_correlation
import mmtbx.utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
import iotbx.pdb
from libtbx import group_args
from cctbx.array_family import flex
from cStringIO import StringIO
import mdb_utils
import pdb_utils

def broadcast(m, log = sys.stdout) :
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

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

def get_fmodel(pdb_file,mtz_file,log=sys.stderr) :
  params = mmtbx.real_space_correlation.master_params().extract()
  params.pdb_file_name = pdb_file
  params.reflection_file_name = mtz_file
  broadcast(m='Data labels : %s' % params.data_labels)
  broadcast(m="Input PDB file name: %s"%params.pdb_file_name, log=log)
  pdbo = pdb_to_xrs(pdb_file_name=params.pdb_file_name,
    scattering_table=params.scattering_table)
  pdbo.xray_structure.show_summary(f=log, prefix="  ")
  broadcast(m="Input reflection file name: %s"%params.reflection_file_name,
    log=log)
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

def get_rscc_mdb_residues(pdb_code,log=None) :
  pdb_files_dict = pdb_utils.get_pdb_files(pdb_code)
  pdb_file = pdb_files_dict["pdb"]
  reflection_file = pdb_files_dict["hklmtz"]
  assert os.path.exists(pdb_file)
  assert os.path.exists(reflection_file)
  rsc_params = real_space_correlation.master_params().extract()
  fmodel = get_fmodel(pdb_file,reflection_file)
  pdb_hierarchy = pdb_utils.get_pdb_hierarchy(pdb_file)
  rsc = real_space_correlation.simple(
    fmodel=fmodel,
    pdb_hierarchy=pdb_hierarchy,
    params=rsc_params,
    log=sys.stderr)
  mdb_residues = {}
  atomd = None
  broadcastdetail = True
  for i, result_ in enumerate(rsc) :
#   print dir(result_);exit()
    chain_id    = result_.id_str[:2].strip()
    altloc      = result_.id_str[3].strip()
    resname     = result_.id_str[5:8].strip()
    resseq      = result_.id_str[9:13].strip()
    icode       = result_.id_str[14].strip()
    resd = {'pdb_id'  : pdb_code,
            'model_id' : None, 
            'chain_id' : chain_id,
            'icode'    : icode,
            'resseq'   : resseq,
            'altloc'   : altloc,
            'resname'  : resname}
    if 'residue' in dir(result_) :
      if broadcastdetail :
        broadcast('detail : residue')
        broadcastdetail = False
      detail = 'residue'
    elif 'atom' in dir(result_) :
      if broadcastdetail :
        broadcast('detail : atom')
        broadcastdetail = False
      detail = 'atom'
#'rscc','twoFo_DFc_value','Fo_DFc_value'
      atomd = {'name':result_.atom.name,
               'adp':result_.atom.b,
               'xyz':result_.atom.xyz,
               'rscc':result_.cc,
               'occ':result_.atom.occ}
    else :
      raise RuntimeError("Could not resolve detail")
    MDBRes = mdb_utils.MDBResidue(**resd)
    reskey = MDBRes.get_residue_key()
    if reskey not in mdb_residues.keys() : mdb_residues[reskey] = MDBRes
    if detail == 'atom' :
      mdb_residues[reskey].deposit_atom(mdb_utils.MDBAtom(**atomd))
      #MDBRes.deposit_atom(mdb_utils.MDBAtom(**atomd))
    #if log : mdb_utils.print_json_pretty(MDBRes.get_residue_mongodoc(),log)
  # mdb_residues is a dict where the keys are unique strings formed by 
  # concatenating pdb_id, model_id, chain_id, icode, resseq, altloc and resname
  # and values are mdb_utils.MDBResidue.

  if log :
    keys = mdb_residues.keys()
    keys.sort()
    for k in keys :
      mdb_utils.print_json_pretty(mdb_residues[k].get_residue_mongodoc(),log)
  # clean up
  print >> sys.stderr, "Clening up ..."
  for t,fn in pdb_files_dict.items() :
    os.remove(fn)

  return mdb_residues

def run(args) :
  desc = "Run real_space_correlation given a pdb id"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_id', help='A pdb code')
  parser.add_argument('-o','--write_out_file',help='write out file',
                      action='store_true')
  parser.add_argument('-d','--out_directory',help='out directory')
  args = parser.parse_args()
  assert len(args.pdb_id) == 4
  if args.write_out_file :
    assert args.out_directory, "Must specify an out directory - see help (-h)."
    err ="Specified out directory doesn't exist"
    assert os.path.exists(args.out_directory), err
    outfn = os.path.join(args.out_directory,"%s.rsc" % args.pdb_id)
    fle = open(outfn,'w')
  else : fle = None
  mdb_residues = get_rscc_mdb_residues(pdb_code = args.pdb_id,log=fle)
  if args.write_out_file :
    fle.close()
    print >> sys.stderr, "%s written." % outfn

if __name__ == '__main__' :
  run(sys.argv[1:])

