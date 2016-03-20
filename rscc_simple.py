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
from utils import mdb_utils
from utils import pdb_utils
from utils import rscc_utils
from utils import utils

def get_rscc_mdb_residues(pdb_code,log=None) :
  # crreate and switch to temp directory
  cwd = os.getcwd()
  tempdir = "%s%i" % (pdb_code,sys.maxint)
  os.makedirs(tempdir)
  tempdir = os.path.abspath(tempdir)
  os.chdir(tempdir)

  moddata = pdb_utils.ModelData(pdb_code)
  pdb_file = moddata.pdb_file
  reflection_file = moddata.mtz_file
  resolution = moddata.d_min
  assert os.path.exists(pdb_file)
  assert os.path.exists(reflection_file)

  rsc = rscc_utils.get_rscc_diff(
    pdb_file=pdb_file,
    reflection_file=reflection_file)

  mdb_residues = {}
  atomd = None
  broadcastdetail = True
  for i, result_ in enumerate(rsc) :
    resd = {'pdb_id'     : pdb_code,
            'model_id'   : None, 
            'chain_id'   : result_.chain_id,
            'icode'      : result_.ins_code,
            'resseq'     : result_.res_num,
            'altloc'     : result_.alt_loc,
            'resname'    : result_.res_name,
            'resolution' : resolution}
    if 'residue' in dir(result_) :
      if broadcastdetail :
        utils.broadcast('detail : residue')
        broadcastdetail = False
      detail = 'residue'
    elif 'atom' in dir(result_) :
      if broadcastdetail :
        utils.broadcast('detail : atom')
        broadcastdetail = False
      detail = 'atom'
      # map_value_1 = Fc
      # map_value_2 = 2mFo-DFc
      # map_value_3 = mFo-DFc
      atomd = {'name':result_.atom.name,
               'adp':result_.atom.b,
               'xyz':result_.atom.xyz,
               'rscc':result_.cc,
               'twoFo_DFc_value':result_.map_value_2,
               'Fo_DFc_value':result_.map_value_3,
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

  # add worst_bb, worse_sc, and worse_all
  if log :
    keys = mdb_residues.keys()
    keys.sort()
    for k in keys :
      mdb_utils.print_json_pretty(mdb_residues[k].get_residue_mongodoc(),log)
  # clean up
  print >> sys.stderr, "Clening up ..."
  os.chdir(cwd)
  for fn in os.listdir(tempdir) :
    os.remove(os.path.join(tempdir,fn))
  os.rmdir(tempdir)

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
    absoutpath = os.path.abspath(args.out_directory)
    outfn = os.path.join(absoutpath,"%s.rsc" % args.pdb_id)
    fle = open(outfn,'w')
  else : fle = sys.stdout
  mdb_residues = get_rscc_mdb_residues(pdb_code = args.pdb_id,log=fle)
  if args.write_out_file :
    fle.close()
    print >> sys.stderr, "%s written." % outfn

if __name__ == '__main__' :
  run(sys.argv[1:])

