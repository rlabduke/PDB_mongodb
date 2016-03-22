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

class RSCCvalidation(object) :

  def __init__(self,pdb_file,hklmtz_file,pdb_code,meta_data) :
    self.pdb_file    = pdb_file
    self.hklmtz_file = hklmtz_file
    self.meta_data   = meta_data
    self.pdb_code    = pdb_code
    self.run_validation()
    assert os.path.exists(self.pdb_file)
    assert os.path.exists(self.hklmtz_file)

  def run_validation(self) :
    rsc = rscc_utils.get_rscc_diff(
      pdb_file=self.pdb_file,
      reflection_file=self.hklmtz_file)

    # self.mdb_residues have MDBRes.get_residue_key() keys and values are lists
    # of MDBAtom[s].
    self.mdb_residues = {}
    atomd = None
    broadcastdetail = True
    for i, result_ in enumerate(rsc) :
      resd = {'pdb_id'     : self.pdb_code,
              'model_id'   : None, 
              'chain_id'   : result_.chain_id,
              'icode'      : result_.ins_code,
              'resseq'     : result_.res_num,
              'altloc'     : result_.alt_loc,
              'resname'    : result_.res_name}
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
      if reskey not in self.mdb_residues.keys() :
        #self.mdb_residues[reskey] = MDBRes
        self.mdb_residues[reskey] = []
      if detail == 'atom' :
        #self.mdb_residues[reskey].deposit_atom(mdb_utils.MDBAtom(**atomd))
        self.mdb_residues[reskey].append(atomd)
