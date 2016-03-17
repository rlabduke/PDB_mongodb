import sys
from mmtbx.validation.clashscore import clashscore
from iotbx import file_reader
import iotbx.pdb
from utils.pdb_utils import MDB_PDB_validation
import utils.mdb_utils

class CLASHSCOREvalidation(object) :

  def __init__(self,pdb_file,detail,mdb_document) :
    self.pdb_file    = pdb_file
    self.detail      = detail
    self.mdb_document= mdb_document
    self.run_validation()
    if self.detail == 'file' : self.add_file()
    elif self.detail == 'residue' : self.add_residue()

  def run_validation(self) :
    out = sys.stderr
    keep_hydrogens=True
    nuclear=False
    verbose=True
    quiet=False
    pdb_in = file_reader.any_file(self.pdb_file)
    pdb_in.check_file_type('pdb')
    input_pdb = iotbx.pdb.input(file_name=self.pdb_file)
    pdb_hierarchy = input_pdb.construct_hierarchy()
    out=sys.stdout
    self.result = clashscore(
      pdb_hierarchy=pdb_hierarchy,
      keep_hydrogens=keep_hydrogens,
      nuclear=nuclear,
      out=out,
      verbose=verbose and not quiet)
 
  def add_file(self) :
    self.mdb_document['clashscore'] = self.result.get_clashscore()
  
  def add_residue(self) :
    residues = [None,None]
    for clash in self.result.results :
      for i,atom in enumerate(clash.atoms_info) :
        resd = {'pdb_id'     : self.mdb_document['_id'],
                'model_id'   : None,
                'chain_id'   : atom.chain_id,
                'icode'      : atom.icode,
                'resseq'     : atom.resseq,
                'altloc'     : atom.altloc,
                'resname'    : atom.resname,
                'resolution' : self.mdb_document['Resolution']}
        res = utils.mdb_utils.MDBResidue(**resd)
        residues[i] = res

