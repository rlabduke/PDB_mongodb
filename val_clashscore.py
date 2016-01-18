import sys
from mmtbx.validation.clashscore import clashscore
from iotbx import file_reader
import iotbx.pdb
from pdb_utils import MDB_PDB_validation

class CLASHSCOREvalidation(MDB_PDB_validation) :

  def __init__(self,pdb_file,mdb_document=None) :
    MDB_PDB_validation.__init__(self,pdb_file,mdb_document=mdb_document)
    self.run_validation()

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
    pass
