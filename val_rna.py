from mmtbx.validation.rna_validate import rna_validation
from mmtbx.monomer_library import server
from mmtbx.monomer_library import pdb_interpretation
import iotbx.pdb
from utils.pdb_utils import MDB_PDB_validation

class RNAvalidation(object) :

  def __init__(self,pdb_file,detail,mdb_document) :
    self.pdb_file    = pdb_file
    self.detail      = detail
    self.mdb_document= mdb_document
    self.run_validation()
    if self.detail == 'file' : self.add_file()
    elif self.detail == 'residue' : self.add_residue()

  def run_validation(self) :
    # get required elements to run:
    #   - pdb_hierarchy and geometry
    pdb_in = iotbx.pdb.input(file_name=self.pdb_file)
    mon_lib_srv = server.server()
    ener_lib = server.ener_lib()
    processed_pdb_file = pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      pdb_inp=pdb_in,
      substitute_non_crystallographic_unit_cell_if_necessary=True)
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    geometry = processed_pdb_file.geometry_restraints_manager()
    #  - params
    import mmtbx.command_line.rna_validate
    params = mmtbx.command_line.rna_validate.get_master_phil().extract()
  
    # run rna_validation
    self.result = rna_validation(
      pdb_hierarchy=pdb_hierarchy,
      geometry_restraints_manager=geometry,
      params=params,
      outliers_only=False)
  
  def add_file(self) :
    self.mdb_document["rna_validate"] = {}
    d = {"n_triaged":self.result.suites.n_triaged,
         "n_outliers":self.result.suites.n_outliers,
         "average_suiteness":self.result.suites.average_suiteness}
    self.mdb_document["rna_validate"]["suites"] = d
    d = {"n_total":self.result.puckers.n_total,
         "n_outliers":self.result.puckers.n_outliers}
    self.mdb_document["rna_validate"]["puckers"] = d
    d = {"n_total":self.result.angles.n_total,
         "n_outliers":self.result.angles.n_outliers}
    self.mdb_document["rna_validate"]["angles"] = d
    d = {"n_total":self.result.bonds.n_total,
         "n_outliers":self.result.bonds.n_outliers}
    self.mdb_document["rna_validate"]["bonds"] = d
  
  def add_residue(self) :
    pass
