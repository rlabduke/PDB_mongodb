import os,sys
from mmtbx.command_line import fetch_pdb
from iotbx.file_reader import any_file
from iotbx import file_reader
import iotbx.pdb
import datetime
import json

validation_types = ['all','rna','clashscore']

class MDB_PDB_validation(object) :
 
  __slots__ = ['set_mdb_document','run_validation','add_file','add_residue'] 
  __slots__+= ['mdb_document','result']
  def __init__(self,pdb_file,validation_type) :
    self.pdb_file = pdb_file
    assert validation_type in validation_types, validation_type
    self.set_mdb_document()

  def set_mdb_document(self) :
    pdb_in = file_reader.any_file(self.pdb_file)
    pdb_in.check_file_type('pdb')
    self.mdb_document = get_pdb_meta_data(self.pdb_file)

  def write_pretty_mdb_document(self,log=sys.stdout) :
    sep = (',', ': ')
    s=json.dumps(self.mdb_document,sort_keys=True,indent=4,separators=sep)
    print >> log, s

def get_pdb_files(pdb_code,pdbcif=False) :
  # note that mtz and pdb cif can't be done simultaneously
  if pdbcif : args=["-c", pdb_code]
  else : args=["--mtz", pdb_code]
  data_files = fetch_pdb.run2(args)
  if pdbcif : return {"pdbcif":data_files}
  d = {}
  for df in  data_files :
    df_input_file = any_file(df)
    if df_input_file.file_type == "hkl" :
      if df.endswith('.mtz') : d["hklmtz"] = df
      elif df.endswith('.cif') : d["hklcif"] = df
    elif df_input_file.file_type == "pdb":
      d["pdb"] = df
    elif df_input_file.file_type == "seq":
      d["seq"] = df
  return d

def get_pdb_hierarchy(pdb_file) :
  pdb_in = file_reader.any_file(pdb_file)
  pdb_in.check_file_type('pdb')
  input_pdb = iotbx.pdb.input(file_name=pdb_file)
  return input_pdb.construct_hierarchy()

# mmCIF utils
def get_deposit_date(block) :
  deposit_date,release_date = None,None
  keys = ['_database_PDB_rev.num','_database_PDB_rev.date']
  keys+= ['_database_PDB_rev.date_original','_database_PDB_rev.status']
  keys+= ['_database_PDB_rev.replaces','_database_PDB_rev.mod_type']
  lists = {}
  if type(block["_database_PDB_rev.num"]) is str :
    dd = block["_database_PDB_rev.date_original"]
    rd = block["_database_PDB_rev.date"]
    return dd,rd
  lngth = -1
  for k in keys : 
    assert k in block.keys()
    lists[k] = [e for e in block[k]]
    if lngth == -1 : lngth = len(lists[k])
    assert len(lists[k]) == lngth, 'Unexpected length'
  # make a list of dictionaries where keys are _database_PDB_rev keys
  lofd = []
  for i in range(len(lists['_database_PDB_rev.num'])) :
    d = {}
    for k in keys :
      d[k] = lists[k][i]
    lofd.append(d)
    if d['_database_PDB_rev.mod_type'] == '0' :
      dd = d['_database_PDB_rev.date_original']
      deposit_date = datetime.datetime.strptime(dd, "%Y-%m-%d")
      rd = d['_database_PDB_rev.date']
      release_date = datetime.datetime.strptime(rd, "%Y-%m-%d")
  return dd,rd
  return deposit_date,release_date

def get_software(loop) :
  software = {}
  for row in loop.iterrows() :
    software[row['_software.classification']] = row['_software.name']
  return software

def get_pdb_meta_data(pdbcif_fn,pdb_code=None) :
  import iotbx.cif
  assert os.path.exists(pdbcif_fn)
  cif = iotbx.cif.reader(pdbcif_fn, strict=False).model()
  if not pdb_code :
    ks = cif.keys()
    assert len(ks) == 1,"trouble assigning pdb_code"
    pdb_code = ks[0]
  pdb_code = cif[pdb_code]["_entry.id"]
  pdb_meta_data = {"_id":pdb_code}
  pdb_meta_data['Experimental Method'] = cif[pdb_code]["_exptl.method"]
  pdb_meta_data['Resolution'] = cif[pdb_code]["_refine.ls_d_res_high"]
  # get dates
  deposit_date,release_date = get_deposit_date(cif[pdb_code])
  pdb_meta_data['Deposition Date'] = deposit_date
  pdb_meta_data['Release Date'] = release_date
  software = get_software(cif[pdb_code].get_loop("_software"))
  pdb_meta_data['software'] = software
  return pdb_meta_data
