import os,sys
from mmtbx.command_line import fetch_pdb
from iotbx.file_reader import any_file
from iotbx import file_reader
import iotbx.pdb
import datetime
import json

validation_types = ['all','rna','clashscore']

class MDB_PDB_validation(object) :
 
  __slots__ = ['pdb_file', 'set_mdb_document','run_validation','add_file'] 
  __slots__+= ['add_residue','mdb_document','result','mdb_document']
  def __init__(self,pdb_file,mdb_document=None) :
    self.pdb_file = pdb_file
    self.set_mdb_document(mdb_document)

  def set_mdb_document(self,mdb_document) :
    if mdb_document is not None : self.mdb_document = mdb_document;return
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

def get_pdb_summary(pdb_file) :
  hrcy = get_pdb_hierarchy(pdb_file)
# for e in dir(hrcy) : print e
  d = {}
  for k,v in hrcy.overall_counts().resname_classes.items() :
    d[k] = v
  # add to d as needed
  d["contains_nucleic_acid"] = hrcy.contains_nucleic_acid()
  d["contains_rna"] = hrcy.contains_rna()
  d["contains_protein"] = hrcy.contains_protein()
  if hrcy.contains_nucleic_acid() and not hrcy.contains_rna() :
    d["contains_dna"] = True
  return d

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

def get_software(block) :
  software = {}
  loop = block.get_loop("_software")
  if loop is not None:
    for row in loop.iterrows() :
      software[row['_software.classification']] = row['_software.name']
  else :
    if '_software.classification' not in block.keys() or '_software.name' :
      return None
    software[block['_software.classification']] = block['_software.name']
  return software

def get_computing(block) :
  keys = [ k for k in block.keys() if k.startswith("_computing")]
  computing = {}
  for k in keys :
    keyl = k.split(".")
    m = "k should be of the format '_xxxx.yyyy'\n k = %s"
    assert len(keyl) == 2, m % k 
    if block[k] == "?" : continue
    computing[keyl[1]] = block[k]
  if len(computing) : return
  return computing

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
  # get counts
  pdb_meta_data['summary'] = get_pdb_summary(pdbcif_fn)
  # get summary
  pdb_meta_data['Experimental Method'] = cif[pdb_code]["_exptl.method"]
  if "_refine.ls_d_res_high" in cif[pdb_code].keys() :
    pdb_meta_data['Resolution'] = cif[pdb_code]["_refine.ls_d_res_high"]
  # get dates
  deposit_date,release_date = get_deposit_date(cif[pdb_code])
  pdb_meta_data['Deposition Date'] = deposit_date
  pdb_meta_data['Release Date'] = release_date
  # get software
  software = get_software(cif[pdb_code])
  if software : pdb_meta_data['software'] = software
  # get computing
  computing = get_computing(cif[pdb_code])
  if computing : pdb_meta_data['computing'] = computing
  return pdb_meta_data
