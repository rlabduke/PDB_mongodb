import os,sys
import argparse
import pdbe_rest
import mongodb_pdb_connect
import json

def get_current_pdbs() :
  pdblistfn = 'allpdbs.l'
  assert os.path.exists(pdblistfn)
  pdbs = []
  fle = open(pdblistfn,'r')
  for l in fle :
    if not l.strip().endswith('.ent.gz') : continue
    pdbs.append(l[l.rfind('/pdb')+4:l.rfind('.ent.gz')].upper())
    assert len(pdbs[-1])
  return pdbs

def get_json_pdb_ex_doc(pdb_id) :
  experiment = "/pdb/entry/experiment"
  pdbe_doc = pdbe_rest.get_request(experiment,pdb_id,True)
# pdbe_doc = pdbe_rest.get_request(summary,pdb_id,True)
  pdbe_json_doc = json.loads(pdbe_doc)
  return pdbe_json_doc

def run(args) :
  desc = "Updates the db called 'experiment' on the pdb_info db on daneel."
  parser = argparse.ArgumentParser(description=desc)
#  parser.add_argument('pdb_file', help='A pdb file')
#  args = parser.parse_args()
#  assert os.path.exists(args.pdb_file)
  pdbs = get_current_pdbs()
  print >> sys.stderr, 'There are currently %i deposited pdbs.' % len(pdbs)

  # connect to mongo on daneel
  c = 'experiment'
# c = 'file_info'
  mcon = mongodb_pdb_connect.Mongodb_PDB()
  pdbs_in_db = mcon.get_existing_pdb_ids(collection=c)
  s = 'There are currently %i pdbs in the collection "%s".'
  print >> sys.stderr, s % (len(pdbs_in_db),c)

  # get missing
  missing = [ e for e in pdbs if e not in pdbs_in_db]
  s = 'There are currently %i pdbs missing from in the collection "%s".'
  print >> sys.stderr, s % (len(missing),c)


  for pdb_id in missing[:1] :
    print pdb_id
    print get_json_pdb_ex_doc(pdb_id)

if __name__ == '__main__' :
  run(sys.argv[1:])

