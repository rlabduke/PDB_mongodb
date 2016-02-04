import os,sys
import argparse
import pdbe_rest
import mongodb_pdb_connect
import json

def get_current_pdbs() :
  pdblistfn = 'allpdbs.l'
  assert os.path.exists(pdblistfn)
  pdbs = {}
  n = 0
  fle = open(pdblistfn,'r')
  for l in fle :
    if not l.strip().endswith('.ent.gz') : continue
    pdb_id = l[l.rfind('/pdb')+4:l.rfind('.ent.gz')].upper()
    assert len(pdb_id) == 4
    mid = pdb_id[1:3]
    if mid not in pdbs.keys() : pdbs[mid] = []
    pdbs[mid].append(pdb_id)
    n += 1
  return n,pdbs

def get_json_pdb_ex_doc(pdb_id) :
  experiment = "/pdb/entry/experiment"
  pdbe_doc = pdbe_rest.get_request(experiment,pdb_id,True)
# pdbe_doc = pdbe_rest.get_request(summary,pdb_id,True)
  pdbe_json_doc = json.loads(pdbe_doc)
  return pdbe_json_doc

def get_missing(pdbs,pdbs_in_db) :
  missing = []
  for mid,lst in pdbs.items() :
    for pdb_id in lst :
      if mid not in pdbs_in_db.keys() or pdb_id not in pdbs_in_db[mid] :
         missing.append(pdb_id)
      #if pdb_id not in pdbs_in_db[mid] : missing.append(pdb_id)
      #if pdb_id in pdbs_in_db[mid] : print pdb_id,"in"
      #	else : print pdb_id,"out" 
  return missing

def print_json_pretty(d,log=sys.stdout) :
  print >> log, json.dumps(d,indent=4, separators=(',', ': '))

def run(args) :
  desc = "Updates the db called 'experiment' on the pdb_info db on daneel."
  parser = argparse.ArgumentParser(description=desc)
  n_pdbs,pdbs = get_current_pdbs()
  print >> sys.stderr, 'There are currently %i deposited pdbs.' % n_pdbs

  # connect to mongo on daneel
  c = 'experiment'
# c = 'file_info'
  mcon = mongodb_pdb_connect.Mongodb_PDB()
  n_pdbs_in_db,pdbs_in_db = mcon.get_existing_pdb_ids(collection=c)
  s = 'There are currently %i pdbs in the collection "%s".'
  print >> sys.stderr, s % (n_pdbs_in_db,c)

  # get missing
  missing = get_missing(pdbs,pdbs_in_db)
  s = 'There are currently %i pdbs missing from in the collection "%s".'
  print >> sys.stderr, s % (len(missing),c)


  for pdb_id in missing[:1] :
    print pdb_id
    doc = get_json_pdb_ex_doc(pdb_id)
    # I am assuming that the documents are a dict of length 1 and the key
    # is the pdb id with the value being a list of one dict (which is the one
    #  we're after). But I am unsure that this is universal thus the foloowing.
    assert len(doc) == 1, doc
    assert len(doc.keys()) == 1,len(doc.keys())
    k = doc.keys()[0]
    assert k.upper() == pdb_id
    assert type(doc[k]) is list,type(doc[k])
    assert len(doc[k]) == 1,len(doc[k])
    assert type(doc[k][0]) is dict,type(doc[k][0])
    print_json_pretty(doc[k][0])

if __name__ == '__main__' :
  run(sys.argv[1:])

