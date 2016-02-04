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
  if pdbe_doc is None : return None
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
  desc = """Updates the pdb_info.experiment on daneel. The script has the following steps :

  - get a list of currently deposited pdbs (in allpdbs.l) 
  - get a list of pdbs in pdb_info.experiment. 
  - get a list of deposited pdbs missing from pdb_info.experiment
  - iterate trough the missing pdbs :
    - get the given pdb's experiment documenti[s] from the PDBe
    - insert that/those document into pdb_info.experiment 

"""
  parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-v','--verbose',help='Verbose outpur',
                      action='store_true')
  args = parser.parse_args()


  # get a list of currently deposited pdbs (in allpdbs.l) 
  n_pdbs,pdbs = get_current_pdbs()
  print >> sys.stderr, 'There are currently %i deposited pdbs.' % n_pdbs


  # connect to mongo on daneel
  c = 'experiment'
  mcon = mongodb_pdb_connect.Mongodb_PDB()


  # get a list of pdbs in pdb_info.experiment
  n_pdbs_in_db,pdbs_in_db = mcon.get_existing_pdb_ids(collection=c)
  s = 'There are currently %i pdbs in the collection "%s".'
  print >> sys.stderr, s % (n_pdbs_in_db,c)


  # get a list of deposited pdbs missing from pdb_info.experiment
  missing = get_missing(pdbs,pdbs_in_db)
  s = 'There are currently %i pdbs missing from in the collection "%s".'
  print >> sys.stderr, s % (len(missing),c)


  # iterate trough the missing pdbs
  if len(missing) > 2000 : factor = 1000
  else : factor = 100
  factor = 100
  print >> sys.stderr, '\n\nBegin iterating missing pdbs...\n'
  msg = '\n%i records inserted - %.2f %% done.\n'
  for i,pdb_id in enumerate(missing) :
    if args.verbose : print >> sys.stderr, "working on %s..." % pdb_id
    # get the given pdb's experiment document from the PDBe
    doc = get_json_pdb_ex_doc(pdb_id)
    if doc is None :
      s = '\nWARNING: Skipping %s -- PDBe request returned None\n' % pdb_id
      print >> sys.stderr, s
      continue
    # I am assuming that the documents are a dict of length 1 and the key
    # is the pdb id with the value being a list of one or more dict (which 
    # is/are the one[s] we're after). But I am unsure that this is universal 
    # thus the foloowing.
    #print_json_pretty(doc)
    assert len(doc) == 1, doc
    assert len(doc.keys()) == 1,len(doc.keys())
    k = doc.keys()[0]
    assert k.upper() == pdb_id
    assert type(doc[k]) is list,type(doc[k])
    assert type(doc[k][0]) is dict,type(doc[k][0])
    #print_json_pretty(doc[k][0])

    # insert that/those document into pdb_info.experiment
    for j,mdoc in enumerate(doc[k]) :
      mdoc["_id"] = {"pdb_id":pdb_id,"n":j}
      mcon.upsert_document(mdoc,c)
      if args.verbose :
        print >> sys.stderr, '  %s inserted into pdb_info.%s' % (pdb_id,c)
      if i%factor == 0 : print >> sys.stderr, msg % (i,(i*100.0)/len(missing))
  print >> sys.stderr, '\n\n%i total records inserted\n\n' % i

if __name__ == '__main__' :
  run(sys.argv[1:])

