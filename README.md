# PDB_mongodb

This repository contains code meant to run on mmcif PDB files (not PDB format). The base class for all runs is pdb_utils.MDB_PDB_validation. This class has a method called 'write_pretty_mdb_document' which will print out a jsom type document. The strait dictionary is pdb_utils.MDB_PDB_validation.mdb_document.

To run, source python :

  $ python run_mongo_validation.py -h
  $ python run_mongo_validation.py 4fe5 --validation_type clashscore
  $ python run_mongo_validation.py 4fe5 --validation_type rna 

