import os,sys
import argparse
from utils import pdb_utils

log = sys.stderr

def get_args() :
  desc = "This script runs various validation programs on a given PDB"
  desc+= " or mmcif file and returns a JSON document(s) having the various"
  desc+= " validation criteria for the given pdb."
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_code', help='A pdb code')
  parser.add_argument('--pdb_cif', dest='cif_fn', help='A pdb mmcif')
  parser.add_argument('-d','--detail', dest='detail', help='file or residue')
  parser.add_argument('--write_out_file',
    dest='write_out_file', help='write file', action='store_const', const=True)
  vth = 'validation_type can be one of the following:\n%s  ' 
  parser.add_argument("--validation_type", dest='validation_type',
    help=vth % ', '.join(pdb_utils.validation_types))
  hs = 'A directory where the output document will be written'
  parser.add_argument('--outdir', dest='outdir', help=hs)
  hs = 'don\'t cleanup downloaded files'
  parser.add_argument('--dont-cleanup', dest='dont_cleanup', help=hs,
                      action='store_const', const=True)
  args = parser.parse_args()


  if not args.validation_type : args.validation_type = 'all'
  if not args.detail : args.detail = 'file'
  assert args.detail in ['file','residue'],args.detail
  return args

def run (out=sys.stdout, quiet=False) :
  args = get_args()

  if not args.cif_fn :
    assert len(args.pdb_code) == 4
    pdb_files = pdb_utils.get_pdb_files(args.pdb_code,pdbcif=True)
    pdbfn = pdb_files['pdbcif']
    print >> log, '%s downloaded' % pdbfn
  elif args.cif_fn :
    pdbfn = args.cif_fn
    args.dont_cleanup = True
  else : RuntimeError("Must provide a pdb code at the least")
  setattr(args,'pdb_file_path',pdbfn)

  if not args.outdir : outdir = os.getcwd()
  else : outdir = args.outdir
  assert os.path.exists(outdir)

  # MDB_PDB_validation will get meta data from the pdb_file. If the detail
  # is file, it will initiate the mdb_document (an attribute of the class)
  # with the meta data. Meta datfa here refers to resoluton, deposition date,
  # Experimental Method, and summary info on molecule contents (aa, hoh, na,
  # protein, and rna)
  validation_class = pdb_utils.MDB_PDB_validation(pdb_file = args.pdb_file_path,
                                                  detail   = args.detail)
  meta_data = validation_class.meta_data
  if args.validation_type in ['rna','all'] :
    if meta_data['summary']['contains_rna'] :
      validation_class.run_rna_validation()

  if args.validation_type in ['clashscore','all'] :
    validation_class.run_clashscore_validation()



  if args.write_out_file :
    if args.outdir :
      bd = os.path.join(outdir,args.pdb_code[1:3])
      if not os.path.exists(bd) : os.makedirs(bd)
    else : bd = outdir
    fn = os.path.join(bd,'%s.validate' % args.pdb_code)
    fle = open(fn,'w')
    validation_class.write_pretty_mdb_document(log=fle)
    fle.close()
    print >> log, '%s written' % fn
  else : validation_class.write_pretty_mdb_document()
  if os.path.exists(args.pdb_file_path) and not args.dont_cleanup :
    for k,fn in pdb_files.items() : os.remove(fn)
    print >> log, 'Cleaned up.'

if (__name__ == "__main__") :
  run()
