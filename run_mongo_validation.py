import os,sys
import argparse
from utils import pdb_utils

log = sys.stderr

def get_args() :
  desc = "Run RNA validate on a given PDB code or a given mmCIF file"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_code', help='A pdb code')
  parser.add_argument('--pdb_cif', dest='cif_fn', help='A pdb mmcif')
  parser.add_argument('--detail', dest='detail', help='file or residue')
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

  validation_class = pdb_utils.MDB_PDB_validation(args.pdb_file_path)
  if args.validation_type in ['rna','all'] :
    doc = validation_class.mdb_document
    if doc['summary']['contains_rna'] :
      from val_rna import RNAvalidation
      validation_class = RNAvalidation(args.pdb_file_path,
                                       mdb_document=doc)
      if args.detail == 'file' : validation_class.add_file()
      else : validation_class.add_residue()

  if args.validation_type in ['clashscore','all'] :
    doc = validation_class.mdb_document
    from val_clashscore import CLASHSCOREvalidation
    validation_class = CLASHSCOREvalidation(args.pdb_file_path,
                                       mdb_document=doc)
    if args.detail == 'file' : validation_class.add_file()
    else : validation_class.add_residue()

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
