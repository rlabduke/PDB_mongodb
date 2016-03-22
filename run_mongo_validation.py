import os,sys
import argparse
from utils import pdb_utils
from utils import mdb_utils

log = sys.stderr

def get_args() :
  desc = "This script runs various validation programs on a given PDB"
  desc+= " or mmcif file and returns a JSON document(s) having the various"
  desc+= " validation criteria for the given pdb."
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_code', help='A pdb code')
  parser.add_argument('--pdb_cif', dest='cif_fn', help='A pdb mmcif')
  parser.add_argument('--hkl_mtz',dest='hklmtz_fn',help='An hkl mtz')
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
  if not args.hklmtz_fn :
    pdb_files = pdb_utils.get_pdb_files(args.pdb_code)
    pdb_files['pdbcif'] = args.pdb_file_path
    hklmtzfn = pdb_files['hklmtz']
    print >> log, '%s downloaded' % hklmtzfn
  elif args.hklmtz_fn :
    hklmtzfn = args.hklmtz_fn
    args.dont_cleanup = True
  setattr(args,'hklmtz_file_path',hklmtzfn)

  if not args.outdir : outdir = os.getcwd()
  else : outdir = args.outdir
  assert os.path.exists(outdir)

  # MDB_PDB_validation will get meta data from the pdb_file, the attribute is
  # called meta_data. If the detail is file, it will initiate the mdb_document
  # (an attribute of the class) with meta_data. Meta data here refers to
  # resoluton, deposition date, experimental Method, and summary info 
  # on molecule contents (aa, hoh, na, protein, and rna)
  # If detail is residue then MDB_PDB_validation will have a residues object 
  # which is a dict with keys being a res id and values being MDBResidue
  # objects.
  validation_class = pdb_utils.MDB_PDB_validation(pdb_file = args.pdb_file_path,
                                            hklmtz_file = args.hklmtz_file_path,
                                                 detail   = args.detail,
                                                 pdb_code = args.pdb_code)
  meta_data = validation_class.meta_data
  print >> log, '*'*79 + '\nSummary:'
  mdb_utils.print_json_pretty(meta_data,log)
  print >> log, '*'*79
  if args.validation_type in ['rna','all'] :
    print >> log, 'Running rna validation...\n'
    if meta_data['summary']['contains_rna'] :
      validation_class.run_rna_validation()

  if args.validation_type in ['clashscore','all'] :
    print >> log, 'Running clashscore...\n'
    validation_class.run_clashscore_validation()

  if args.validation_type in ['rotalyze','all'] :
    print >> log, 'Running rotalyze...\n'
    validation_class.run_rotalyze()

  if args.validation_type in ['ramalyze','all'] :
    print >> log, 'Running ramalyze...\n'
    validation_class.run_ramalyze()

  if args.validation_type in ['omegalyze','all'] :
    print >> log, 'Running omegalyze...\n'
    validation_class.run_omegalyze()

  if args.validation_type in ['cablam','all'] :
    print >> log, 'Running cablam...\n'
    validation_class.run_cablam()

  if args.validation_type in ['rscc','all'] :
    print >> log, 'Running real-space correlation...\n'
    validation_class.run_rscc()


  if args.write_out_file :
    if args.outdir :
      bd = os.path.join(outdir,args.pdb_code[1:3])
      if not os.path.exists(bd) : os.makedirs(bd)
    else : bd = outdir
    fn = os.path.join(bd,'%s.validate' % args.pdb_code)
    fle = open(fn,'w')
  else : fle = sys.stdout
  if args.detail == 'file': 
    validation_class.write_pretty_mdb_document(log=fle)
  else :
    validation_class.write_pretty_residue_mdb_documents(log=fle)
  if args.write_out_file :
    fle.close()
    print >> log, '%s written' % fn

  if os.path.exists(args.pdb_file_path) and not args.dont_cleanup :
    for k,fn in pdb_files.items() : os.remove(fn)
    print >> log, 'Cleaned up.'

if (__name__ == "__main__") :
  run()
