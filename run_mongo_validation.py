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
  parser.add_argument('--pdb_file', dest='pdb_fn', help='A pdb file')
  parser.add_argument('--hkl_mtz',dest='hklmtz_fn',help='An hkl mtz')
  parser.add_argument('-d','--detail', dest='detail', help='file or residue')
  parser.add_argument('--write_out_file',
    dest='write_out_file', help='write file', action='store_const', const=True)
  vth = 'validation_type can be one of the following:\n%s  ' 
  parser.add_argument("--validation_type", nargs='+',
    help=vth % ', '.join(pdb_utils.validation_types))
  hs = 'A directory where the output document will be written'
  parser.add_argument('--outdir', dest='outdir', help=hs)
  hs = 'don\'t cleanup downloaded files'
  parser.add_argument('--dont-cleanup', dest='dont_cleanup', help=hs,
                      action='store_const', const=True)
  parser.add_argument('--do_flips', help='do flips when calculating clashes',
                      action='store_const', const=True)
  args = parser.parse_args()

  if not args.validation_type : args.validation_type = 'all'
  if not args.detail : args.detail = 'file'
  assert args.detail in ['file','residue'],args.detail
  if args.detail == 'file' : args.validation_type = 'clashscore'
  return args

def run (out=sys.stdout, quiet=False) :
  args = get_args()

  get_meta = True
  pdb_files = {}
  if not args.cif_fn and not args.pdb_fn :
    assert len(args.pdb_code) == 4
    pdb_files = pdb_utils.get_pdb_files(args.pdb_code,pdbcif=True)
    pdbfn = pdb_files['pdbcif']
    print >> log, '%s downloaded' % pdbfn
  elif args.cif_fn :
    pdbfn = args.cif_fn
    args.dont_cleanup = True
  elif args.pdb_fn :
    assert args.pdb_code
    pdbfn = args.pdb_fn
    args.dont_cleanup = True
    get_meta = False
  else : RuntimeError("Must provide a pdb code at the least")
  setattr(args,'pdb_file_path',pdbfn)

  #get and print meta data
  high_resolution = None
  if get_meta :
    meta_data = pdb_utils.get_pdb_meta_data(args.pdb_file_path)
    print >> log, '*'*79 + '\nSummary:'
    mdb_utils.print_json_pretty(meta_data,log)
    print >> log, '*'*79
    if "Resolution" in meta_data.keys() :
      high_resolution = meta_data["Resolution"]
  else : print >> log, '*'*79 + '\nNo Summary available.\n' + '*'*79
  if get_meta : 
    is_xray = meta_data["Experimental Method"] == "X-RAY DIFFRACTION"
  elif args.hklmtz_fn :
    is_xray = True
  allrscc = 'all' in args.validation_type or 'rscc' in args.validation_type
  if is_xray and allrscc :
    if not args.hklmtz_fn :
      modeldata = pdb_utils.ModelData(args=[args.pdb_code])
      if 'pdbcif' not in pdb_files.keys() :
        pdb_files = pdb_utils.get_pdb_files(args.pdb_code,pdbcif=True)
      pdb_files['hklmtz'] = modeldata.mtz_file
      hklmtzfn = pdb_files['hklmtz']
      print >> log, '%s downloaded' % hklmtzfn
    elif args.hklmtz_fn :
      hklmtzfn = args.hklmtz_fn
      args.dont_cleanup = True
    setattr(args,'hklmtz_file_path',hklmtzfn)
  else : setattr(args,'hklmtz_file_path',None)

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
  if args.do_flips is None : args.do_flips = False
  validation_class = pdb_utils.MDB_PDB_validation(
                                            pdb_file    = args.pdb_file_path,
                                            hklmtz_file = args.hklmtz_file_path,
                                            detail      = args.detail,
                                            pdb_code    = args.pdb_code,
                                            high_resolution = high_resolution,
                                            do_flips        = args.do_flips)
  getRNA = False
  if type(args.validation_type) == str :
    args.validation_type = [args.validation_type] 
  for validation_type in args.validation_type :
    if validation_type in ['rna','all'] :
      if get_meta and meta_data['summary']['contains_rna'] : metaRNA = True
      else : metaRNA = False
      if validation_type == 'rna' or metaRNA : getRNA = True
      if getRNA :
        print >> log, 'Running rna validation...\n'
        validation_class.run_rna_validation()

    if validation_type in ['clashscore','all'] :
      print >> log, 'Running clashscore...\n'
      validation_class.run_clashscore_validation()

    if validation_type in ['rotalyze','all'] :
      print >> log, 'Running rotalyze...\n'
      validation_class.run_rotalyze()

    if validation_type in ['ramalyze','all'] :
      print >> log, 'Running ramalyze...\n'
      validation_class.run_ramalyze()

    if validation_type in ['omegalyze','all'] :
      print >> log, 'Running omegalyze...\n'
      validation_class.run_omegalyze()

    if validation_type in ['cablam','all'] :
      print >> log, 'Running cablam...\n'
      validation_class.run_cablam()

    if validation_type in ['rscc','all'] :
      print >> log, 'Running real-space correlation...\n'
      validation_class.run_rscc()

    if validation_type in ['geometry','all'] :
      print >> log, 'Running geometry validation...\n'
      validation_class.run_geometry_validation()

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
