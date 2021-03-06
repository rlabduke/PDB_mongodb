import os,sys
from mmtbx.command_line import load_model_and_data
from mmtbx.command_line import generate_master_phil_with_inputs
from phenix.automation.refinement import refinement_base, refinement_callback
from iotbx.file_reader import any_file
from iotbx import file_reader
import iotbx.pdb
import datetime
import json
import mdb_utils
import re

# here is all the validation types. to add a validation type, add the name
# of the type to this list, add a function to run the analysis below (adding
# the info into the object as appropriate), and then add the option to run it
# in run_mongo_validation.py.
validation_types = ['all','rna','clashscore','rotalyze','ramalyze']
validation_types+= ['omegalyze','cablam','rscc','geometry']

class MDB_PDB_validation(object) :
 
  __slots__ = ['pdb_file','hklmtz_file', 'set_mdb_document','run_validation']
  __slots__+= ['add_file','add_residue','mdb_document','result','mdb_document']
  __slots__+= ['residues','meta_data','detail','hierarchy','pdb_code']
  __slots__+= ['do_flips','high_resolution','cmdline']
  def __init__(self,pdb_file,hklmtz_file,
               detail,high_resolution=None,mdb_document=None,pdb_code=None,
               do_flips=False) :
    assert detail in ['file','residue'],detail
    assert type(do_flips) == bool
    self.pdb_file = pdb_file
    self.hklmtz_file = hklmtz_file
    self.detail = detail
    self.pdb_code = pdb_code
    self.high_resolution = high_resolution
    self.do_flips = do_flips
    if not pdb_code : self.pdb_code = 'N/A'
    pdb_in = file_reader.any_file(pdb_file)
    self.hierarchy = pdb_in.file_object.hierarchy
    args = [self.pdb_file]
    if self.hklmtz_file : args.append(self.hklmtz_file)
    self.cmdline = load_model_and_data(
      args=args,
      master_phil=generate_master_phil_with_inputs(""),
      require_data=False,
      create_fmodel=True,
      process_pdb_file=True,
      prefer_anomalous=True)
    # keys are res ids and values are MDBResidue objects.
    if self.detail == 'residue' :
      self.initiate_residues()
    self.set_mdb_document(mdb_document)

  # This initiates the residues dictionary with 'all' the residues in the given
  # pdb. Because we do this step, subsequent reidue-level criteria will already
  # have a MDBResidue to wich the criteria can be added to.
  def initiate_residues(self) :
    print >> sys.stderr, 'initializing residues...\n'
    self.residues = {}
    for chain in self.cmdline.pdb_hierarchy.chains():
      for residue_group in chain.residue_groups():
        for conformer in residue_group.conformers():
          for residue in conformer.residues():
            resd = {'pdb_id'     : self.pdb_code,
                    'model_id'   : None,
                    'chain_id'   : chain.id,
                    'icode'      : residue.icode,
                    'resseq'     : residue.resseq_as_int(),
                    'altloc'     : conformer.altloc,
                    'resname'    : residue.resname}
            MDBRes = mdb_utils.MDBResidue(**resd)
            reskey = MDBRes.get_residue_key()
            #print >> sys.stdout, reskey
            self.residues[reskey] = MDBRes
    #reskeys = self.residues.keys()
    #reskeys.sort()
    #for k in reskeys :
    #  print k
    #  if k.find('A79') != -1 : print k
    #exit()

  def set_mdb_document(self,mdb_document) :
    if mdb_document is not None : self.mdb_document = mdb_document;return
    pdb_in = file_reader.any_file(self.pdb_file)
    pdb_in.check_file_type('pdb')
    try : self.meta_data = get_pdb_meta_data(self.pdb_file)
    except : 
     self.meta_data = None
     return
    # if the detail of the output is 'file' then of course the meta data is
    # relevent to the output.
    if self.detail == 'file' : self.mdb_document = self.meta_data.copy()
    else : self.mdb_document = None

  def write_pretty_mdb_document(self,log=sys.stdout) :
    sep = (',', ': ')
    s=json.dumps(self.mdb_document,sort_keys=True,indent=4,separators=sep)
    print >> log, s

  def write_pretty_residue_mdb_documents(self,log=sys.stdout) :
    sep = (',', ': ')
    residue_keys = self.residues.keys()
    residue_keys.sort()
    for rkey in residue_keys :
      res = self.residues[rkey]
      s=json.dumps(res.get_residue_mongodoc(),
                   sort_keys=True,indent=4,separators=sep)
      print >> log, s

  def get_alternate_keys(self,resd) :
    if resd['model_id'] == None : resd['model_id'] = ''
    reskey = [resd['pdb_id'],resd['model_id'],resd['chain_id']]
    reskey+= [resd['icode'].strip(),'%i'%resd['resseq'],'[a-zA-Z]']
    reskey+= [resd['resname']]
    rk = ':'.join(reskey)
    mo = re.compile(rk)
    reskeys = []
    for k in self.residues.keys() :
      if mo.match(k) : reskeys.append(k)
    return reskeys

  def get_res_key_and_atom(self, atom_info, alt) :
    resd  = mdb_utils.get_resd(self.pdb_code,atom_info)
    if alt : resd['altloc'] = alt[resd['resseq']]
    MDBRes = mdb_utils.MDBResidue(**resd)
    reskey = MDBRes.get_residue_key()
    #assert reskey in self.residues.keys(), reskey+'  '+atom_info.name
    if reskey not in self.residues.keys() : # the residue has alts
      reskeys = self.get_alternate_keys(resd)
    else : reskeys = [reskey]
    return reskeys, atom_info.name 

  def run_clashscore_validation(self) :
    from val_clashscore import CLASHSCOREvalidation
    try :
      vc = CLASHSCOREvalidation(self.pdb_file,self.detail,self.mdb_document,
             self.pdb_code,self.do_flips)
    except : pass
    else :
      self.mdb_document = vc.mdb_document
      if self.detail == 'residue' :
        for clash in vc.result.results :
          clashatoms = []
          for i,atom in enumerate(clash.atoms_info) :
            resd = mdb_utils.get_resd(self.pdb_code,atom)
            MDBRes = mdb_utils.MDBResidue(**resd)
            reskey = MDBRes.get_residue_key()
            # When dealing with alts in only a portion of the residue, reskey 
            # (as just defined) will not be found in self.residues.keys().
            # The reason is that initiate_residues takes a residue with n alts
            # and splits into n separate residues. e.g. if a residue has alt
            # A and B in one residue atom for just sidchain atoms, there will
            # exist two residues, A and B, both having the atoms that don't have
            # alts and each having their respective alt atoms. When this is the
            # case we remedy by putting the clash in both A and B. This is what
            # the folllowing code is doing.
            if reskey not in self.residues.keys():
              reskeys = self.get_alternate_keys(resd)
              reskey = [k for k in reskeys]
            clashatoms.append({'targ_reskey':reskey,
                               'targname':atom.name,
                               'overlap':clash.overlap})
          assert len(clashatoms) == 2
          # targ and src can be confusing here. The following puts targ reskey
          # targname in the coresponding src residue. Take the following clash :
          #    A  72  ARG  HG2  A  72  ARG  O   :-1.038
          # clashatoms[0] = {'targname': ' HG2', 'targ_reskey': '1ubqA72ARG'}
          # clashatoms[1] = {'targname': ' O  ', 'targ_reskey': '1ubqA72ARG'}
          # we then put the src name in next :
          clashatoms[0]['srcname'] = clashatoms[1]['targname']
          clashatoms[1]['srcname'] = clashatoms[0]['targname']
          # next we add clashatoms[1] to the residue object corresponding to
          # clashatoms[0] and vise versa.
          if type(clashatoms[0]['targ_reskey']) == list :
            for tk in clashatoms[0]['targ_reskey'] :
              if type(clashatoms[1]['targ_reskey']) == list :
                for sk in clashatoms[1]['targ_reskey'] :
                  newdict = {'targ_reskey':sk,
                             'targname':clashatoms[1]['targname'],
                             'overlap':clashatoms[1]['overlap']}
                  self.residues[tk].add_clash(newdict)
              else :
                self.residues[tk].add_clash(clashatoms[1])
          else:
            self.residues[clashatoms[0]['targ_reskey']].add_clash(clashatoms[1])
          if type(clashatoms[1]['targ_reskey']) == list :
            for tk in clashatoms[1]['targ_reskey'] :
              if type(clashatoms[0]['targ_reskey']) == list :
                for sk in clashatoms[0]['targ_reskey'] :
                  newdict = {'targ_reskey':sk,
                             'targname':clashatoms[0]['targname'],
                             'overlap':clashatoms[0]['overlap']}
                  self.residues[tk].add_clash(newdict)
              else :
                self.residues[tk].add_clash(clashatoms[0])
          else:
            self.residues[clashatoms[1]['targ_reskey']].add_clash(clashatoms[0])

  def run_rna_validation(self) :
    from val_rna import RNAvalidation
    vc = RNAvalidation(self.pdb_file,self.detail,self.mdb_document)
    self.mdb_document = vc.mdb_document

  def get_alt_if_any(self,atoms_info) :
    alt = None
    for atom_info in atoms_info :
      resd  = mdb_utils.get_resd(self.pdb_code,atom_info)
      resseq,altloc = resd['resseq'],resd['altloc']
      if alt is None : alt = {resseq:altloc}
      # This function is used for seeing if there are alternates in a list of 
      # atoms involved in an angle or length measure so if there is an altloc
      # they should be the same if they have the same resseq.
      else :
        if resseq in alt.keys():
          if alt[resseq] != altloc and alt[resseq] == '' : alt[resseq] = altloc
          if altloc != '': assert alt[resseq]==altloc, altloc + '   ' + str(alt)
        else : alt[resd['resseq']] = resd['altloc']
    return alt

  def run_geometry_validation(self) :
    from mmtbx.validation import restraints
    restraints = restraints.combined(
      pdb_hierarchy=self.cmdline.pdb_hierarchy,
      xray_structure=self.cmdline.xray_structure,
      geometry_restraints_manager=self.cmdline.geometry,
      ignore_hd=True,
      cdl=False)
    restraints.angles.show()
    # get bond legth outlier data
    for result in restraints.bonds.results :
      #print >> sys.stderr, dir(result)
      # Since cases can exist where one of the atoms have an alt and the other
      # does not, wee need to check if any of the participating atoms have
      # alts.
      alt = self.get_alt_if_any(result.atoms_info)
      # Since this is a bond length there are two atoms involved 
      reskeys0,name0 = self.get_res_key_and_atom(result.atoms_info[0],alt)
      reskeys1,name1 = self.get_res_key_and_atom(result.atoms_info[1],alt)
      for reskey0 in reskeys0 :
        for reskey1 in reskeys1 :
          self.residues[reskey0].add_bondlength_result(result,reskey0,name0,
                                                              reskey1,name1)
      if reskey0 != reskey1 :
        self.residues[reskey1].add_bondlength_result(result,reskey1,name1,
                                                            reskey0,name0)
    for result in restraints.angles.results :
      #print >> sys.stderr, dir(result)
      alt = self.get_alt_if_any(result.atoms_info)
      # Since this is a angle there are three atoms involved 
      reskeys0,name0 = self.get_res_key_and_atom(result.atoms_info[0],alt)
      reskeys1,name1 = self.get_res_key_and_atom(result.atoms_info[1],alt)
      reskeys2,name2 = self.get_res_key_and_atom(result.atoms_info[2],alt)
      # The logic that follows inserts the angle outlier info into the
      # two or one residue[s] involved. At least two of the three residues
      # will be the same residue.
      for reskey0 in reskeys0 :
        for reskey1 in reskeys1 :
          for reskey2 in reskeys2 :
            self.residues[reskey0].add_angle_result(result,reskey0,name0,
                                                           reskey1,name1,
                                                           reskey2,name2)
            if reskey0 != reskey1 :
              self.residues[reskey1].add_angle_result(result,reskey1,name1,
                                                             reskey0,name0,
                                                             reskey2,name2)
            elif reskey0 != reskey2 :
              self.residues[reskey2].add_angle_result(result,reskey2,name2,
                                                             reskey1,name1,
                                                             reskey0,name0)


  def run_rotalyze(self) :
    from mmtbx.validation import rotalyze
    rotalyze_result = rotalyze.rotalyze(self.hierarchy)
    for result in rotalyze_result.results : 
      resd = mdb_utils.get_resd(self.pdb_code,result)
      MDBRes = mdb_utils.MDBResidue(**resd)
      reskey = MDBRes.get_residue_key()
      if reskey not in self.residues.keys(): # alternates likely exist
        reskeys = self.get_alternate_keys(resd)
        for k in reskeys :
          self.residues[k].add_rotalyze_result(result)
      else : # No alternates
        self.residues[reskey].add_rotalyze_result(result)

  def run_ramalyze(self) :
    from mmtbx.validation import ramalyze
    ramalyze_result = ramalyze.ramalyze(self.hierarchy)
    for result in ramalyze_result.results :
      resd = mdb_utils.get_resd(self.pdb_code,result)
      MDBRes = mdb_utils.MDBResidue(**resd)
      reskey = MDBRes.get_residue_key()
      if reskey not in self.residues.keys(): # alternates likely exist
        reskeys = self.get_alternate_keys(resd)
        for k in reskeys :
          self.residues[k].add_ramalyze_result(result)
      else : # No alternates
        self.residues[reskey].add_ramalyze_result(result)

  def run_omegalyze(self) :
    from mmtbx.validation import omegalyze
    omegalyze_result = omegalyze.omegalyze(
                               pdb_hierarchy = self.hierarchy,
                               nontrans_only = False,
                               out           = sys.stdout,
                               quiet         = False)
    for result in omegalyze_result.results :
      resd = mdb_utils.get_resd(self.pdb_code,result)
      MDBRes = mdb_utils.MDBResidue(**resd)
      reskey = MDBRes.get_residue_key()
      if reskey not in self.residues.keys(): # alternates likely exist
        reskeys = self.get_alternate_keys(resd)
        for k in reskeys :
          self.residues[k].add_omegalyze_result(result)
      else : # No alternates
        self.residues[reskey].add_omegalyze_result(result)

  def run_cablam(self) :
    from mmtbx.validation import cablam
    cablam_result = cablam.cablamalyze(
                               pdb_hierarchy = self.hierarchy,
                               outliers_only = False,
                               out           = sys.stdout,
                               quiet         = False)
    t = True
    for result in cablam_result.results :
      resd = mdb_utils.get_resd(self.pdb_code,result)
      MDBRes = mdb_utils.MDBResidue(**resd)
      reskey = MDBRes.get_residue_key()
      if not result.prevres : continue
      if reskey in self.residues.keys() :
        self.residues[reskey].add_cablam_result(result)
      elif MDBRes.altloc != '' :
        # alts exist in reskey but the actual residue has no alt
        #assert MDBRes.altloc != '','"%s"' % MDBRes.altloc
        newkey = MDBRes.get_residue_key(no_alt=True)
        if newkey not in self.residues.keys() :
          print >> sys.stderr, 'WARNING : trouble finding %s' % reskey
        else :
          #print '  ' + newkey,newkey in self.residues.keys()
          cal = result.altloc
          result.altloc = ' '
          self.residues[newkey].add_cablam_result(result,cablam_altloc=cal)
      else : # Side chain has alternates but bb does not.
        for l in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' :
          newkey = MDBRes.get_residue_key(assign_alt=l)
          if newkey in self.residues.keys() :
            result.altloc = l
            self.residues[newkey].add_cablam_result(result)
          newkey = MDBRes.get_residue_key(assign_alt=l.lower())
          if newkey in self.residues.keys() :
            result.altloc = l.lower()
            result.altloc = self.residues[newkey].add_cablam_result(result)

  def run_rscc(self) :
    if not self.hklmtz_file :
      print >> sys.stderr, '*'*79 + '\nSkipping rscc - no hkl file found'
      return
    import val_rscc
    rscob = val_rscc.RSCCvalidation(pdb_file        = self.pdb_file,
                                    hklmtz_file     = self.hklmtz_file,
                                    pdb_code        = self.pdb_code,
                                    high_resolution = self.high_resolution)
    for reskey,atomlist in rscob.mdb_residues.items() :
      for atomd in atomlist :
        self.residues[reskey].deposit_atom(mdb_utils.MDBAtom(**atomd))

class ModelData(refinement_base) :

  def __init__(self, args=None, params=None, out=None,
               tmp_dir=None, dest_dir=None, extra_args=()) :
    super(ModelData, self).__init__(
      args=args, params=params, out=out, tmp_dir=tmp_dir, dest_dir=dest_dir)

  def run(self):
    self.get_inputs()
    self.extract_model()
    self.extract_data()

def get_pdb_files(pdb_code,pdbcif=False) :
  from mmtbx.command_line import fetch_pdb
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
    pdb_meta_data['Resolution'] = float(cif[pdb_code]["_refine.ls_d_res_high"])
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
