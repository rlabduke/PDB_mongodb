import json
import sys
from elbow.formats.protein_sequence_setup import protein_sequence_to_three
from elbow.formats.dna_rna_sequence_parser import dna_rna_sequence_to_three
from elbow.formats.protein_sequence_setup import protein_sequence_to_names
from libtbx import group_args
import re

# protein
l = ["ACE","NME"]
reslist = [e for e in protein_sequence_to_three.values() if e not in l]
# result list is ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE',
#                 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN', 'PRO',
#                 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']

# DNA/RNA
reslist_na = [e for e in dna_rna_sequence_to_three.values() if e != None]
# reslist_na list is ['ADE', 'CYT', 'GUA', 'URI', 'THY']

# get bb and sc atoms
bb_atoms = ["N", "CA", "C", "O", "CB"]
sc_atoms = ["CB"]
for oneaa,atomlist in protein_sequence_to_names.items() :
  if oneaa in ['B','Z'] : continue
  for an in atomlist :
    if an not in bb_atoms and an not in sc_atoms : sc_atoms.append(an)

def print_json_pretty(d,log=sys.stdout) :
  print >> log, json.dumps(d,indent=4, separators=(',', ': '))

def get_resd(pdb_code,result) :
  return {'pdb_id'     : pdb_code,
          'model_id'   : None,
          'chain_id'   : result.chain_id,
          'icode'      : result.icode,
          'resseq'     : int(result.resseq),
          'altloc'     : result.altloc,
          'resname'    : result.resname}

class MDBAtom(object) :

  __slots__ = ['name','xyz','adp','rscc','twoFo_DFc_value']
  __slots__+= ['Fo_DFc_value','occ']

  def __init__(self,**kwargs) :
    for key, value in kwargs.iteritems():
      setattr(self, key, value)

  def aet_atom_dict(self) :
    d = {}
    for k in self.__slots__:
      try :
        d[k] = getattr(self,k)
      except AttributeError:
        d[k] = None
    return d

class MDBResidue(object) :

  __slots__ = ['pdb_id','model_id','chain_id','icode']
  __slots__+= ['resseq','altloc','resname','atoms','resolution','restype']
  __slots__+= ['rotalyze','ramalyze','omegalyze','cablam','clashes']

  def __init__(self,**kwargs) :
    for key, value in kwargs.iteritems():
      setattr(self, key, value)
      if hasattr(self,'pdb_id') : self.pdb_id = self.pdb_id.lower()
    self.atoms= []

  def deposit_atom(self,atom) :
    assert type(atom) == MDBAtom
    self.atoms.append(atom)

  def get_max(self,l) :
    return max(l, key=lambda d: d['value'])

  def get_min(self,l) :
    return min(l, key=lambda d: d['value'])

  def get_worst(self) :
    worst_bb  = {}
    worst_sc  = {}
    worst_all = {}
    # make a list of adp, rscc, 2fo-fc, and fo-fc for all, bb, and sc
    adps = []
    rsccs = []
    twoFo_DFc_values = []
    Fo_DFc_values = []
    adps_bb = []
    rsccs_bb = []
    twoFo_DFc_values_bb = []
    Fo_DFc_values_bb = []
    adps_sc = []
    rsccs_sc = []
    twoFo_DFc_values_sc = []
    Fo_DFc_values_sc = []
    for a in self.atoms :
      # all
      adps.append({'name':a.name,'value':a.adp})
      rsccs.append({'name':a.name,'value':a.rscc})
      twoFo_DFc_values.append({'name':a.name,'value':a.twoFo_DFc_value})
      Fo_DFc_values.append({'name':a.name,'value':a.Fo_DFc_value})
      # bb
      if a.name.strip() in bb_atoms :
        adps_bb.append({'name':a.name,'value':a.adp})
        rsccs_bb.append({'name':a.name,'value':a.rscc})
        twoFo_DFc_values_bb.append({'name':a.name,'value':a.twoFo_DFc_value})
        Fo_DFc_values_bb.append({'name':a.name,'value':a.Fo_DFc_value})
      # sc
      if a.name.strip() in sc_atoms :
        adps_sc.append({'name':a.name,'value':a.adp})
        rsccs_sc.append({'name':a.name,'value':a.rscc})
        twoFo_DFc_values_sc.append({'name':a.name,'value':a.twoFo_DFc_value})
        Fo_DFc_values_sc.append({'name':a.name,'value':a.Fo_DFc_value})
    # all
    if len(adps) > 0 : worst_all["adp"] = self.get_max(adps)
    if len(rsccs) > 0 : worst_all["rscc"] = self.get_min(rsccs)
    if len(twoFo_DFc_values) > 0 :
      worst_all["twoFo_DFc_value"] = self.get_min(twoFo_DFc_values)
    if len(Fo_DFc_values) > 0 :
      worst_all["Fo_DFc_value_max"] = self.get_max(Fo_DFc_values)
      worst_all["Fo_DFc_value_min"] = self.get_min(Fo_DFc_values)
    # bb
    if len(adps_bb) > 0 : worst_bb["adp"] = self.get_max(adps_bb)
    if len(rsccs_bb) > 0 : worst_bb["rscc"] = self.get_min(rsccs_bb)
    if len(twoFo_DFc_values_bb) > 0 :
      worst_bb["twoFo_DFc_value"] = self.get_min(twoFo_DFc_values_bb)
    if len(Fo_DFc_values_bb) > 0 :
      worst_bb["Fo_DFc_value_max"] = self.get_max(Fo_DFc_values_bb)
      worst_bb["Fo_DFc_value_min"] = self.get_min(Fo_DFc_values_bb)
    # sc
    if len(adps_sc) > 0 : worst_sc["adp"] = self.get_max(adps_sc)
    if len(rsccs_sc) > 0 : worst_sc["rscc"] = self.get_min(rsccs_sc)
    if len(twoFo_DFc_values_sc) > 0 :
      worst_sc["twoFo_DFc_value"] = self.get_min(twoFo_DFc_values_sc)
    if len(Fo_DFc_values_sc) > 0 :
      worst_sc["Fo_DFc_value_max"] = self.get_max(Fo_DFc_values_sc)
      worst_sc["Fo_DFc_value_min"] = self.get_min(Fo_DFc_values_sc)

    return worst_bb,worst_sc,worst_all

  def add_rotalyze_result(self,result) :
    #print dir(result);exit()
    chi_angles = result.chi_angles
    while None in chi_angles : chi_angles.remove(None)
    self.rotalyze = group_args(
        is_outlier          = result.is_outlier(),
        rotalyze_evaluation = result.evaluation,
        chi_angles          = chi_angles,
        rotamer             = result.rotamer_name,
        score               = result.score)

  def add_ramalyze_result(self,result) :
    #print dir(result);exit()
    self.ramalyze = group_args(
        is_outlier    = result.is_outlier(),
        type          = result.ramalyze_type(),
        phi           = result.phi,
        psi           = result.psi,
        score         = result.score)

  def add_omegalyze_result(self,result) :
    #print dir(result);exit()
    self.omegalyze = group_args(
        is_outlier    = result.is_outlier(),
        omega         = result.omega,
        type          = result.omegalyze_type())

  def add_cablam_result(self,result) :
    if result.nextres :
      resd = get_resd(self.pdb_id,result.nextres)
      MDBRes = MDBResidue(**resd)
      nextkey = MDBRes.get_residue_key()
    else : nextkey = None
    if result.prevres :
      resd = get_resd(self.pdb_id,result.prevres)
      MDBRes = MDBResidue(**resd)
      prevkey = MDBRes.get_residue_key()
    else : prevkey = None
    #print result.nextres.id_str()
    #print result.prevres.id_str()
    #print dir(result.nextres)
    #print dir(result);exit()
    self.cablam = group_args(
         next                 = nextkey,
         prev                 = prevkey,
         mu_in                = result.measures.mu_in,
         mu_out               = result.measures.mu_out,
         nu                   = result.measures.nu,
         ca_virtual           = result.measures.ca_virtual,
         cablam               = result.scores.cablam,
         c_alpha_geom         = result.scores.c_alpha_geom,
         alpha_score          = result.scores.alpha,
         beta_score           = result.scores.beta,
         threeten_score       = result.scores.threeten,
         alpha_feedback       = result.feedback.alpha,
         beta_feedback        = result.feedback.beta,
         threeten_feedback    = result.feedback.threeten,
         structure_suggestion = result.find_single_structure_suggestion(),
         cablam_outlier       = result.feedback.cablam_outlier,
         cablam_disfavored    = result.feedback.cablam_disfavored,
         c_alpha_geom_outlier = result.feedback.c_alpha_geom_outlier)

  def add_clash(self,resd) :
    if not hasattr(self,'clashes') : self.clashes = []
    self.clashes.append(resd)

  def is_protein(self) :
    # Bases soley on resname
    return self.resname.upper() in reslist

  def is_nucleic_acid(self) :
    # Bases soley on resname
    return self.resname.upper() in reslist_na

  def get_module_dict(self,module) :
    # module should point to a group_args object. This puts the contents
    # thereof into a dict.
    assert module in ['rotalyze','ramalyze','omegalyze','cablam']
    regexob = re.compile('^__.*__')
    nd = {}
    modob = getattr(self,module)
    for k in modob.__dict__ :
      if regexob.match(k) : continue
      nd[k] = getattr(modob,k)
    return nd

  def get_residue_mongodoc(self) :
#   al = [a.aet_atom_dict() for a in self.atoms]
    al = {}
    for a in self.atoms :
      ad = a.aet_atom_dict()
      atomname = ad['name'].strip()
      ad.pop('name')
      al[atomname] = ad
    if self.is_protein : restyp = "protein"
    elif self.is_nucleic_acid : restyp = "na"
    else : retyp = None
    d = {'_id':self.get_residue_key(),
         'pdb_id':self.pdb_id,
         'model_id':self.model_id,
         'chain_id':self.chain_id,
         'icode':self.icode,
         'resseq':self.resseq,
         'altloc':self.altloc,
         'resname':self.resname,
         'restype':restyp}
    if hasattr(self,'resolution') : d['resolution'] = self.resolution
    if len(al) > 0 :
      d['atoms']  = al
      if self.resname.upper() in reslist+reslist_na :
        d['worst_bb'],d['worst_sc'],d['worst_all'] = self.get_worst()
    regexob = re.compile('^__.*__')
    if hasattr(self,'rotalyze') :
      d['rotalyze'] = self.get_module_dict('rotalyze')
    if hasattr(self,'ramalyze') :
      d['ramalyze'] = self.get_module_dict('ramalyze')
    if hasattr(self,'omegalyze') :
      d['omegalyze'] = self.get_module_dict('omegalyze')
    if hasattr(self,'cablam') :
      d['cablam'] = self.get_module_dict('cablam')
    if hasattr(self,'clashes') :
      d['clashes'] = self.clashes
    return d

  def get_residue_key(self) :
    # this is suppose to be unique r each residue and
    # will be used as the _id in mongo
   key = ''
   l = ['pdb_id','model_id','chain_id','icode']
   l+= ['resseq','altloc','resname','atoms','resolution']
   for attr in l :
     if attr in ['atoms','resolution'] : continue
     s = getattr(self,attr)
     if s is None : continue 
     if type(s) == str :   key += s.strip()
     elif type(s) == int : key += '%i' % s
   return key
