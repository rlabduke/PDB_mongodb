import json
import sys

def print_json_pretty(d,log=sys.stdout) :
  print >> log, json.dumps(d,indent=4, separators=(',', ': '))

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
  __slots__+= ['resseq','altloc','resname','atoms']

  def __init__(self,**kwargs) :
    for key, value in kwargs.iteritems():
      setattr(self, key, value)
    self.atoms= []

  def deposit_atom(self,atom) :
    assert type(atom) == MDBAtom
    self.atoms.append(atom)

  def get_residue_mongodoc(self) :
    al = [a.aet_atom_dict() for a in self.atoms]
    d = {'_id':self.get_residue_key(),
         'pdb_id':self.pdb_id,
         'model_id':self.model_id,
         'chain_id':self.chain_id,
         'icode':self.icode,
         'resseq':self.resseq,
         'altloc':self.altloc,
         'resname':self.resname,
         'atoms':al}
    return d

  def get_residue_key(self) :
    # this is suppose to be unique r each residue and
    # will be used as the _id in mongo
   key = ''
   for attr in self.__slots__ :
     if attr == 'atoms' : continue
     s = getattr(self,attr)
     if s is None : continue 
     key += s.strip()
   return key
