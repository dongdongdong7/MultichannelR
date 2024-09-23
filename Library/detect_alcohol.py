from rdkit import Chem
'''
smiles = "C(O)CC[C@H](O)C(O)=O"
smiles = "N[C@@H](CC1=CC=C(O)C(I)=C1)C(O)=O"
patt = Chem.MolFromSmarts("[C;!$(C=O)][OX2H]")
mol = Chem.MolFromSmiles(smiles)
match = mol.GetSubstructMatches(patt)
len = len(match)
'''
def whether_COH(smiles):
  patt = Chem.MolFromSmarts("[C;!$(C=O)][OX2H]")
  mol = Chem.MolFromSmiles(smiles)
  match = mol.GetSubstructMatches(patt)
  n = len(match)
  if(n!=0):
    return(True)
  else:
    return(False)
