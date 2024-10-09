from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pandas as pd

def where_der_NH_(mol):
    # 找到所有的能反应的NH的位置
    patt = Chem.MolFromSmarts('[#6;!$(C=*)]-[NH]-[#6;!$(C=*)]')
    match_list = list(mol.GetSubstructMatches(patt))
    return(match_list)

def where_der_NH2_(mol):
    # 找到所有的能反应的NH2的位置
    patt = Chem.MolFromSmarts('[#6;!$(C=*)]-[NH2]')
    match_list = list(mol.GetSubstructMatches(patt))
    return(match_list)

def whether_NH12_(smiles):
    mol = Chem.MolFromSmiles(smiles)
    n1 = len(where_der_NH_(mol))
    n2 = len(where_der_NH2_(mol))
    if(n1 != 0 or n2 != 0):
        return(True)
    else:
        return(False)

def whether_cOH_(smiles):
    mol = Chem.MolFromSmiles(smiles)
    patt = Chem.MolFromSmarts('c-[OH]')
    match_list = list(mol.GetSubstructMatches(patt))
    n = len(match_list)
    if(n != 0):
        return(True)
    else:
        return(False)

def whether_COH_(smiles):
    mol = Chem.MolFromSmiles(smiles)
    patt = Chem.MolFromSmarts('[C;!$(C=O)]-[OH]')
    match_list = list(mol.GetSubstructMatches(patt))
    n = len(match_list)
    if(n != 0):
        return(True)
    else:
        return(False)
