from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pandas as pd

def where_NH(mol):
    # 找到所有的NH的位置
    patt = Chem.MolFromSmarts('[NH]')
    match_list = list(mol.GetSubstructMatches(patt))
    return(match_list)

def where_NH2(mol):
    # 找到所有的NH2的位置
    patt = Chem.MolFromSmarts('[NH2]')
    match_list = list(mol.GetSubstructMatches(patt))
    return(match_list)
  
def where_der_NH(mol):
    # 找到所有的能反应的NH的位置
    patt = Chem.MolFromSmarts('[#6;!$(C=*)]-[NH]-[#6;!$(C=*)]')
    match_list = list(mol.GetSubstructMatches(patt))
    return(match_list)

def where_der_NH2(mol):
    # 找到所有的能反应的NH2的位置
    patt = Chem.MolFromSmarts('[#6;!$(C=*)]-[NH2]')
    match_list = list(mol.GetSubstructMatches(patt))
    return(match_list)

def whether_NH12(mol):
    n1 = len(where_der_NH(mol))
    print("n1: ", n1)
    n2 = len(where_der_NH2(mol))
    print("n2: ", n2)
    print(n1 != 0)
    print(n2 != 0)
    if(n1 != 0 or n2 != 0):
        return(True)
    else:
        return(False)

def CheckTwoMatchList(match_list1, match_list):
    # 检查两个list的元组中, 后者的哪一个元组的元素属于前者中任意一个元组
    i = 0
    for tupleTmp in match_list:
        j = 0
        for tupleTmp1 in match_list1:
            if(tupleTmp[0] in tupleTmp1):
                return(i)
            j = j + 1
        i = i + 1
    return("")

def derNH12(mol):
    repl = Chem.MolFromSmiles('NS(=O)(=O)c1cccc2c(N(CC)CC)cccc12')
    match_list1 = where_der_NH(mol)
    match_list2 = where_der_NH2(mol)
    n1 = len(match_list1) # 可以反应的NH1的数量
    n2 = len(match_list2) # 可以反应的NH2的数量
    print("n1 = ", n1)
    print("n2 = ", n2)
    if(n1 == 0 and n2 == 0):
        # 既没有NH1也没有NH2(可反应的)
        return("")
    if(n1 != 0):
        # 有可以反应的NH1
        patt = Chem.MolFromSmarts('[NH]')
        match_list = where_NH(mol)
        n = len(match_list) #所有 NH1 的数量
        print("NH1 number:",n)
        if(n1 == n): # 可以反应的 NH1 等于所有的 NH1
            print("n1 == n")
            rms = AllChem.ReplaceSubstructs(mol,patt,repl, replaceAll = True)
            mol = rms[0]
            match_list1 = where_der_NH(mol)
            n1 = len(match_list1) # 可以反应的NH1的数量
            match_list = where_NH(mol)
            n = len(match_list) #所有 NH1 的数量
        elif(n1 < n): # 可以反应的 NH1 小于所有的 NH1
            print("n1 < n")
            while(n1 != 0):
                i = CheckTwoMatchList(match_list1, match_list)
                print("i = ",i)
                rms = AllChem.ReplaceSubstructs(mol,patt,repl, replaceAll = False)
                mol = rms[i]
                smiles = Chem.MolToSmiles(mol) # 必须做一步这个变化, 要不然mol的H会有问题
                mol = Chem.MolFromSmiles(smiles)
                #Draw.ShowMol(mol, size = (600,600), kekulize = True)
                match_list1 = where_der_NH(mol)
                n1 = len(match_list1) # 可以反应的NH1的数量
                match_list = where_NH(mol)
                n = len(match_list) #所有 NH1 的数量
    if(n2 != 0):
        # 有可以反应的NH2
        patt = Chem.MolFromSmarts('[NH2]')
        match_list = where_NH2(mol)
        n = len(match_list) # 所有的 NH2 的数量
        print("NH2 number:",n)
        if(n2 == n): # 可以反应的 NH2 等于所有的NH2
            print("n2 == n")
            rms = AllChem.ReplaceSubstructs(mol,patt,repl, replaceAll = True)
            mol = rms[0]
            match_list2 = where_der_NH(mol)
            n2 = len(match_list2) # 可以反应的NH2的数量
            match_list = where_NH2(mol)
            n = len(match_list) # 所有的 NH2 的数量
        elif(n2 < n): # 可以反应的NH2 数量小于所有的 NH2
            print("n2 < n")
            while(n2 != 0):
                i = CheckTwoMatchList(match_list2, match_list)
                if(i == ""):
                    break
                rms = AllChem.ReplaceSubstructs(mol,patt,repl, replaceAll = False)
                mol = rms[i]
                smiles = Chem.MolToSmiles(mol)
                mol = Chem.MolFromSmiles(smiles)
                match_list2 = where_der_NH2(mol)
                n2 = len(match_list2) # 可以反应的NH2的数量
                match_list = where_NH2(mol)
                n = len(match_list) #所有 NH2 的数量
    
    print("n1 = ", n1)
    print("n2 = ", n2)
    return(Chem.MolToSmiles(mol))

def whether_cOH(mol):
    patt = Chem.MolFromSmarts('c-[OH]')
    match_list = list(mol.GetSubstructMatches(patt))
    n = len(match_list)
    if(n != 0):
        return(True)
    else:
        return(False)
    
def where_cOH(mol):
    patt = Chem.MolFromSmarts('c-[OH]')
    match_list = list(mol.GetSubstructMatches(patt))
    return(match_list)

def dercOH(mol):
    patt = Chem.MolFromSmarts('c-[OH]')
    if(mol.HasSubstructMatch(patt)):
        repl = Chem.MolFromSmiles('COS(=O)(=O)c1cccc2c(N(CC)CC)cccc12')
        rms = AllChem.ReplaceSubstructs(mol,patt,repl, replaceAll = True)
        return(Chem.MolToSmiles(rms[0]))
    else:
        return('')

def Der_P(smiles):
  mol = Chem.MolFromSmiles(smiles)
  if(not whether_NH12(mol) and not whether_cOH(mol)):
      print(1)
      return("")
  if(whether_cOH(mol)):
      print(3)
      smiles = dercOH(mol)
      mol = Chem.MolFromSmiles(smiles)
  return(Chem.MolToSmiles(mol))

def Der_AP(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if(not whether_NH12(mol) and not whether_cOH(mol)):
        print(1)
        return("")
    if(whether_NH12(mol)):
        print(2)
        smiles = derNH12(mol)
        mol = Chem.MolFromSmiles(smiles)
    if(whether_cOH(mol)):
        print(3)
        smiles = dercOH(mol)
        mol = Chem.MolFromSmiles(smiles)
    return(Chem.MolToSmiles(mol))

def whether_COH(mol):
    patt = Chem.MolFromSmarts('[C;!$(C=O)]-[OH]')
    match_list = list(mol.GetSubstructMatches(patt))
    n = len(match_list)
    if(n != 0):
        return(True)
    else:
        return(False)
    
def where_COH(mol):
    patt = Chem.MolFromSmarts('[C;!$(C=O)]-[OH]')
    match_list = list(mol.GetSubstructMatches(patt))
    return(match_list)

def derCOH(mol):
    patt = Chem.MolFromSmarts('[C;!$(C=O)]-[OH]')
    if(mol.HasSubstructMatch(patt)):
        repl = Chem.MolFromSmiles('COS(=O)(=O)c1cccc2c(N(CC)CC)cccc12')
        rms = AllChem.ReplaceSubstructs(mol,patt,repl, replaceAll = False)
        return(Chem.MolToSmiles(rms[0]))
    else:
        return('')

def Der_Hy(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if(whether_COH(mol)):
        smiles = derCOH(mol)
        mol = Chem.MolFromSmiles(smiles)
    return(Chem.MolToSmiles(mol))
