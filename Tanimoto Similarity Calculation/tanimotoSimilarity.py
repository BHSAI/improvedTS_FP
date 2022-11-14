from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_structure_pipeline import standardizer

def standarize_chembl(smiles_list):
  count = len(smiles_list)
  smiles_std_list = [None] * count
  for index in range(0,count):
    smiles = smiles_list[index]
    if smiles is not None:
      molecule = Chem.MolFromSmiles(smiles)
      if molecule is not None:
        m_no_salt = standardizer.get_parent_mol(molecule)
        molecule_std = standardizer.standardize_mol(m_no_salt[0])
        if molecule_std is not None:
          smiles_std = Chem.MolToSmiles(molecule_std)
          smiles_std_list[index] = smiles_std
  return smiles_std_list


def calculateFP(smiles_list, radius):
  count = len(smiles_list)
  fingerprint = [None] * count
  for index in range(0,count):
    smiles = smiles_list[index]
    if smiles is not None:
      molecule = Chem.MolFromSmiles(smiles)
      if molecule is not None:
        fp = AllChem.GetMorganFingerprint(molecule,radius)
        dict_element = fp.GetNonzeroElements()
        fingerprint[index] = set(dict_element.keys())
  return fingerprint


def calculateTDMatrix(fpSet1, fpSet2):
  len_set1 = len(fpSet1)
  len_set2 = len(fpSet2)
  td_matrix = [None] * len_set1
  for set1Index in range(0,len_set1):
    fp1 = fpSet1[set1Index]
    td_row = [None] * len_set2
    for set2Index in range(0,len_set2):
      fp2 = fpSet2[set2Index]
      td = calculateTD(fp1, fp2)
      td_row[set2Index] = td
    td_matrix[set1Index] = td_row
  return td_matrix

def calculateTD(fpSet1, fpSet2):
  a = len(fpSet1)
  b = len(fpSet2)
  #the common
  intersection = fpSet1.intersection(fpSet2)
  c = len(intersection)
  td = c / ((a + b) - c)
  print(fpSet1)
  print(fpSet2)
  print(intersection)
  print(str(a) + ' ' + str(b) + ' ' + str(c) + ' '+ str(td))
  return td


#example
querySMILES = ["C1CCC1","CC(=O)CC"]
referenceSMILES = ["CC(=O)OC1=CC=CC=C1C(=O)O","OC(c1cnccc1)=O"]

querySMILES_std = standarize_chembl(querySMILES)
referenceSMILES_std = standarize_chembl(referenceSMILES)

radius = 2
querySMILES_fp = calculateFP(querySMILES_std, radius)
referenceSMILES_fp = calculateFP(referenceSMILES_std, radius)

tdMatrix = calculateTDMatrix(querySMILES_fp, referenceSMILES_fp)
tdMatrix
