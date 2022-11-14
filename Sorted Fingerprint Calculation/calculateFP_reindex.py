import pickle
import pandas as pd
from collections import Counter
from rdkit import Chem
from rdkit.Chem import AllChem
import os

def calculateFP_reindex(smilesList, weights, nlen, radius):
  FP_result = [None] * len(smilesList)
  for row in range(0,len(smilesList)):
    smiles = smilesList[row]
    if smiles is not None:
      molecule = Chem.MolFromSmiles(smiles)
      if molecule is not None:
        fp = AllChem.GetMorganFingerprint(molecule,radius)
        dict_element = fp.GetNonzeroElements()
        counts = [0] * nlen
        FP_result[row] = counts
        for idx, c_ in weights:
          count=dict_element.get(idx)
          if count is not None:
            newidx=idx%nlen
            counts[newidx] += dict_element[idx]
  return FP_result

smiles = ["C1CCC1","CC(=O)CC"]
#the file savedcounts can be found under the data folder
#it is compiled by extracting features from 400k compounds
countfilePath = "data/savedcounts"
weights = pickle.load( open( countfilePath, "rb" ) ).most_common()
nlen = 2048
radius = 2
fp = calculateFP_reindex(smiles, weights, nlen, radius)
