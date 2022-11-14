This page provides computer codes and datasets described in the manuscript Methods to improve the performance of topological fingerprints in machine learning and molecular similarity calculations 

## Whats included:
1. Two datasets used for evaluation
2. Manuscript that describes the theory in detail
3. Tanimoto Similarity Calculation (Set implementation).
4. Sorted Fingerprint Calculation

The latter two offer in both python and KNIME implementations. 

## Tanimoto Similarity Calculation (Set implementation)
The program takes in two lists of SMILES and do the followings in order: 
1. Standarize
2. Calculate fingerprint with RDKit fingerprint
   
   The fingerprint is a set of numbers which represent the indexes of features. The index of a feature is determined by the RDKit.
   
3. Calculate the tanimoto distance by comparing two sets instead of two arrays.

## Sorted Fingerprint Calculation
The program takes in a list of SMILES and do the followings in order:
1. Standarize
2. Calculate fingerprint with RDKit fingerprint
3. Sort the fingerprint based on a weight file
   
   The weight file was calculated by running RDKit Morgan fingerprint on a larget set of compounds and reocrding the occurance of features.
   
4. Fold the fingerprint into the desired length

## Using Python Scripts
Dependencies
1. Python Verssion 3.0+
2. Rdkit module [Link](https://www.rdkit.org/docs/Install.html)
   
   conda create -c conda-forge -n my-rdkit-env rdkit
   
3. CheMBL module [Link](https://github.com/chembl/ChEMBL_Structure_Pipeline)
   
   conda install -c chembl chembl_structure_pipeline
   
The python code are goruped into functions with sample code to run in the end.

## Using Knime Workflows
1. install KNIME [Link](https://www.knime.com/installation)
2. setup python integration [Link](https://docs.knime.com/2018-12/python_installation_guide/index.html#_introduction)
   and the dependencies mentioned in the previous section
3. import workflows
4. enter data into starting nodes
5. execute a selected workflow

The KNIME editor allows mixing code of different languages with ease. The implementation uses both python and java. 
The Java code uses Hashset when it comes to the Tanimoto Distance calculation. 
The KNIME edtior is version 4.5+ and requires java 8.
