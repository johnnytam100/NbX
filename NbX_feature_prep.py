# Import library
from pyrosetta import *
import os
import glob
import pandas as pd
import re

from pyrosetta.rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.protocols.antibody import *
from pyrosetta.rosetta.protocols.antibody.design import *
from pyrosetta.rosetta.utility import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

from biopandas.pdb import PandasPdb

from argparse import ArgumentParser

# Init
init()

parser = ArgumentParser()
parser.add_argument('--antigen_chain', type=str)
parser.add_argument('--antibody_chain', type=str)
parser.add_argument('--native', type=str)
args = parser.parse_args()

### Define path to python
path_to_python = "/home/cltam/anaconda3/envs/nbx/bin/python"

### Define path to DockQ
path_to_dockq = "/data/cltam/script/DockQ/"

### Define correct pose .pdb name
if args.native != None:
    correct_pose_pdb = args.native

### Define path to FoldX
path_to_foldx = "/data/cltam/script/FoldX/foldx_20221231" 

### Define antigen chain and antibody chain
antigen_chain = args.antigen_chain
antibody_chain = args.antibody_chain

### Define desired epitope to be targeted (i.e. positive epitope)
epitope_positive_start_residue = 9998
epitope_positive_end_residue = 9999

### Define epitope to be avoided (i.e. negative epitope)
epitope_negative_start_residue = 9998
epitope_negative_end_residue = 9999

### Define CDR start residue and end residue

CDR1_start_residue = 24
CDR1_end_residue = 42

CDR2_start_residue = 57
CDR2_end_residue = 69

CDR3_start_residue = 107
CDR3_end_residue = 138


### Define a function that create a list with chain and resdiue numbers
def create_chain_res_list(chain, start_res, end_res):
  empty_list = []
  for i in range(start_res, end_res+1):
    empty_list.append(chain+str(i))
  return empty_list


### Define a function that extract records from a dataframe according to match in a certain column
def extract_record(df, col_name, match):
  return df.loc[df[col_name] == match]


### Define a function that multiplies residue count with aadescriptors and sum the product for 20 aa 
def multiply_residue_count_with_aadescriptor_sum(residue_count_df, aa_descriptor_df, chain):

  result = pd.DataFrame(0, index = range(1), columns=aa_descriptor_df.columns)

  for i in range(residue_count_df.shape[1]):
    multiply_df = pd.DataFrame(residue_count_df.iloc[:,i:i+1].values * aa_descriptor_df.iloc[i:i+1,:].values)
    result = pd.DataFrame(result.values +  multiply_df.values)

  new_col_name = []
  for elem2 in aa_descriptor_df.iloc[i:i+1,:].columns:
    new_col_name.append("chain" + chain + "_" + elem2)

  result.columns = new_col_name
  
  return result


### Create lists containing every residues of every CDRs for counting
CDR1 = create_chain_res_list(antibody_chain,CDR1_start_residue,CDR1_end_residue)
CDR2 = create_chain_res_list(antibody_chain,CDR2_start_residue,CDR2_end_residue)
CDR3 = create_chain_res_list(antibody_chain,CDR3_start_residue,CDR3_end_residue)
epitope_positive = create_chain_res_list(antigen_chain,epitope_positive_start_residue,epitope_positive_end_residue)
epitope_negative = create_chain_res_list(antigen_chain,epitope_negative_start_residue,epitope_negative_end_residue)

### Create lists containing wild card for antigen and antibody to count all interacting residues
antigen_chain_wildcard = re.compile("." + antigen_chain + ".+")
antibody_chain_wildcard = re.compile("." + antibody_chain + ".+")

### Define residue list for counting individual residue
residue_list = ['G','A','V','P','I','L','F','W','S','T','C','E','D','Q','N','H','K','R','Y','M']
residue_list.sort()


### Create empty dataframe to collect all features from all .pdb
all_pdb_summary_df = pd.DataFrame()


### Collect all .pdb paths in the current directory into a list, sort the list
filepath_list = []

for filepath in glob.iglob('./*.pdb'):
  filepath_list.append(filepath)

filepath_list.sort()

### Loop through all *.pdb in current directory
for filepath in filepath_list:

  # try:

      ### Grep ATOM
      os.system(str("grep 'ATOM' " + filepath + " > " + filepath[2:-4] + "_grep-ATOM.pdb"))

      ### Calculate DockQ score

      if args.native != None:

          os.system(str(path_to_python + ' ' + path_to_dockq + "DockQ.py" + ' ' + filepath[2:-4] + "_grep-ATOM.pdb" + ' ' + correct_pose_pdb + ' > ' + filepath[2:-4] + "_grep-ATOM.pdb" + '.DockQ')) 

          filein = open(filepath[2:-4] + "_grep-ATOM.pdb" + '.DockQ', 'r', encoding='utf-8')
          readline = filein.readlines()

          dockq = []

          for i in readline[-5:]:
            dockq.append(i.replace('\n',''))

          dockq_df = pd.DataFrame([sub.split(" ") for sub in dockq])
          dockq_df = dockq_df.set_index(dockq_df.columns[0])
          dockq_df = dockq_df.iloc[:,0:1].T

          dockq_df.reset_index(drop=True, inplace=True)


      # Import a pose to PyRosetta
      pose = pose_from_pdb(filepath[2:-4] + "_grep-ATOM.pdb")

      # Run PyRosetta InterfaceAnalyzerMover
      iam = InterfaceAnalyzerMover()
      iam.set_pack_separated(True)
      iam.apply(pose)

      score = pose.scores
      score_df = pd.DataFrame.from_dict(score, orient='index')
      score_df_transposed = score_df.T
      score_df_transposed = score_df_transposed.iloc[:,:20]

      ### Run FoldX AnalyzeComplex
      os.system(str(path_to_foldx + " --command=AnalyseComplex --pdb=" + filepath[2:-4] + "_grep-ATOM.pdb" + " --analyseComplexChains=A,H"))
      

      ### Collect results from FoldX AnalyzeComplex
      filein = open('Indiv_energies_' + filepath[2:-4] + "_grep-ATOM" + '_AC.fxout', 'r', encoding='utf-8')
      readline = filein.readlines()

      energy1 = readline[-1].replace('\n','')
      energy1_split = energy1.split('\t')
      energy2 = readline[-2].replace('\n','')
      energy2_split = energy2.split('\t')

      energy1_df = pd.DataFrame([energy1_split])
      energy2_df = pd.DataFrame([energy2_split])

      group1 = energy1_split[1]
      group2 = energy2_split[1]

      header = readline[-3].replace('\n','')
      header_split = header.split('\t')
      header_split_concat1 = ['chain' + group1 + ' - ' + elem for elem in header_split]
      header_split_concat2 = ['chain' + group2 + ' - ' + elem for elem in header_split]

      energy1_df.columns = header_split_concat1
      energy2_df.columns = header_split_concat2

      FoldX_AnalyzeComplex_energy_df = pd.concat([energy2_df.iloc[:,1:], energy1_df.iloc[:,1:]], axis = 1)
      FoldX_AnalyzeComplex_energy_df.pop('chainH - Group')
      FoldX_AnalyzeComplex_energy_df.pop('chainA - Group')


      ### Read interface residue output from FoldX AnalyzeComplex
      filein = open('Interface_Residues_' + filepath[2:-4] + "_grep-ATOM" + '_AC.fxout', 'r', encoding='utf-8')
      readline = filein.readlines()
      interface_residue = readline[-1].replace('\n','').split('\t')

      ### Count interacting residues in 1) CDR1, 2) CDR2, 3) CDR3, 4) positive epitope, 5) negative epitope and 6) paratope

      interface_residue_reduced = [i[1:] for i in interface_residue]

      CDR1_count = len(set(interface_residue_reduced).intersection(CDR1))
      CDR2_count = len(set(interface_residue_reduced).intersection(CDR2))
      CDR3_count = len(set(interface_residue_reduced).intersection(CDR3))

      epitope_positive_count = len(set(interface_residue_reduced).intersection(epitope_positive))
      epitope_negative_count = len(set(interface_residue_reduced).intersection(epitope_negative))

      antigen_count = len(list(filter(antigen_chain_wildcard.match,interface_residue)))
      antibody_count = len(list(filter(antibody_chain_wildcard.match,interface_residue)))

      ### convert .pdb to dataframe by biopandas
      ppdb = PandasPdb()
      ppdb.read_pdb(filepath[2:-4] + "_grep-ATOM.pdb")
      pdb_df = extract_record(ppdb.df['ATOM'], 'chain_id', antibody_chain)

      ### count the number of residues
      residue_set = set()

      for i in pdb_df['residue_number']:
        residue_set.add(i)
      
      ### count full length of CDR1, CDR2 and CDR3
      CDR1_full_list = []
      CDR2_full_list = []
      CDR3_full_list = []

      for i in residue_set:
        if i >= CDR1_start_residue and i <= CDR1_end_residue:
          CDR1_full_list.append(i)
          CDR1_full_number = len(CDR1_full_list)

      for i in residue_set:
        if i >= CDR2_start_residue and i <= CDR2_end_residue:
          CDR2_full_list.append(i)
          CDR2_full_number = len(CDR2_full_list)

      for i in residue_set:
        if i >= CDR3_start_residue and i <= CDR3_end_residue:
          CDR3_full_list.append(i)
          CDR3_full_number = len(CDR3_full_list)


      ### Create a dataframe to collect all the epitope and CDR count results
      epitope_CDR_count_df = pd.DataFrame(
          {
              
          'epitope_positive_count' : epitope_positive_count,
          'epitope_negative_count' : epitope_negative_count,

          'CDR1_count': CDR1_count,
          'CDR2_count': CDR2_count,
          'CDR3_count': CDR3_count,
          
          'antigen_total_count' : antigen_count,
          'antibody_total_count' : antibody_count,
          
          'CDR1_FL' : int(CDR1_full_number),
          'CDR2_FL' : int(CDR2_full_number),
          'CDR3_FL' : int(CDR3_full_number)
          
          
          }, index = [0])
      
      epitope_CDR_count_df.reset_index(drop=True, inplace=True)


      ### Calculate ratio over total counts
      epitope_CDR_count_df['epitope_positive_count/antigen_total'] = epitope_CDR_count_df['epitope_positive_count']/epitope_CDR_count_df['antigen_total_count']
      epitope_CDR_count_df['epitope_negative_count/antigen_total'] = epitope_CDR_count_df['epitope_negative_count']/epitope_CDR_count_df['antigen_total_count']
      epitope_CDR_count_df['CDRs_count/antibody_total'] = (epitope_CDR_count_df['CDR1_count'] + epitope_CDR_count_df['CDR2_count'] + epitope_CDR_count_df['CDR3_count']) /epitope_CDR_count_df['antibody_total_count']

      epitope_CDR_count_df['CDR1_count/CDR1_FL'] = epitope_CDR_count_df['CDR1_count']/epitope_CDR_count_df['CDR1_FL']
      epitope_CDR_count_df['CDR2_count/CDR2_FL'] = epitope_CDR_count_df['CDR2_count']/epitope_CDR_count_df['CDR2_FL']
      epitope_CDR_count_df['CDR3_count/CDR3_FL'] = epitope_CDR_count_df['CDR3_count']/epitope_CDR_count_df['CDR3_FL']


      antigen_individual_residue_count_df = pd.DataFrame()

      for i in residue_list:
        r = re.compile(i + antigen_chain + ".+")
        residue_count = len(list(filter(r.match,interface_residue)))
        col_name = 'chain' + antigen_chain + '_' + i

        antigen_individual_residue_count_df[col_name] = [residue_count]


      antibody_individual_residue_count_df = pd.DataFrame()

      for i in residue_list:
        r = re.compile(i + antibody_chain + ".+")
        residue_count = len(list(filter(r.match,interface_residue)))
        col_name = 'chain' + antibody_chain + '_' + i

        antibody_individual_residue_count_df[col_name] = [residue_count]


      ### Multiply the individual interface residue counts with aadescriptors
      aa_descriptor = pd.read_csv('aaDescriptors.csv', index_col=0)
      antigen_aadescriptor_df = multiply_residue_count_with_aadescriptor_sum(antigen_individual_residue_count_df, aa_descriptor, antigen_chain)
      antibody_aadescriptor_df = multiply_residue_count_with_aadescriptor_sum(antibody_individual_residue_count_df, aa_descriptor, antibody_chain)


      if args.native != None:

          each_pdb_summary_df = pd.concat([dockq_df,
                                          epitope_CDR_count_df, 
                                          antigen_individual_residue_count_df,
                                          antibody_individual_residue_count_df,
                                          score_df_transposed, 
                                          FoldX_AnalyzeComplex_energy_df,
                                          antigen_aadescriptor_df,
                                          antibody_aadescriptor_df], 
                                          
                                          axis = 1)

      else:

          each_pdb_summary_df = pd.concat([epitope_CDR_count_df, 
                                    antigen_individual_residue_count_df,
                                    antibody_individual_residue_count_df,
                                    score_df_transposed, 
                                    FoldX_AnalyzeComplex_energy_df,
                                    antigen_aadescriptor_df,
                                    antibody_aadescriptor_df], 
                                    
                                    axis = 1)
      
      each_pdb_summary_df.insert(0, 'PDB', filepath[2:-4])


      # Append results to the empty dataframe
      all_pdb_summary_df = all_pdb_summary_df.append(each_pdb_summary_df)

      # Delete files from FoldX AnalyzeComplex
      os.system(str("rm -rf *" + filepath[2:-4] + "_grep-ATOM*")) 

  # except:

      # os.system(str("rm -rf *" + filepath[2:-4] + "_grep-ATOM*")) 


### Save the summary feature dataframe to a .csv
all_pdb_summary_df.to_csv('NbX_feature.csv', index=False)