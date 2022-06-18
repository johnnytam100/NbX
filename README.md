# NbX

---

This repository contains the NbX models for the re-ranking of nanobody–antigen binding poses.

Tam, C.; Kumar, A.; Zhang, K.Y.J. NbX: Machine Learning-Guided Re-Ranking of Nanobody–Antigen Binding Poses. Pharmaceuticals 2021, 14, 968. https://doi.org/10.3390/ph14100968

---

# How to run NbX

## Step 1 : install pre-requisites

PyRosetta (https://www.pyrosetta.org/downloads)

FoldX (https://foldxsuite.crg.eu/)

DockQ (optional, https://github.com/bjornwallner/DockQ)

Rosetta (optional, https://new.rosettacommons.org/demos/latest/tutorials/install_build/install_build)

`pip install biopandas==0.4.1`

`pip install xgboost==0.90`

`pip install scikit-learn==0.22.2.post1`

`pip install joblib==1.1.0`

`pip install dill==0.3.5.1`

## Step 2 : change paths

Change the following paths inside `NbX_feature_prep.py`

`path_to_python = "/home/cltam/anaconda3/envs/nbx/bin/python"`

`path_to_foldx = "/data/cltam/script/FoldX/foldx_20221231"`

`path_to_dockq = "/data/cltam/script/DockQ/" (optional)`


## Step 3 : renumber nanobody
Before any docking, please renumber your nanobody with PyIgClassify (http://dunbrack2.fccc.edu/pyigclassify/)

If your Nb (or Nb-Ag) structure is confidential and you don't want to submit to a webserver: 

modify the CDRs start and end residue numbers (search `"CDR1_start_residue"`) inside `NbX_feature_prep.py`.


## Step 4 : change directory, copy

`cd run_NbX`

`cp -r ../model ../NbX_feature_prep.py ../aaDescriptors.csv ../NbX_predict.py ./`


## Step 5 : feature prep

### (option 1,  without DockQ)

`python NbX_feature_prep.py --antigen_chain A --antibody_chain H`

### (option 2,  with DockQ)

`python NbX_feature_prep.py --antigen_chain A --antibody_chain H --native 6oq8_complex.pdb`


## Step 6 : predict

`python NbX_predict.py`


## Step 7 : analyse results

**Descendingly sort the `mean_predicted_CAPRI_binary_proba` in `NbX_prediction.csv`**, we get the following results

![投影片1](https://user-images.githubusercontent.com/51283097/174423865-865a8b73-d382-4080-b080-8fa49e5b2a44.PNG)

---

# RosettaDock (optional)

To mimic the NbX benchmark setting, you can perform RosettaDock refinement of your Nb-Ag complex structures before feature prediction.

`cd run_NbX`
`sh RosettaDock.sh`

# Limitations of NbX

1) **"Garbage in, garbage out". NbX is not a docking but a re-ranking method, which completely depends on the quality of the input Nb-Ag complex structures to suggest native-like solutions.**

    **Take action:** 
    
    * Use a docking algorithm that is well-tested on predicting native-like Nb-Ab complex structures, no matter how the docking method ranks them. 

    * We used Nb-Ag complex structures from ClusPro -> RosettaDock full-atom refinement (or equivalent/ better docking methods) to benchmark NbX.

2) **NbX was largely unable to model a single classification threshold that can generally applied to all tested Nb-Ag complexes to distinguish non-native-like (0) or native-like (1) Nb-Ag complex structures.**


    **Take action:** 
    
    * Descendingly sort the `mean_predicted_CAPRI_binary_proba` in `NbX_prediction.csv` (i.e. the mean native-like probablilty of the 5-fold validated NbX models). This is our NbX re-rank for you.

    * Do consider top ranks as more probable native-like Nb-Ag complex structures compared to the lower ranks.

    * Avoid using a single classification threshold of the absolute value of the probablity as non-native-like (0) or native-like (1) among unrelated Nb-Ag pairs.

    * Avoid comparing the absolute values of the probablity among unrelated Nb-Ag pairs.

    * 202206 Update: Use the following distributions of "mean_predicted_CAPRI_binary_proba" among 1) crystal 2) native-like and 3) non-native like Nb-Ag complex      structures to guide your selection of native-like Nb-Ag complex structures.
