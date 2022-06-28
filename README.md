# NbX

This repository contains the NbX models for the re-ranking of nanobody–antigen binding poses.

Tam, C.; Kumar, A.; Zhang, K.Y.J. NbX: Machine Learning-Guided Re-Ranking of Nanobody–Antigen Binding Poses. Pharmaceuticals 2021, 14, 968. https://doi.org/10.3390/ph14100968

![postersample](https://user-images.githubusercontent.com/51283097/174424802-a14cc780-1a59-4527-bc64-1321f63bb98c.png)

---

# Model

The five models inside the `model` folder are the  XGBoost models shown on Figure 1 of our NbX paper.

(i.e. the 5-fold validated models trained and tested on Nb-Ag complexes where pairwise Ag structural alignment quality score < 0.9 to minimize train-test information leakage).

| Model | PR-AUC |
| ------------- | ------------- |
| model_0001 | 0.229 |
| model_0002 | 0.169 |
| model_0003 | 0.349 (best) |
| model_0004 | 0.205 |
| model_0005 | 0.276 |

---

# How to run NbX

To run an example: do [step1](https://github.com/johnnytam100/NbX#step-1--install-pre-requisites), [step2](https://github.com/johnnytam100/NbX#step-2--change-paths) then go directly to [step5](https://github.com/johnnytam100/NbX#step-5--change-directory).

To run with your Nb-Ag complex structures: start from [step1](https://github.com/johnnytam100/NbX#step-1--install-pre-requisites).

## Step 1 : clone and create environment

`git clone https://github.com/johnnytam100/NbX.git`

`conda env create -f NbX/environment.yml`

`conda activate nbx`

which should have the following libraries installed

`pip install biopandas==0.4.1`

`pip install xgboost==0.90`

`pip install scikit-learn==0.22.2.post1`

`pip install joblib==1.1.0`

`pip install dill==0.3.5.1`



Then, manually install the following

PyRosetta (https://www.pyrosetta.org/downloads,

NbX was tested with installtion with `pyrosetta-2022.23+release.f1e0f6d7bf7-cp38-cp38-linux_x86_64.whl`)

FoldX (https://foldxsuite.crg.eu/)

DockQ (optional, https://github.com/bjornwallner/DockQ)

Rosetta (optional, https://new.rosettacommons.org/demos/latest/tutorials/install_build/install_build)


## Step 2 : change paths

Change the following paths inside `NbX_feature_prep.py`

`path_to_python = "/home/cltam/anaconda3/envs/nbx/bin/python"`

`path_to_foldx = "/data/cltam/script/FoldX/foldx_20221231"`

`path_to_dockq = "/data/cltam/script/DockQ/"` (optional, be careful not to omit the last `/` in this path)


## Step 3 : renumber nanobody

Before any docking, please renumber your nanobody with PyIgClassify (http://dunbrack2.fccc.edu/pyigclassify/)

If your Nb (or Nb-Ag complex) structure is confidential and you don't want to submit to a webserver: 

modify the CDRs start and end residue numbers (search `"CDR1_start_residue"`) inside `NbX_feature_prep.py`.

## Step 4 : copy your Nb-Ag complex structures to `run_NbX` folder

`cp (path to your Nb-Ag complex structures .pdb) ./run_NbX`

## Step 5 : change directory

`cd run_NbX`

`cp -r ../model ../NbX_feature_prep.py ../aaDescriptors.csv ../NbX_predict.py ./`

## Step 6 : feature preparation

### (option 1,  without DockQ)

`python NbX_feature_prep.py --antigen_chain A --antibody_chain H`

### (option 2,  with DockQ)

`python NbX_feature_prep.py --antigen_chain A --antibody_chain H --native 6oq8_complex.pdb`


## Step 7 : predict

`python NbX_predict.py`


## Step 8 : analyze results

**Important concept in NbX re-ranking:**

Descendingly sort the `mean_predicted_CAPRI_binary_proba` in `NbX_prediction.csv`, we get the following results

![投影片1](https://user-images.githubusercontent.com/51283097/174423865-865a8b73-d382-4080-b080-8fa49e5b2a44.PNG)


# RosettaDock (optional)

To mimic the NbX benchmark setting, you can perform RosettaDock refinement of your Nb-Ag complex structures before feature preparation.

step 1 : change `ROSETTA_PATH` in `RosettaDock.sh` to your Rosetta path.

step 2 : `sh RosettaDock.sh`

---

# Limitations of NbX

1) **"Garbage in, garbage out". NbX is not a docking but a re-ranking method, which completely depends on the quality of the input Nb-Ag complex structures to suggest native-like solutions.**

    **Take action:** 
    
    * Use a docking algorithm that is well-tested on predicting native-like Nb-Ab complex structures, no matter how the docking method ranks them. 

    * We used Nb-Ag complex structures from ClusPro -> RosettaDock full-atom refinement to benchmark NbX. Please use equivalent or better docking methods.

2) **NbX is largely unable to model a single classification threshold that generally applies to all tested Nb-Ag complexes to distinguish non-native-like (0) or native-like (1) Nb-Ag complex structures.**


    **Take action:** 
    
    * Descendingly sort the `mean_predicted_CAPRI_binary_proba` in `NbX_prediction.csv` (i.e. the mean native-like probablilty of the 5-fold validated NbX models). This is our NbX re-rank for you.

    * Do consider top ranks as more probable native-like Nb-Ag complex structures compared to the lower ranks.

    * Avoid applying a single classification threshold on the absolute value of the probability to distinguish non-native-like (0) or native-like (1) among unrelated Nb-Ag pairs.

    * Avoid comparing the absolute value of the probability predicted among unrelated Nb-Ag pairs.

    * 202206 Update: Use the following distributions of `mean_predicted_CAPRI_binary_proba` among 1) crystal 2) native-like and 3) non-native-like Nb-Ag complex      structures to guide your selection of native-like Nb-Ag complex structures.

   ![image](https://user-images.githubusercontent.com/51283097/174424674-a5d30058-64aa-460d-b2a5-e6346901170d.png)

