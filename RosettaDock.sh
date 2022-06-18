#!/bin/bash

ROSETTA_PATH="/data/cltam/script/rosetta_src_2021.16.61629_bundle/main"
ROSETTA_BIN=$ROSETTA_PATH"/source/bin"
ROSETTA_DB=$ROSETTA_PATH"/database"
ROSETTA_NSTRUCT=100

# RosettaDock
cd input
ls -1 *.pdb | sort -V > pdb_list
cp ../RosettaDock_input/* ./
$ROSETTA_BIN/rosetta_scripts.linuxgccrelease --database $ROSETTA_DB -nstruct $ROSETTA_NSTRUCT @flags