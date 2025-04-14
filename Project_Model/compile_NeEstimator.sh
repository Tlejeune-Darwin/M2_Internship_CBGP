#!/bin/bash

# Cloner le dépôt NeEstimator
git clone https://github.com/molecularfisherieslab/NeEstimator2.X.git
cd NeEstimator2.X

# Compiler NeEstimator
make

# Copier l'exécutable dans le dossier bin
cp Ne2L ../bin/NeEstimator_linux
