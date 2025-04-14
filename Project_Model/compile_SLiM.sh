#!/bin/bash

# Cloner le dépôt SLiM
git clone https://github.com/MesserLab/SLiM.git
cd SLiM

# Créer un dossier de build
mkdir build && cd build

# Générer les fichiers de build avec CMake
cmake ..

# Compiler SLiM
make -j$(nproc)

# Copier l'exécutable dans le dossier bin
cp slim ../../bin/slim
