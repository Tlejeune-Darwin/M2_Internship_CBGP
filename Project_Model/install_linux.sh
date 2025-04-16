#!/bin/bash

# Mettre à jour les paquets
sudo apt update && sudo apt upgrade -y

# Installer les dépendances nécessaires
sudo apt install -y build-essential cmake qtbase5-dev zlib1g-dev git python3 python3-pip

# Installer les packages Python requis
pip3 install -r requirements.txt
