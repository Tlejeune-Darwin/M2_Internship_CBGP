#!/bin/bash

# Vérifier que les exécutables existent
if [ ! -f bin/slim_linux ]; then
    echo "SLiM n'est pas compilé. Exécutez ./compile_slim.sh"
    exit 1
fi

if [ ! -f bin/NeEstimator_linux ]; then
    echo "NeEstimator n'est pas compilé. Exécutez ./compile_neestimator.sh"
    exit 1
fi

# Lancer la simulation SLiM
./bin/slim_linux slim_scripts/Model_microsat.slim

# Lancer NeEstimator
./bin/NeEstimator_linux config/option
