#!/bin/bash

# Vérifier que les exécutables existent
if [ ! -f Bin/slim ]; then
    echo "SLiM n'est pas compilé. Exécutez ./compile_slim.sh"
    exit 1
fi

if [ ! -f Bin/Ne2x ]; then
    echo "NeEstimator n'est pas compilé. Exécutez ./compile_neestimator.sh"
    exit 1
fi

# Lancer la simulation SLiM
./bin/slim_linux slim_scripts/Model_microsat.slim

# Lancer NeEstimator
./bin/NeEstimator_linux config/option
