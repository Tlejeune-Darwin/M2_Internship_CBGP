#!/bin/bash

echo "🔎 Vérification des exécutables dans ./bin/"

# Vérifie SLiM
if [ ! -f bin/slim_linux ]; then
    echo "❌ SLiM (slim_linux) manquant dans bin/"
    exit 1
else
    chmod +x bin/slim_linux
    echo "✅ SLiM prêt"
fi

# Vérifie NeEstimator
if [ ! -f bin/NeEstimator_linux ]; then
    echo "❌ NeEstimator (NeEstimator_linux) manquant dans bin/"
    exit 1
else
    chmod +x bin/NeEstimator_linux
    echo "✅ NeEstimator prêt"
fi

echo "📁 Les exécutables sont prêts à être utilisés."
