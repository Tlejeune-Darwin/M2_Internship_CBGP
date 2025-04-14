#!/bin/bash

echo "🔎 Vérification des exécutables dans ./bin/"

# Vérifie SLiM
if [ ! -f bin/slim ]; then
    echo "❌ SLiM (slim) manquant dans bin/"
    exit 1
else
    chmod +x bin/slim_linux
    echo "✅ SLiM prêt"
fi

# Vérifie NeEstimator
if [ ! -f bin/Ne2x ]; then
    echo "❌ NeEstimator (Ne2x) manquant dans bin/"
    exit 1
else
    chmod +x bin/NeEstimator_linux
    echo "✅ NeEstimator prêt"
fi

echo "📁 Les exécutables sont prêts à être utilisés."
