#!/bin/bash

# === Étape 1 : Vérifier Python ===
echo "🔍 Vérification de Python..."
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 n'est pas installé. Veuillez l'installer d'abord."
    exit 1
fi

# === Étape 2 : Installer les dépendances Python ===
echo "📦 Installation des dépendances Python..."
pip3 install -r requirements.txt

# === Étape 3 : Donner les droits d'exécution aux exécutables ===
echo "⚙️ Préparation des exécutables..."
chmod +x bin/slim
chmod +x bin/NeEstimator

# === Étape 4 : Exécuter le script principal ===
echo "🚀 Lancement de la simulation..."
python3 python_scripts/main_test.py

# === Étape 5 : Nettoyage (optionnel) ===
echo "🧹 Nettoyage des fichiers temporaires..."
# Uncomment if needed
# python3 python_scripts/clean_outputs.py

