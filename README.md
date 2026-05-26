# MRGPRX2 Ligand Predictor

A web-based research platform for predicting MRGPRX2 ligand activity using chemical structure information, PubChem data, RDKit molecular descriptors, machine learning models, and structure-based docking workflows.

## Project Overview

This project aims to build a structure-based prediction system for MRGPRX2 ligands.

The final goal is to allow users to enter a chemical name or SMILES string, retrieve its molecular structure from PubChem, calculate molecular descriptors using RDKit, and predict whether the compound is likely to act as an MRGPRX2 agonist, antagonist, inhibitor, or nonbinder.

The system is designed as a full-stack application with a Next.js frontend and a FastAPI backend.

## Main Features

- Search compounds by chemical name or SMILES
- Retrieve compound information from PubChem
- Display 2D molecular structures
- Calculate molecular descriptors using RDKit
- Predict MRGPRX2 ligand activity
- Compare compounds against known agonist and antagonist lists
- Support machine learning-based classification
- Integrate docking scores from AutoDock Vina
- Manage curated MRGPRX2 ligand datasets
