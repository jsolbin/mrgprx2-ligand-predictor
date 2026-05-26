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

## Repository Structure

```text
frontend/  Next.js UI for compound search and prediction results
backend/   FastAPI API for compound lookup and prediction logic
```

## Current Status

This repository is now initialized with a working full-stack skeleton:

- `frontend/`: Next.js App Router project with a compound search page, prediction panel, and probability chart
- `backend/`: FastAPI app with `GET /health`, `GET /compound/search`, and `POST /predict`
- `backend/data/mrgprx2_ligands.sample.csv`: starter dataset template
- Mock compound and prediction data so the UI can be wired before PubChem, RDKit, and model training are added

The backend does not yet call PubChem or RDKit. It returns deterministic placeholder data for MVP integration.

## Local Setup

Recommended runtime:

- Node.js `22.x` LTS
- Python `3.11` or `3.12`

### 1. Frontend

```bash
cd frontend
cp .env.example .env.local
npm install
npm run dev
```

Frontend runs at `http://localhost:3000`.

If you are using Node `25.x`, the experimental Web Storage API can break Next.js dev mode with errors like `localStorage.getItem is not a function`. The frontend scripts already disable that feature, but Node `22.x` LTS is still the recommended version.

Useful frontend commands:

```bash
npm run lint
npm run lint:fix
npm run format
npm run format:check
```

### 2. Backend

```bash
cd backend
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env
uvicorn app.main:app --reload
```

Backend runs at `http://localhost:8000`.

## Implemented API

### `GET /health`

Returns API health status.

### `GET /compound/search?name=...`

Returns:

- compound name
- PubChem CID placeholder
- SMILES
- InChIKey
- molecular weight
- LogP
- TPSA
- SVG structure placeholder

### `POST /predict`

Request:

```json
{
  "compound_name": "Ciprofloxacin",
  "receptor": "MRGPRX2"
}
```

Returns:

- final prediction label
- agonist / antagonist / nonbinder probabilities
- descriptor snapshot
- docking score placeholder
- interpretation

## Recommended Next Steps

1. Replace mock compound lookup with PubChem PUG REST integration.
2. Add RDKit descriptor calculation and real 2D structure rendering.
3. Add a curated ligand dataset in `backend/data/`.
4. Implement a rule-based baseline using fingerprint similarity.
5. Train and load a first `RandomForest` classifier.
6. Add docking jobs as a separate worker pipeline instead of synchronous API calls.
