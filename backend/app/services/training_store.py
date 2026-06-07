"""Lightweight persistence for user-submitted agonist/antagonist labels.

When a user tells the app "this SMILES is an agonist/antagonist", we record
it here (structure + label). Future predictions load this store and treat
each entry as an additional reference ligand for similarity scoring, so the
classifier's structural pattern-matching improves as more compounds are
labeled - without needing a full retraining pipeline.
"""

from __future__ import annotations

import json
import threading
from datetime import datetime, UTC
from pathlib import Path

from rdkit import Chem

STORE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "user_labels.json"
VALID_LABELS = {"agonist", "antagonist", "nonbinder"}

_lock = threading.Lock()


def _rdkit_mol(smiles: str) -> Chem.Mol | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return molecule
    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    if molecule is None:
        return None
    Chem.SanitizeMol(molecule, catchErrors=True)
    return molecule


def _read_all() -> list[dict]:
    if not STORE_PATH.exists():
        return []
    try:
        with STORE_PATH.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return []
    return data if isinstance(data, list) else []


def _write_all(entries: list[dict]) -> None:
    STORE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with STORE_PATH.open("w", encoding="utf-8") as handle:
        json.dump(entries, handle, indent=2)


def add_label(smiles: str, label: str, name: str | None = None, note: str | None = None) -> dict:
    normalized_smiles = smiles.strip()
    normalized_label = label.strip().lower()

    if normalized_label not in VALID_LABELS:
        raise ValueError(f"Label must be one of {sorted(VALID_LABELS)}, got '{label}'.")

    molecule = _rdkit_mol(normalized_smiles)
    if molecule is None:
        raise ValueError(f"'{smiles}' is not a parseable SMILES string.")

    canonical_smiles = Chem.MolToSmiles(molecule)

    entry = {
        "name": (name or "User-labeled compound").strip() or "User-labeled compound",
        "smiles": canonical_smiles,
        "label": normalized_label,
        "note": (note or "").strip() or None,
        "source": "user",
        "submitted_at": datetime.now(UTC).isoformat(),
    }

    with _lock:
        entries = _read_all()
        entries = [item for item in entries if item.get("smiles") != canonical_smiles]
        entries.append(entry)
        _write_all(entries)

    return entry


def list_labels() -> list[dict]:
    with _lock:
        return list(_read_all())


def reference_ligands_from_labels() -> list[dict]:
    """Adapt stored user labels into the same shape as the curated
    KNOWN_LIGANDS reference list, so they can be merged for similarity
    matching and structure-based score nudging.
    """
    ligands = []
    for entry in _read_all():
        if entry.get("label") not in VALID_LABELS:
            continue
        ligands.append(
            {
                "name": entry["name"],
                "label": entry["label"],
                "smiles": entry["smiles"],
                "rationale": entry.get("note")
                or f"User-labeled {entry['label']} - structural pattern learned from prior submissions.",
                "source": "user",
            }
        )
    return ligands
