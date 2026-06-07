"""Persistence for user-submitted agonist/antagonist labels.

When a user tells the app "this SMILES is an agonist/antagonist", we record
it here (structure + label). Future predictions load this store and treat
each entry as an additional reference ligand for similarity scoring, so the
classifier's structural pattern-matching improves as more compounds are
labeled - without needing a full retraining pipeline.

Backed by a real database (see `app.db`) rather than a flat JSON file, so
concurrent writes are transactionally safe and the data can live in a
hosted Postgres instance that survives redeploys in production.
"""

from __future__ import annotations

from datetime import datetime, UTC

from rdkit import Chem
from sqlalchemy import select
from sqlalchemy.exc import NoResultFound

from app.db import get_session
from app.models import LabeledCompound

VALID_LABELS = {"agonist", "antagonist", "nonbinder"}


def _rdkit_mol(smiles: str) -> Chem.Mol | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return molecule
    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    if molecule is None:
        return None
    Chem.SanitizeMol(molecule, catchErrors=True)
    return molecule


def add_label(smiles: str, label: str, name: str | None = None, note: str | None = None) -> dict:
    normalized_smiles = smiles.strip()
    normalized_label = label.strip().lower()

    if normalized_label not in VALID_LABELS:
        raise ValueError(f"Label must be one of {sorted(VALID_LABELS)}, got '{label}'.")

    molecule = _rdkit_mol(normalized_smiles)
    if molecule is None:
        raise ValueError(f"'{smiles}' is not a parseable SMILES string.")

    canonical_smiles = Chem.MolToSmiles(molecule)

    with get_session() as session:
        try:
            entry = session.execute(
                select(LabeledCompound).where(LabeledCompound.smiles == canonical_smiles)
            ).scalar_one()
        except NoResultFound:
            entry = LabeledCompound(smiles=canonical_smiles)
            session.add(entry)

        entry.name = (name or "User-labeled compound").strip() or "User-labeled compound"
        entry.label = normalized_label
        entry.note = (note or "").strip() or None
        entry.source = "user"
        entry.submitted_at = datetime.now(UTC)

        session.commit()
        session.refresh(entry)
        return entry.to_dict()


def list_labels() -> list[dict]:
    with get_session() as session:
        entries = session.execute(
            select(LabeledCompound).order_by(LabeledCompound.submitted_at.asc())
        ).scalars().all()
        return [entry.to_dict() for entry in entries]


def reference_ligands_from_labels() -> list[dict]:
    """Adapt stored user labels into the same shape as the curated
    KNOWN_LIGANDS reference list, so they can be merged for similarity
    matching and structure-based score nudging.
    """
    ligands = []
    for entry in list_labels():
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
