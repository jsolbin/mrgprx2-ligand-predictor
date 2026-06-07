"""Trained ML classifier for the MRGPRX2 nonbinder/agonist/antagonist call.

The rest of the prediction pipeline (factor analysis, MRGPRX1 selectivity,
drug-likeness, ...) is transparent rule-based scoring, which is fine as
*supporting* evidence - but the final verdict was previously a hand-tuned
if/else heuristic dressed up as a "model". This module instead fits a real
scikit-learn classifier on every labeled structure available (curated
reference ligands + every compound a user has labeled via "Teach the
Classifier"), across all three classes (nonbinder / agonist / antagonist -
see `LABEL_TO_CLASS`), and reports its prediction together with the
evidence behind it: the nearest labeled neighbours and which molecular
features the model leaned on.

Because a user can add labels at any time, the classifier is refit
on-demand whenever the labeled set changes (cheap: a few dozen molecules,
sub-second to fit) rather than trained once and frozen.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from sklearn.ensemble import RandomForestClassifier

FINGERPRINT_RADIUS = 2
FINGERPRINT_BITS = 512

DESCRIPTOR_NAMES = [
    "Molecular weight",
    "LogP (lipophilicity)",
    "Topological polar surface area",
    "Hydrogen-bond donor count",
    "Hydrogen-bond acceptor count",
    "Rotatable bond count",
    "Ring count",
    "Basic nitrogen count",
    "Aromatic ring count",
    "Positively-charged nitrogen count",
]

_BASIC_NITROGEN = Chem.MolFromSmarts("[NX3;!$(N[C,S]=[O,S])]")
_POSITIVE_NITROGEN = Chem.MolFromSmarts("[N+]")

# Three-way split requested for the MRGPRX2 call: a structure is either a
# Nonbinder (doesn't engage the receptor at all), an Agonist, or an
# Antagonist - a binary agonist-vs-antagonist model can't express "neither",
# so it has to be a real class the model is trained to recognise.
LABEL_TO_CLASS = {"nonbinder": 0, "agonist": 1, "antagonist": 2}
CLASS_TO_LABEL = {value: key for key, value in LABEL_TO_CLASS.items()}

_cache: dict | None = None


def _rdkit_mol(smiles: str) -> Chem.Mol | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return molecule
    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    if molecule is None:
        return None
    Chem.SanitizeMol(molecule, catchErrors=True)
    return molecule


def _morgan_fp(molecule: Chem.Mol):
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(
        molecule, FINGERPRINT_RADIUS, nBits=FINGERPRINT_BITS
    )


def _fingerprint_array(fingerprint) -> np.ndarray:
    array = np.zeros((FINGERPRINT_BITS,), dtype=np.float64)
    DataStructs.ConvertToNumpyArray(fingerprint, array)
    return array


def _descriptor_array(molecule: Chem.Mol) -> np.ndarray:
    return np.array(
        [
            Descriptors.MolWt(molecule),
            Crippen.MolLogP(molecule),
            rdMolDescriptors.CalcTPSA(molecule),
            Lipinski.NumHDonors(molecule),
            Lipinski.NumHAcceptors(molecule),
            Lipinski.NumRotatableBonds(molecule),
            rdMolDescriptors.CalcNumRings(molecule),
            len(molecule.GetSubstructMatches(_BASIC_NITROGEN)),
            rdMolDescriptors.CalcNumAromaticRings(molecule),
            len(molecule.GetSubstructMatches(_POSITIVE_NITROGEN)),
        ],
        dtype=np.float64,
    )


def _feature_vector(molecule: Chem.Mol, fingerprint) -> np.ndarray:
    return np.concatenate([_fingerprint_array(fingerprint), _descriptor_array(molecule)])


def _training_set_key(reference_compounds: list[dict]) -> tuple:
    return tuple(sorted((entry["smiles"], entry["label"]) for entry in reference_compounds))


def _fit(reference_compounds: list[dict]) -> dict | None:
    rows: list[np.ndarray] = []
    labels: list[int] = []
    indexed_compounds: list[dict] = []

    for entry in reference_compounds:
        class_index = LABEL_TO_CLASS.get(entry.get("label"))
        if class_index is None:
            continue
        molecule = _rdkit_mol(entry["smiles"])
        if molecule is None:
            continue
        fingerprint = _morgan_fp(molecule)
        rows.append(_feature_vector(molecule, fingerprint))
        labels.append(class_index)
        indexed_compounds.append({**entry, "_fingerprint": fingerprint})

    # A classifier needs every class represented to learn decision
    # boundaries between all three - with too little data it would just
    # memorise noise rather than learn a generalisable rule.
    if len(set(labels)) < 3 or len(labels) < 9:
        return None

    features = np.vstack(rows)
    targets = np.array(labels)

    model = RandomForestClassifier(
        n_estimators=300,
        max_depth=5,
        min_samples_leaf=2,
        random_state=7,
    )
    model.fit(features, targets)

    return {"model": model, "compounds": indexed_compounds}


def _get_fitted(reference_compounds: list[dict]) -> dict | None:
    global _cache
    key = _training_set_key(reference_compounds)
    if _cache is None or _cache["key"] != key:
        fitted = _fit(reference_compounds)
        _cache = {"key": key, "fitted": fitted}
    return _cache["fitted"]


def classify(smiles: str, reference_compounds: list[dict]) -> dict | None:
    """Runs the trained classifier on `smiles` and returns its verdict plus
    the evidence behind it (nearest labeled neighbours, most influential
    features), or None if there isn't enough labeled data to train on yet
    or the structure can't be parsed."""
    molecule = _rdkit_mol(smiles)
    if molecule is None:
        return None

    fitted = _get_fitted(reference_compounds)
    if fitted is None:
        return None

    model: RandomForestClassifier = fitted["model"]
    compounds = fitted["compounds"]

    query_fingerprint = _morgan_fp(molecule)
    features = _feature_vector(molecule, query_fingerprint).reshape(1, -1)

    probabilities = model.predict_proba(features)[0]
    class_position = {label: index for index, label in enumerate(model.classes_)}

    def _probability_for(label: str) -> float:
        position = class_position.get(LABEL_TO_CLASS[label])
        return float(probabilities[position]) * 100.0 if position is not None else 0.0

    nonbinder_probability = _probability_for("nonbinder")
    agonist_probability = _probability_for("agonist")
    antagonist_probability = _probability_for("antagonist")

    neighbors = sorted(
        (
            (DataStructs.TanimotoSimilarity(query_fingerprint, entry["_fingerprint"]), entry)
            for entry in compounds
        ),
        key=lambda item: item[0],
        reverse=True,
    )[:3]
    nearest_neighbors = [
        {
            "name": entry["name"],
            "label": entry["label"],
            "similarity": round(similarity * 100, 1),
            "source": entry.get("source", "curated"),
        }
        for similarity, entry in neighbors
    ]

    importances = model.feature_importances_
    descriptor_importances = importances[FINGERPRINT_BITS:]
    fingerprint_importance_total = float(importances[:FINGERPRINT_BITS].sum()) * 100.0
    ranked_descriptor_indices = np.argsort(descriptor_importances)[::-1]
    top_features = [
        {
            "description": DESCRIPTOR_NAMES[index],
            "importance": round(float(descriptor_importances[index]) * 100.0, 1),
        }
        for index in ranked_descriptor_indices[:3]
        if descriptor_importances[index] > 0
    ]
    if fingerprint_importance_total > 0:
        top_features.append(
            {
                "description": "Substructure fingerprint patterns (combined)",
                "importance": round(fingerprint_importance_total, 1),
            }
        )
    top_features.sort(key=lambda item: item["importance"], reverse=True)

    agonist_count = sum(1 for entry in compounds if entry["label"] == "agonist")
    antagonist_count = sum(1 for entry in compounds if entry["label"] == "antagonist")
    nonbinder_count = sum(1 for entry in compounds if entry["label"] == "nonbinder")
    user_count = sum(1 for entry in compounds if entry.get("source") == "user")

    return {
        "model_type": "Random Forest classifier (scikit-learn)",
        "trained_on": len(compounds),
        "nonbinder_probability": round(nonbinder_probability, 1),
        "agonist_probability": round(agonist_probability, 1),
        "antagonist_probability": round(antagonist_probability, 1),
        "nearest_neighbors": nearest_neighbors,
        "top_features": top_features,
        "summary": (
            f"Trained on {len(compounds)} labeled MRGPRX2 compounds "
            f"({agonist_count} agonist, {antagonist_count} antagonist, "
            f"{nonbinder_count} nonbinder"
            + (f", {user_count} from your own labels" if user_count else "")
            + ") using Morgan-fingerprint and physicochemical features across "
            "all three classes. This call is the model's learned probability "
            "for the structure - grounded in the nearest labeled neighbours "
            "and feature importances shown below, not a hand-written rule."
        ),
    }
