import re
from datetime import datetime, UTC

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

from app.schemas import (
    BindingAnalysis,
    DescriptorSnapshot,
    DrugLikenessCheck,
    ExperimentalAdjustment,
    FactorScore,
    FeatureGroup,
    ModelEvidence,
    ModelFeatureContribution,
    ModelNeighbor,
    Mrgprx1Comparison,
    PredictRequest,
    PredictResponse,
    PredictionProbabilities,
    SimilarLigand,
    WeightedFactor,
)
from app.services.compound_service import search_compound
from app.services import ml_classifier
from app.services.mrgprx1_reference import compare_with_mrgprx1
from app.services.mrgprx2_reference import all_reference_ligands, structural_reference_ligands
from app.services.training_store import reference_ligands_from_labels

KNOWN_LIGANDS = [
    {
        "name": "Compound 48/80",
        "label": "agonist",
        "smiles": "CNCCC1=CC(=C(C=C1)OC)CC2=CC(=CC(=C2OC)CC3=C(C=CC(=C3)CCNC)OC)CCNC",
        "rationale": "Classic MRGPRX2 agonist polyamine oligomer (PubChem CID 2855) with multiple basic secondary amines and methoxy-aromatic rings characteristic of mast-cell degranulators.",
        "source": "curated",
    },
    {
        "name": "Substance P fragment",
        "label": "agonist",
        "smiles": "CC(=O)NC(CCCCN)C(=O)NC(CC1=CC=CC=C1)C(=O)N",
        "rationale": "Peptidic agonist fragment with cationic amine and multiple polar contacts consistent with activation.",
        "source": "curated",
    },
]


def _reference_ligands() -> list[dict]:
    """Curated references (seed list + the broader literature-derived MRGPRX2
    database) plus any user-submitted structure/label pairs. Merging all three
    means every new label a user submits immediately becomes part of the
    structural similarity search used below, alongside the growing curated set."""
    return KNOWN_LIGANDS + structural_reference_ligands() + reference_ligands_from_labels()

# Below ~45% Tanimoto similarity to the closest curated reference, the model
# has no real structural anchor for an MRGPRX2 agonist/antagonist call - any
# verdict at that point would be a guess dressed up as a confident result.
SIMILARITY_CUTOFF_PCT = 45.0

BASIC_NITROGEN = Chem.MolFromSmarts("[NX3;!$(N[C,S]=[O,S])]")
POSITIVELY_CHARGED = Chem.MolFromSmarts("[N+]")
AMIDE = Chem.MolFromSmarts("C(=O)N")
PHENOL = Chem.MolFromSmarts("c[OX2H]")


def _clamp(value: float, low: float = 0.0, high: float = 100.0) -> float:
    return max(low, min(high, value))


def _rdkit_mol(smiles: str) -> Chem.Mol | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return molecule

    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    if molecule is None:
        return None

    Chem.SanitizeMol(molecule, catchErrors=True)
    return molecule


def _descriptor_snapshot(molecule: Chem.Mol) -> DescriptorSnapshot:
    return DescriptorSnapshot(
        molecular_weight=float(Descriptors.MolWt(molecule)),
        logp=float(Crippen.MolLogP(molecule)),
        tpsa=float(rdMolDescriptors.CalcTPSA(molecule)),
        h_bond_donors=int(Lipinski.NumHDonors(molecule)),
        h_bond_acceptors=int(Lipinski.NumHAcceptors(molecule)),
        rotatable_bonds=int(Lipinski.NumRotatableBonds(molecule)),
        ring_count=int(rdMolDescriptors.CalcNumRings(molecule)),
    )


def _structure_score(descriptors: DescriptorSnapshot) -> float:
    score = 55.0

    if 120 <= descriptors.molecular_weight <= 650:
        score += 10
    else:
        score -= min(15, abs(descriptors.molecular_weight - 385) / 18)

    if -1 <= descriptors.logp <= 4.5:
        score += 10
    else:
        score -= min(15, abs(descriptors.logp - 1.8) * 4)

    if 20 <= descriptors.tpsa <= 180:
        score += 8
    else:
        score -= min(12, abs(descriptors.tpsa - 90) / 10)

    if descriptors.rotatable_bonds <= 12:
        score += 6
    if 1 <= descriptors.ring_count <= 4:
        score += 6

    return round(_clamp(score), 1)


def _pharmacophore_groups(molecule: Chem.Mol, descriptors: DescriptorSnapshot) -> list[FeatureGroup]:
    has_basic_n = molecule.HasSubstructMatch(BASIC_NITROGEN)
    has_positive_charge = molecule.HasSubstructMatch(POSITIVELY_CHARGED) or has_basic_n
    aromatic_ring = rdMolDescriptors.CalcNumAromaticRings(molecule) > 0
    amphipathic = aromatic_ring and (has_basic_n or descriptors.tpsa >= 40)
    hydrophobic_core = descriptors.logp >= 0 or descriptors.ring_count > 0

    groups: list[FeatureGroup] = []
    critical = [item for item, present in [
        ("Basic Nitrogen", has_basic_n),
        ("Cationic Charge", has_positive_charge),
    ] if present]
    important = [item for item, present in [
        ("Aromatic Ring", aromatic_ring),
        ("Amphipathic Nature", amphipathic),
    ] if present]
    supportive = [item for item, present in [
        ("Hydrophobic Core", hydrophobic_core),
    ] if present]

    if critical:
        groups.append(FeatureGroup(title="Critical Features", items=critical))
    if important:
        groups.append(FeatureGroup(title="Important Features", items=important))
    if supportive:
        groups.append(FeatureGroup(title="Supportive", items=supportive))
    return groups


def _pharmacophore_score(molecule: Chem.Mol, descriptors: DescriptorSnapshot) -> float:
    score = 20.0
    if molecule.HasSubstructMatch(BASIC_NITROGEN):
        score += 30
    if molecule.HasSubstructMatch(POSITIVELY_CHARGED) or molecule.HasSubstructMatch(BASIC_NITROGEN):
        score += 20
    if rdMolDescriptors.CalcNumAromaticRings(molecule) > 0:
        score += 15
    if rdMolDescriptors.CalcNumAromaticRings(molecule) > 0 and descriptors.tpsa >= 40:
        score += 10
    if descriptors.logp >= 0 or descriptors.ring_count > 0:
        score += 10
    if molecule.HasSubstructMatch(PHENOL):
        score += 5
    return round(_clamp(score), 1)


# Mirrors the residue tiers shown in the "Binding Site Residues" picker
# (Asp184/Glu164 anchor the ligand's cationic group via salt bridges - the
# best-documented MRGPRX2 contact; Tyr279/Phe170 form the aromatic subpocket;
# His259 sits in the more peripheral allosteric-adjacent region).
_CRITICAL_RESIDUES = {"asp184": "Asp184", "glu164": "Glu164"}
_IMPORTANT_RESIDUES = {"tyr279": "Tyr279", "phe170": "Phe170"}
_MODERATE_RESIDUES = {"his259": "His259"}
_KNOWN_RESIDUES = {**_CRITICAL_RESIDUES, **_IMPORTANT_RESIDUES, **_MODERATE_RESIDUES}


def _binding_site_adjustment(payload: PredictRequest) -> tuple[str, float] | None:
    """Bounded, proportional read of how strongly the user-asserted binding
    site evidence (selected pocket residues + binding mode) points to the
    compound actually occupying the MRGPRX2 pocket - not which way.

    Like a docking pose or an expression swing, "this ligand contacts
    Asp184/Glu164" is consistent with either an agonist or an antagonist
    locking into the same anchor point, so it's folded into the same
    engagement-strength calibration as the other experimental evidence (see
    `_apply_experimental_evidence`) rather than scored as its own free-floating
    factor disconnected from the model's actual agonist/antagonist call.
    """
    selected = set(payload.selected_residues)
    matched = selected & _KNOWN_RESIDUES.keys()

    magnitude = (
        len(selected & _CRITICAL_RESIDUES.keys()) * 2.0
        + len(selected & _IMPORTANT_RESIDUES.keys()) * 1.2
        + len(selected & _MODERATE_RESIDUES.keys()) * 0.8
    )
    if payload.binding_mode == "Orthosteric":
        magnitude += 1.5
    elif payload.binding_mode == "Allosteric":
        magnitude += 0.8

    magnitude = round(min(6.0, magnitude), 1)
    if magnitude <= 0:
        return None

    detail_parts = []
    if matched:
        detail_parts.append("/".join(_KNOWN_RESIDUES[residue] for residue in sorted(matched)))
    if payload.binding_mode and payload.binding_mode != "Unknown":
        detail_parts.append(payload.binding_mode.lower())
    detail = f" ({', '.join(detail_parts)})" if detail_parts else ""

    return f"Asserted MRGPRX2 pocket contacts{detail}", magnitude


def _composite_score(factors: list[FactorScore]) -> float:
    total_weight = sum(factor.weight for factor in factors)
    weighted = sum(factor.score * factor.weight for factor in factors)
    return round(weighted / total_weight, 1) if total_weight else 0.0


def _weighted_factor_breakdown(factors: list[FactorScore]) -> list[WeightedFactor]:
    """Each factor's share of the composite score, as a percentage - this is
    what lets the UI render the Composite Score as a weighted stacked bar."""
    contributions = [factor.score * factor.weight for factor in factors]
    total = sum(contributions)
    return [
        WeightedFactor(
            label=factor.label,
            weight=factor.weight,
            score=factor.score,
            contribution_pct=round((contribution / total) * 100, 1) if total else 0.0,
        )
        for factor, contribution in zip(factors, contributions)
    ]


def _structure_probabilities(model_evidence: dict) -> PredictionProbabilities:
    """The trained classifier's three-way class probabilities (the
    structure-grounded "prior"), as a chartable distribution.

    The classifier is trained on all three classes directly (see
    `ml_classifier.LABEL_TO_CLASS`: nonbinder / agonist / antagonist), so
    "nonbinder-like" is a real learned prediction here - not a fallback for
    when the model can't decide between the other two.
    """
    return PredictionProbabilities(
        nonbinder=round(model_evidence["nonbinder_probability"] / 100.0, 4),
        agonist=round(model_evidence["agonist_probability"] / 100.0, 4),
        antagonist=round(model_evidence["antagonist_probability"] / 100.0, 4),
    )


def _rank_prediction(probabilities: PredictionProbabilities) -> tuple[str, float]:
    """Turns a three-way probability distribution into the top label and a
    confidence score (the top class's probability, as a 0-100 percentage)."""
    ranked = sorted(
        [
            ("nonbinder-like", probabilities.nonbinder),
            ("agonist-like", probabilities.agonist),
            ("antagonist-like", probabilities.antagonist),
        ],
        key=lambda item: item[1],
        reverse=True,
    )
    label, top_probability = ranked[0]
    return label, round(top_probability * 100, 1)


# Caps how far experimental evidence can move the structure-grounded
# probabilities, in percentage points. The structure-based ML call stays the
# dominant signal - experimental data calibrates it within a bounded range
# rather than overriding it, so a single plausible-looking docking number
# can't make a clearly-nonbinder structure read as a strong binder.
MAX_EXPERIMENTAL_ADJUSTMENT_PCT = 20.0


def _assay_context_suffix(
    method: str | None, cell_line: str | None, concentration: float | None, time_hours: float | None
) -> str:
    """Renders the assay conditions (method, cell line, concentration,
    duration) a cell-based expression readout was produced under, as a
    trailing description suffix - e.g. "(4.2x, qPCR, HMC-1, 50µM, 24h)".

    These fields don't carry their own directional agonist/antagonist signal
    (a cell line name or a timepoint isn't "evidence of X"), so rather than
    inventing a numeric weight for e.g. "HMC-1 vs LAD2" - which would just be
    a fabricated relationship dressed up as domain knowledge - they're surfaced
    as visible, traceable provenance attached to the measurement they describe.
    Concentration is the one condition with a well-established, non-fabricated
    pharmacological consequence, so it additionally feeds into
    `_specificity_scale` below.
    """
    parts = [part for part in (method, cell_line) if part]
    if concentration is not None:
        parts.append(f"{concentration:.0f}µM")
    if time_hours is not None:
        parts.append(f"{time_hours:.0f}h")
    return f", {', '.join(parts)}" if parts else ""


def _specificity_scale(concentration: float | None) -> tuple[float, str | None]:
    """Discounts expression-based evidence recorded at high test
    concentrations, and explains why (or returns no discount/note).

    This is a standard, widely-taught GPCR pharmacology caveat rather than a
    guess: MRGPRX2 ligand studies typically run in the ~0.1-100µM range, and
    responses recorded well above that increasingly reflect non-specific,
    promiscuous membrane effects rather than MRGPRX2-specific engagement - so
    a large fold-change measured at e.g. 400µM is weaker evidence of genuine
    receptor engagement than the same fold-change at 10µM. The discount ramps
    linearly from full strength at <=100µM to half-strength at >=300µM.
    """
    if concentration is None or concentration <= 100.0:
        return 1.0, None
    scale = round(_clamp(1.0 - (concentration - 100.0) / 400.0, 0.5, 1.0), 2)
    note = (
        f"Tested at a high concentration ({concentration:.0f}µM) - results in "
        "this range can reflect non-specific membrane effects rather than "
        f"MRGPRX2-specific signalling, so this expression evidence is discounted "
        f"to {scale * 100:.0f}% strength."
    )
    return scale, note


def _experimental_adjustment_components(
    payload: PredictRequest,
) -> tuple[list[tuple[str, float]], str | None]:
    """Bounded, additive nudges (description, +percentage points) derived
    from experimental evidence (AutoDock docking score, MRGPRX2 mRNA/protein
    expression fold-change, with assay conditions attached as context), plus
    an optional caveat note about how those conditions affected the read.

    These describe how *strongly the compound engages MRGPRX2 at all* - they
    don't carry directional agonist-vs-antagonist information on their own
    (a strong docking pose, a large expression swing, or a confirmed pocket
    contact is consistent with either an activator or a blocker), so they're
    applied as reinforcement of whichever non-nonbinder class the
    structure-based prior already favours (see `_apply_experimental_evidence`).
    """
    components: list[tuple[str, float]] = []
    specificity_note: str | None = None

    binding_site_component = _binding_site_adjustment(payload)
    if binding_site_component is not None:
        components.append(binding_site_component)

    data = payload.experimental_data
    if data is None:
        return components, specificity_note

    if data.docking_score is not None:
        magnitude = _docking_adjustment_pct(data.docking_score)
        if magnitude > 0:
            components.append((f"AutoDock affinity ({data.docking_score:.1f} kcal/mol)", magnitude))

    # Concentration, cell line, and duration describe the cell-based assay
    # that produced the expression readouts below - they don't apply to the
    # docking score above, which is a purely computational result.
    scale, note = _specificity_scale(data.concentration)

    if data.mrna_fold_change is not None:
        magnitude = round(_expression_adjustment_pct(data.mrna_fold_change) * scale, 1)
        if magnitude > 0:
            context = _assay_context_suffix(data.mrna_method, data.cell_line, data.concentration, data.time_hours)
            components.append((f"MRGPRX2 mRNA fold-change ({data.mrna_fold_change:.1f}x{context})", magnitude))
            specificity_note = note

    if data.protein_fold_change is not None:
        magnitude = round(_expression_adjustment_pct(data.protein_fold_change) * scale, 1)
        if magnitude > 0:
            context = _assay_context_suffix(data.protein_method, data.cell_line, data.concentration, data.time_hours)
            components.append((f"MRGPRX2 protein fold-change ({data.protein_fold_change:.1f}x{context})", magnitude))
            specificity_note = note

    return components, specificity_note


def _docking_adjustment_pct(docking_score: float) -> float:
    """Continuous, proportional read of AutoDock affinity, in place of a
    tiered scale - a tiered cutoff (e.g. "-7.0 = Moderate") creates an
    artificial cliff where -6.99 and -7.01 land in different buckets despite
    being practically identical. Instead this ramps linearly: scores weaker
    than -5.0 kcal/mol carry no meaningful MRGPRX2-engagement signal (0pp),
    the adjustment grows proportionally from there, and caps at +8.0pp
    around -10.0 kcal/mol and beyond (very strong affinity)."""
    if docking_score >= -5.0:
        return 0.0
    return round(_clamp((-docking_score - 5.0) * 1.6, 0.0, 8.0), 1)


def _expression_adjustment_pct(fold_change: float) -> float:
    """Continuous, proportional read of mRNA/protein expression fold-change:
    no change (1.0x) carries no signal (0pp), the adjustment ramps linearly
    from there, and caps at +5.0pp around 4.3x and beyond - the same shape
    as the docking ramp, just scaled to a smaller maximum since expression
    fold-change is one step further removed from direct receptor engagement."""
    if fold_change <= 1.0:
        return 0.0
    return round(_clamp((fold_change - 1.0) * 1.5, 0.0, 5.0), 1)


def _apply_experimental_evidence(
    probabilities: PredictionProbabilities, payload: PredictRequest
) -> tuple[PredictionProbabilities, dict | None]:
    """Calibrates the structure-based probabilities with experimental
    evidence, within a bounded range, and returns a transparent breakdown of
    what moved and by how much (or `None` when there's nothing to apply).

    Docking/expression data says "this compound engages MRGPRX2" - it doesn't
    say *how* (an antagonist binds just as tightly as an agonist), so the
    adjustment reinforces whichever of agonist/antagonist the structure-based
    model already leans toward, borrowing probability mass from "nonbinder".
    The agonist-vs-antagonist call itself stays purely structure-grounded.
    """
    components, specificity_note = _experimental_adjustment_components(payload)
    if not components:
        return probabilities, None

    target_label = "antagonist" if probabilities.antagonist >= probabilities.agonist else "agonist"
    raw = {"agonist": probabilities.agonist, "antagonist": probabilities.antagonist, "nonbinder": probabilities.nonbinder}

    requested_delta_pct = min(MAX_EXPERIMENTAL_ADJUSTMENT_PCT, sum(delta for _, delta in components))
    requested_delta = requested_delta_pct / 100.0
    # Can't take more probability mass from "nonbinder" than it actually has.
    applied_delta = min(requested_delta, raw["nonbinder"])

    raw[target_label] += applied_delta
    raw["nonbinder"] -= applied_delta

    adjusted = PredictionProbabilities(
        agonist=round(raw["agonist"], 4),
        antagonist=round(raw["antagonist"], 4),
        nonbinder=round(raw["nonbinder"], 4),
    )

    breakdown = {
        "target_label": target_label,
        "structure_probability": round(getattr(probabilities, target_label) * 100, 1),
        "adjusted_probability": round(getattr(adjusted, target_label) * 100, 1),
        "applied_delta_pct": round(applied_delta * 100, 1),
        "components": [{"description": description, "delta_pct": delta} for description, delta in components],
        "specificity_note": specificity_note,
    }
    return adjusted, breakdown


def _best_reference_match(molecule: Chem.Mol) -> tuple[dict | None, float]:
    query_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(molecule, radius=2, nBits=2048)
    best = None
    best_score = -1.0

    for ligand in _reference_ligands():
        ligand_mol = _rdkit_mol(ligand["smiles"])
        if ligand_mol is None:
            continue
        ligand_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(ligand_mol, radius=2, nBits=2048)
        score = DataStructs.TanimotoSimilarity(query_fp, ligand_fp)
        if score > best_score:
            best = ligand
            best_score = score

    return best, max(best_score, 0.0)


def _most_similar_ligand(smiles: str) -> SimilarLigand:
    molecule = _rdkit_mol(smiles)
    if molecule is None:
        return SimilarLigand(
            name="Unknown",
            label="reference",
            similarity=0.0,
            rationale="No valid molecular graph was available for similarity scoring.",
            source="curated",
        )

    best, best_score = _best_reference_match(molecule)
    assert best is not None

    # Below this floor the "best" match is just the least-dissimilar entry in
    # a small reference set, not a meaningful structural relative - naming it
    # anyway (e.g. "14.5% similar to Levofloxacin") reads as a confident,
    # data-backed call when it isn't one.
    if best_score * 100 < SIMILARITY_CUTOFF_PCT:
        return SimilarLigand(
            name="No close structural match",
            label="reference",
            similarity=round(best_score * 100, 1),
            rationale=(
                "The nearest curated reference shares only "
                f"{best_score * 100:.1f}% structural similarity ({best['name']}) - "
                "too low to draw a meaningful comparison. The prediction below "
                "relies on the descriptor- and pharmacophore-based factors instead."
            ),
            source="curated",
        )

    return SimilarLigand(
        name=best["name"],
        label=best["label"],
        similarity=round(best_score * 100, 1),
        rationale=best["rationale"],
        source=best.get("source", "curated"),
    )


def _drug_likeness(descriptors: DescriptorSnapshot) -> list[DrugLikenessCheck]:
    return [
        DrugLikenessCheck(label="MW ≤ 500", passed=descriptors.molecular_weight <= 500),
        DrugLikenessCheck(label="LogP ≤ 5", passed=descriptors.logp <= 5),
        DrugLikenessCheck(label="HBD ≤ 5", passed=descriptors.h_bond_donors <= 5),
        DrugLikenessCheck(label="HBA ≤ 10", passed=descriptors.h_bond_acceptors <= 10),
    ]


def _normalize_reference_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _name_matched_reference(name: str) -> dict | None:
    """Exact, punctuation/case-insensitive name match against curated
    references that have no resolvable small-molecule structure.

    Some clinically/literature-documented MRGPRX2 ligands (e.g. EVO756) don't
    have a clean public SMILES to featurize. Forcing a guessed/garbled
    structure into the trained classifier would teach it a fabricated
    structure-activity relationship, so for these a name match against the
    curated reference list is the grounded answer instead - it's a documented
    literature/clinical call, not a structure-based guess.
    """
    target = _normalize_reference_name(name)
    if not target:
        return None
    for ligand in all_reference_ligands():
        if ligand.get("smiles"):
            continue
        if _normalize_reference_name(ligand.get("name", "")) == target:
            return ligand
    return None


def _confirmed_reference_response(payload: PredictRequest, reference: dict, smiles: str) -> PredictResponse:
    """Definitive, name-matched answer for a curated reference ligand whose
    structure can't be resolved into SMILES (see `_name_matched_reference`).

    This bypasses the structure-based ML classifier entirely - the call comes
    straight from the curated literature/clinical reference entry, which is
    the grounded source of truth when no real structure is available to learn
    from.
    """
    label = reference["label"]
    descriptors = DescriptorSnapshot(
        molecular_weight=0.0,
        logp=0.0,
        tpsa=0.0,
        h_bond_donors=0,
        h_bond_acceptors=0,
        rotatable_bonds=0,
        ring_count=0,
    )
    mrgprx1_comparison = Mrgprx1Comparison(**compare_with_mrgprx1(smiles))
    probability_map = {"agonist": 0.0, "antagonist": 0.0, "nonbinder": 0.0}
    probability_map[label] = 1.0
    likelihood_map = {"agonist": 0.0, "antagonist": 0.0}
    if label in likelihood_map:
        likelihood_map[label] = 100.0
    similar_ligand = SimilarLigand(
        name=reference["name"],
        label=label,
        similarity=100.0,
        rationale=reference.get("rationale", ""),
        source="curated",
    )
    docking_score = payload.experimental_data.docking_score if payload.experimental_data else None

    return PredictResponse(
        compound_name=reference["name"],
        receptor=payload.receptor,
        smiles=smiles,
        prediction=f"{label}-like",
        probabilities=PredictionProbabilities(**probability_map),
        descriptors=descriptors,
        docking_score=docking_score,
        interpretation=(
            f"\"{reference['name']}\" is a documented MRGPRX2 {label} from the "
            "curated literature/clinical reference list, identified here by "
            "name rather than by structure: no public small-molecule SMILES "
            "could be resolved for it, so running it through the trained "
            "classifier would mean teaching it a fabricated structure-activity "
            f"relationship from a guessed structure. {reference.get('rationale', '')} "
            "This is a confirmed reference call grounded in curated "
            "documentation, not a structure-based probability estimate."
        ),
        model_version="curated-name-match-v1",
        confidence=100.0,
        factor_analysis=[],
        weighted_factors=[],
        composite_score=0.0,
        binding_analysis=BindingAnalysis(
            agonist_likelihood=likelihood_map["agonist"],
            antagonist_likelihood=likelihood_map["antagonist"],
        ),
        feature_groups=[],
        similar_ligand=similar_ligand,
        mrgprx1_comparison=mrgprx1_comparison,
        drug_likeness=[],
        analyzed_at=datetime.now(UTC).isoformat(),
    )


def _unknown_structure_response(payload: PredictRequest, smiles: str) -> PredictResponse:
    """Graceful result for input that can't be resolved into a real structure.

    Every structure-input path (typed name, pasted SMILES, file/photo upload)
    should still "work" rather than surfacing a raw parse/lookup error - it
    just reports the compound as unknown so the rest of the UI renders
    normally and the user can try a different input.
    """
    descriptors = DescriptorSnapshot(
        molecular_weight=0.0,
        logp=0.0,
        tpsa=0.0,
        h_bond_donors=0,
        h_bond_acceptors=0,
        rotatable_bonds=0,
        ring_count=0,
    )
    mrgprx1_comparison = Mrgprx1Comparison(**compare_with_mrgprx1(smiles))
    similar_ligand = SimilarLigand(
        name="No close structural match",
        label="reference",
        similarity=0.0,
        rationale=(
            "The input couldn't be resolved into a valid molecular structure, "
            "so no structural comparison could be made."
        ),
        source="curated",
    )
    docking_score = payload.experimental_data.docking_score if payload.experimental_data else None

    return PredictResponse(
        compound_name="Unknown compound",
        receptor=payload.receptor,
        smiles=smiles,
        prediction="nonbinder-like",
        probabilities=PredictionProbabilities(agonist=0.0, antagonist=0.0, nonbinder=100.0),
        descriptors=descriptors,
        docking_score=docking_score,
        interpretation=(
            f"\"{payload.compound_name}\" couldn't be recognized as a known "
            "compound or parsed as a valid chemical structure, so no "
            "structure-based analysis could be run. Try a different name, a "
            "valid SMILES string, or a clearer structure file/photo."
        ),
        model_version="heuristic-multifactor-v3",
        confidence=0.0,
        factor_analysis=[],
        weighted_factors=[],
        composite_score=0.0,
        binding_analysis=BindingAnalysis(agonist_likelihood=0.0, antagonist_likelihood=0.0),
        feature_groups=[],
        similar_ligand=similar_ligand,
        mrgprx1_comparison=mrgprx1_comparison,
        drug_likeness=[],
        analyzed_at=datetime.now(UTC).isoformat(),
    )


def predict_activity(payload: PredictRequest) -> PredictResponse:
    try:
        compound = search_compound(payload.compound_name)
    except LookupError:
        reference = _name_matched_reference(payload.compound_name)
        if reference is not None:
            return _confirmed_reference_response(payload, reference, payload.compound_name.strip())
        return _unknown_structure_response(payload, payload.compound_name.strip())

    molecule = _rdkit_mol(compound.smiles)
    if molecule is None:
        reference = _name_matched_reference(payload.compound_name) or _name_matched_reference(compound.name)
        if reference is not None:
            return _confirmed_reference_response(payload, reference, compound.smiles)
        return _unknown_structure_response(payload, compound.smiles)

    descriptors = _descriptor_snapshot(molecule)
    mrgprx1_comparison = Mrgprx1Comparison(**compare_with_mrgprx1(compound.smiles))

    factor_analysis = [
        FactorScore(label="Structure Analysis", weight=1.0, score=_structure_score(descriptors)),
        FactorScore(label="Pharmacophore Match", weight=1.5, score=_pharmacophore_score(molecule, descriptors)),
    ]
    # Experimental evidence (docking score, mRNA/protein expression) and
    # user-asserted binding-site evidence (selected pocket residues, binding
    # mode) no longer get their own free-floating "Experimental Support" /
    # "Binding Alignment" factors here - that produced disconnected stories
    # (e.g. a high heuristic score for a middling docking number, or a
    # "Binding Alignment: 89" that never moved the actual agonist/antagonist
    # call) alongside the structure-grounded ML verdict. They now flow
    # through a single path: `_apply_experimental_evidence` calibrates the
    # trained classifier's probabilities directly (see below and the
    # "Calibrated by your experimental data" panel), so there's one
    # consistent account of what the evidence does to the call.
    factor_analysis.append(
        FactorScore(
            label="MRGPRX2 vs MRGPRX1 Selectivity",
            weight=0.7,
            score=mrgprx1_comparison.mrgprx2_selectivity_score,
        )
    )
    composite_score = _composite_score(factor_analysis)
    weighted_factors = _weighted_factor_breakdown(factor_analysis)
    similar_ligand = _most_similar_ligand(compound.smiles)
    feature_groups = _pharmacophore_groups(molecule, descriptors)
    docking_score = payload.experimental_data.docking_score if payload.experimental_data else None

    # The agonist/antagonist call itself comes from a Random Forest trained on
    # every curated + user-labeled MRGPRX2 ligand (see ml_classifier.py) -
    # this is what makes the verdict a learned, evidence-backed prediction
    # rather than a hand-tuned if/else heuristic. `evidence` is None only when
    # there isn't yet enough labeled data of both classes to fit a model, or
    # the structure can't be featurized - in which case no grounded call can
    # be made and we say so plainly instead of guessing.
    evidence = ml_classifier.classify(compound.smiles, _reference_ligands())

    if evidence is None:
        prediction = "indeterminate"
        confidence = 0.0
        probabilities = PredictionProbabilities(agonist=0.0, antagonist=0.0, nonbinder=0.0)
        interpretation = (
            f"\"{compound.name}\" couldn't be given a grounded Agonist/Antagonist "
            "call: there isn't yet enough labeled MRGPRX2 data covering both "
            "classes for the trained classifier to learn a reliable boundary "
            "for this structure, or the structure couldn't be featurized. "
            "Label this compound (and others like it) under \"Teach the "
            "Classifier\" to help the model learn this region of chemical space."
        )
        return PredictResponse(
            compound_name=compound.name,
            receptor=payload.receptor,
            smiles=compound.smiles,
            prediction=prediction,
            probabilities=probabilities,
            descriptors=descriptors,
            docking_score=docking_score,
            interpretation=interpretation,
            model_version="random-forest-v1",
            confidence=confidence,
            factor_analysis=factor_analysis,
            weighted_factors=weighted_factors,
            composite_score=composite_score,
            binding_analysis=BindingAnalysis(agonist_likelihood=0.0, antagonist_likelihood=0.0),
            feature_groups=feature_groups,
            similar_ligand=similar_ligand,
            mrgprx1_comparison=mrgprx1_comparison,
            drug_likeness=_drug_likeness(descriptors),
            analyzed_at=datetime.now(UTC).isoformat(),
        )

    structure_probabilities = _structure_probabilities(evidence)
    probabilities, experimental_adjustment = _apply_experimental_evidence(structure_probabilities, payload)
    prediction, confidence = _rank_prediction(probabilities)

    model_evidence = ModelEvidence(
        model_type=evidence["model_type"],
        trained_on=evidence["trained_on"],
        nonbinder_probability=evidence["nonbinder_probability"],
        agonist_probability=evidence["agonist_probability"],
        antagonist_probability=evidence["antagonist_probability"],
        nearest_neighbors=[ModelNeighbor(**neighbor) for neighbor in evidence["nearest_neighbors"]],
        top_features=[ModelFeatureContribution(**feature) for feature in evidence["top_features"]],
        summary=evidence["summary"],
        experimental_adjustment=(
            ExperimentalAdjustment(**experimental_adjustment) if experimental_adjustment else None
        ),
    )

    nearest = model_evidence.nearest_neighbors[0] if model_evidence.nearest_neighbors else None
    interpretation = (
        f"The trained classifier puts \"{compound.name}\" at "
        f"{evidence['nonbinder_probability']:.1f}% nonbinder / "
        f"{evidence['agonist_probability']:.1f}% agonist / "
        f"{evidence['antagonist_probability']:.1f}% antagonist likelihood, "
        f"learned from {evidence['trained_on']} labeled MRGPRX2 compounds "
        "(see Model Evidence for the nearest labeled neighbors and the "
        "features driving this call)."
    )
    if nearest is not None:
        interpretation += (
            f" Its closest labeled neighbor is {nearest.name} "
            f"({nearest.similarity:.1f}% similarity, labeled {nearest.label})."
        )
    interpretation += (
        f" Closest curated reference ligand for structural context: "
        f"{similar_ligand.name} ({similar_ligand.similarity:.1f}% similarity, "
        f"{similar_ligand.source} reference)."
    )
    if experimental_adjustment is not None:
        component_summaries = ", ".join(
            f"{component['description']} (+{component['delta_pct']:.1f}pp)"
            for component in experimental_adjustment["components"]
        )
        interpretation += (
            f" Experimental evidence you supplied ({component_summaries}) nudged "
            f"the {experimental_adjustment['target_label']} call from "
            f"{experimental_adjustment['structure_probability']:.1f}% (structure-only) "
            f"to {experimental_adjustment['adjusted_probability']:.1f}% - this calibrates "
            "the structure-grounded probability within a bounded range rather than "
            "overriding it; the agonist-vs-antagonist call itself remains "
            "structure-based."
        )
        if experimental_adjustment["specificity_note"]:
            interpretation += f" {experimental_adjustment['specificity_note']}"

    return PredictResponse(
        compound_name=compound.name,
        receptor=payload.receptor,
        smiles=compound.smiles,
        prediction=prediction,
        probabilities=probabilities,
        descriptors=descriptors,
        docking_score=docking_score,
        interpretation=interpretation,
        model_version="random-forest-v1",
        confidence=confidence,
        factor_analysis=factor_analysis,
        weighted_factors=weighted_factors,
        composite_score=composite_score,
        binding_analysis=BindingAnalysis(
            agonist_likelihood=round(probabilities.agonist * 100, 1),
            antagonist_likelihood=round(probabilities.antagonist * 100, 1),
        ),
        feature_groups=feature_groups,
        similar_ligand=similar_ligand,
        mrgprx1_comparison=mrgprx1_comparison,
        model_evidence=model_evidence,
        drug_likeness=_drug_likeness(descriptors),
        analyzed_at=datetime.now(UTC).isoformat(),
    )
