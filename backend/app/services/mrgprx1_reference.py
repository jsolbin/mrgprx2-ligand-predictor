"""Curated MRGPRX1 reference ligands and MRGPRX1-vs-MRGPRX2 binding-site
differences, used to sharpen MRGPRX2 agonist/antagonist calls by checking
whether a query structure looks like it would prefer the paralogous MRGPRX1
pocket instead.

Residue numbers follow the cryo-EM MRGPRX1 structures (e.g. PDB 8JGG /
Yang et al., Nat Commun 2023): the orthosteric acidic pair Asp177(5.38) /
Glu157(4.60) is conserved with MRGPRX2's Asp184(5.38) / Glu164(4.60), but
MRGPRX1 additionally carries a TM3-TM4 side pocket (Tyr99/Pro100/Trp158)
that recognises the "RF/RY" dipeptide motif found at the C-terminus of
peptidergic ligands such as BAM8-22 - a feature MRGPRX2's broader, more
promiscuous pocket lacks. Small cationic amphiphiles (e.g. Compound 48/80)
activate MRGPRX2 but are poor MRGPRX1 ligands, while rigid bicyclic
amines such as chloroquine are selective MRGPRX1 agonists.
"""

from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors

MRGPRX1_LIGANDS = [
    {
        "name": "Chloroquine",
        "label": "agonist",
        "smiles": "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12",
        "rationale": "Selective MRGPRX1 small-molecule agonist - rigid bicyclic aminoquinoline with a single basic side-chain amine; lacks the dipeptide RF/RY motif MRGPRX1's TM3-TM4 pocket also tolerates, and is a poor MRGPRX2 ligand.",
    },
    {
        "name": "Beta-Alanine",
        "label": "agonist",
        "smiles": "NCCC(=O)O",
        "rationale": "Minimal endogenous MRGPRX1 agonist (zwitterionic amino acid); far below the size/aromaticity profile MRGPRX2 ligands typically need.",
    },
    {
        "name": "BAM8-22 (C-terminal fragment)",
        "label": "agonist",
        "smiles": "NC(Cc1ccccc1)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(CCCNC(N)=N)C(=O)NC(Cc1ccccc1)C(N)=O",
        "rationale": "Peptidergic MRGPRX1 agonist; its C-terminal Phe-Arg/Tyr-Phe 'RF/RY'-like dipeptide motif inserts into the MRGPRX1-specific TM3-TM4 side pocket (Tyr99/Pro100/Trp158) that MRGPRX2 does not possess.",
    },
]

# Conserved orthosteric acidic pair (numbering per Yang et al. 2023 / PDB 8JGG).
SHARED_ORTHOSTERIC_RESIDUES = {
    "mrgprx1": ["Asp177 (5.38)", "Glu157 (4.60)"],
    "mrgprx2": ["Asp184 (5.38)", "Glu164 (4.60)"],
    "note": "Both receptors anchor the ligand's basic/cationic group with this acidic pair, which is why cationic amines drive activity at either paralog.",
}

# Pocket features that diverge between the two receptors and can be read
# off a 2D structure as coarse selectivity signals.
MRGPRX1_SELECTIVE_POCKET = {
    "residues": ["Tyr99 (3.29)", "Pro100 (3.30)", "Trp158 (4.61)"],
    "feature": "TM3-TM4 'RF/RY' side pocket",
    "structural_signature": "C-terminal aromatic-basic dipeptide motif (e.g. Phe-Arg / Tyr-Phe-NH2) typical of peptidergic MRGPRX1 agonists such as BAM8-22",
}

MRGPRX2_SELECTIVE_POCKET = {
    "residues": ["Tyr279 (6.51)", "Phe170", "His259"],
    "feature": "Wide, lipophilic orthosteric vestibule",
    "structural_signature": "Small cationic amphiphile with a flexible basic amine plus an aromatic/phenolic ring (e.g. Compound 48/80) - tolerated by MRGPRX2's more promiscuous, shallower pocket but poorly accommodated by MRGPRX1's narrower RF/RY-shaped pocket",
}

RF_RY_MOTIF = Chem.MolFromSmarts("[CX3](=O)[NX3][CX4][CX3](=O)[NX3]")
GUANIDINE = Chem.MolFromSmarts("NC(N)=N")
BASIC_AMINE = Chem.MolFromSmarts("[NX3;!$(N[C,S]=[O,S]);!$(N=*)]")
AROMATIC = Chem.MolFromSmarts("c1ccccc1")
METHOXY_AROMATIC = Chem.MolFromSmarts("c[OX2][CH3]")

# Compound 48/80's actual structure (polyamine oligomer with methoxy-aromatic
# rings) - the "small cationic amphiphile" signal below should only fire when
# a query genuinely resembles *this*, not just any molecule that happens to
# carry one basic amine and one benzene ring (a near-universal drug-like trait
# that was previously triggering the signal on almost every input).
_COMPOUND_48_80_SMILES = (
    "CNCCC1=CC(=C(C=C1)OC)CC2=CC(=CC(=C2OC)CC3=C(C=CC(=C3)CCNC)OC)CCNC"
)


def _rdkit_mol(smiles: str) -> Chem.Mol | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return molecule
    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    if molecule is None:
        return None
    Chem.SanitizeMol(molecule, catchErrors=True)
    return molecule


def _fingerprint(molecule: Chem.Mol):
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(molecule, radius=2, nBits=2048)


def compare_with_mrgprx1(smiles: str) -> dict:
    """Return a structural read-out of how the query compound's features map
    onto the MRGPRX1 vs MRGPRX2 orthosteric pockets, plus a 0-100 score that
    rewards MRGPRX2-leaning (and penalises MRGPRX1-leaning) structures so it
    can feed into the agonist/antagonist composite.
    """
    molecule = _rdkit_mol(smiles)
    if molecule is None:
        return {
            "mrgprx2_selectivity_score": 50.0,
            "closest_mrgprx1_ligand": None,
            "mrgprx1_similarity": 0.0,
            "shared_residues": SHARED_ORTHOSTERIC_RESIDUES,
            "mrgprx1_pocket": MRGPRX1_SELECTIVE_POCKET,
            "mrgprx2_pocket": MRGPRX2_SELECTIVE_POCKET,
            "signals": [],
            "summary": "No valid molecular graph available for cross-receptor comparison.",
        }

    query_fp = _fingerprint(molecule)
    best = None
    best_score = -1.0
    for ligand in MRGPRX1_LIGANDS:
        ligand_mol = _rdkit_mol(ligand["smiles"])
        if ligand_mol is None:
            continue
        score = DataStructs.TanimotoSimilarity(query_fp, _fingerprint(ligand_mol))
        if score > best_score:
            best, best_score = ligand, score

    has_rf_ry_motif = molecule.HasSubstructMatch(RF_RY_MOTIF) and molecule.HasSubstructMatch(GUANIDINE)
    has_basic_amine = molecule.HasSubstructMatch(BASIC_AMINE)
    has_aromatic = molecule.HasSubstructMatch(AROMATIC)
    amide_bonds = rdMolDescriptors.CalcNumAmideBonds(molecule)
    is_peptidic = amide_bonds >= 3

    # A loose "has *an* amine and *a* benzene ring" test matches almost any
    # drug-like molecule, which was making the Compound-48/80-like signal show
    # up on unrelated structures. Require it to actually look like Compound
    # 48/80: multiple basic amines, multiple methoxy-aromatic rings, and
    # measurable fingerprint similarity to the real reference structure.
    basic_amine_count = len(molecule.GetSubstructMatches(BASIC_AMINE))
    methoxy_aromatic_count = len(molecule.GetSubstructMatches(METHOXY_AROMATIC))
    compound_48_80_mol = _rdkit_mol(_COMPOUND_48_80_SMILES)
    compound_48_80_similarity = (
        DataStructs.TanimotoSimilarity(query_fp, _fingerprint(compound_48_80_mol))
        if compound_48_80_mol is not None
        else 0.0
    )
    is_small_amphiphile = (
        not is_peptidic
        and has_basic_amine
        and has_aromatic
        and basic_amine_count >= 2
        and methoxy_aromatic_count >= 1
        and compound_48_80_similarity >= 0.2
    )

    signals: list[str] = []
    score = 50.0

    if has_rf_ry_motif or (is_peptidic and amide_bonds >= 4):
        signals.append(
            "Peptidic backbone with an aromatic-basic C-terminal motif resembles MRGPRX1's RF/RY-recognition pocket (Tyr99/Pro100/Trp158) more than MRGPRX2's shallower vestibule."
        )
        score -= 16
    if is_small_amphiphile:
        signals.append(
            f"Polyamine, methoxy-aromatic scaffold shares real structural similarity ({compound_48_80_similarity * 100:.1f}%) with Compound 48/80, the prototypical small-molecule MRGPRX2 agonist that MRGPRX2's wider Tyr279/Phe170/His259 vestibule favours over MRGPRX1's narrower pocket."
        )
        score += 14
    if best is not None and best["label"] == "agonist" and best_score >= 0.35:
        signals.append(
            f"Notable similarity ({best_score * 100:.1f}%) to the MRGPRX1-selective agonist {best['name']} - weigh the MRGPRX2 call down slightly unless other evidence (docking, expression) is strong."
        )
        score -= min(12, best_score * 20)
    if not signals:
        signals.append(
            "No strong MRGPRX1-selectivity signature detected; structure is consistent with either paralog at the conserved Asp/Glu orthosteric anchor."
        )

    score = max(0.0, min(100.0, score))

    return {
        "mrgprx2_selectivity_score": round(score, 1),
        "closest_mrgprx1_ligand": best["name"] if best else None,
        "mrgprx1_similarity": round(best_score * 100, 1) if best else 0.0,
        "shared_residues": SHARED_ORTHOSTERIC_RESIDUES,
        "mrgprx1_pocket": MRGPRX1_SELECTIVE_POCKET,
        "mrgprx2_pocket": MRGPRX2_SELECTIVE_POCKET,
        "signals": signals,
        "summary": (
            "Compares the query against curated MRGPRX1 agonists (chloroquine, beta-alanine, BAM8-22) "
            "and the divergent TM3-TM4 'RF/RY' side pocket that distinguishes MRGPRX1 from MRGPRX2, "
            "to flag structures that may be cross-reactive or MRGPRX1-selective rather than MRGPRX2-selective."
        ),
    }
