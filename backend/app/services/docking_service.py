"""AutoDock Vina dual-state docking service for MRGPRX2.

Docks the ligand against TWO receptor conformations and returns ΔΔScore:

  Active state  : PDB 7VDH chain R (Gq-coupled MRGPRX2 + C48/80 fragment,
                  cryo-EM 2.9 Å). Grid centred on the 6IB ligand.
  Inactive state: AlphaFold2 monomer v6 (Q96LB1), which adopts an
                  inactive-like conformation with inward TM6 and closed
                  G-protein interface.

ΔΔScore = affinity_active − affinity_inactive  (both in kcal/mol)

Interpretation:
  ΔΔScore ≪ 0  → ligand prefers active state  → agonist-consistent
  ΔΔScore ≈ 0  → state-indifferent binding    → cannot distinguish
  ΔΔScore ≫ 0  → ligand prefers inactive state → antagonist-consistent

The absolute docking score (active) still correlates with binding affinity
(binder vs nonbinder), while ΔΔScore carries the directional signal.
"""

from __future__ import annotations

import os
import platform
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

_BASE = Path(__file__).resolve().parent.parent.parent

# Active state: PDB 7VDH chain R (Gq-coupled, agonist-bound)
ACTIVE_RECEPTOR_PDBQT = _BASE / "data" / "receptor" / "mrgprx2_receptor.pdbqt"
ACTIVE_GRID_CENTER = (99.35, 65.30, 82.77)   # centroid of 6IB ligand in 7VDH

# Inactive state: AlphaFold2 Q96LB1 (inactive-like, TM6 inward)
INACTIVE_RECEPTOR_PDBQT = _BASE / "data" / "receptor" / "mrgprx2_inactive.pdbqt"
INACTIVE_GRID_CENTER = (8.08, -3.06, -4.32)  # centroid of key pocket residues

GRID_SIZE = (25.0, 25.0, 25.0)
DEFAULT_EXHAUSTIVENESS = 8


def _vina_binary() -> Path:
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "darwin":
        name = "vina_mac_aarch64" if "arm" in machine else "vina_mac_x86_64"
    else:
        name = "vina_linux_x86_64"
    path = _BASE / "bin" / name
    if not path.exists():
        raise FileNotFoundError(f"Vina binary not found: {path}")
    path.chmod(0o755)
    return path


def _smiles_to_pdbqt(smiles: str, tmp_dir: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # Some chemically valid SMILES (e.g. chromone/coumarin/lactone derivatives
        # with aromatic ring notation like c(=O) in a pyranone) fail RDKit's
        # strict sanitization.  Retry with partial sanitization that skips the
        # final property/valence check — enough for most heterocyclic rings.
        mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol_raw is not None:
            try:
                Chem.SanitizeMol(
                    mol_raw,
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
                )
                mol = mol_raw
            except Exception:
                pass
    if mol is None:
        raise ValueError(
            "3D structure preparation failed — RDKit cannot parse this SMILES. "
            "The molecule likely uses non-standard aromatic lactone/chromone notation "
            "(e.g. c(=O) inside an aromatic ring). "
            "Please enter the docking score manually."
        )
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) == -1:
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv2()) == -1:
            raise ValueError(
                "3D conformation generation failed for this SMILES "
                "(conformer embedding returned no valid poses). "
                "Please enter the docking score manually."
            )
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

    prep = MoleculePreparation()
    mol_setups = prep.prepare(mol)
    if not mol_setups:
        raise ValueError("Meeko returned no molecule setup")
    pdbqt_str, _, _ = PDBQTWriterLegacy.write_string(mol_setups[0])

    ligand_path = os.path.join(tmp_dir, "ligand.pdbqt")
    with open(ligand_path, "w") as f:
        f.write(pdbqt_str)
    return ligand_path


def _parse_vina_best(stdout: str) -> float | None:
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                return float(parts[1])
            except ValueError:
                continue
    return None


def _run_vina(
    receptor: Path,
    ligand_path: str,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    out_path: str,
    exhaustiveness: int,
    num_modes: int,
    vina: Path,
) -> tuple[float, int]:
    cx, cy, cz = center
    sx, sy, sz = size
    cmd = [
        str(vina),
        "--receptor", str(receptor),
        "--ligand", ligand_path,
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--out", out_path,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        raise RuntimeError("Vina timed out after 120 s")

    if result.returncode != 0:
        raise RuntimeError(
            f"Vina exited with code {result.returncode}:\n{result.stderr}"
        )

    affinity = _parse_vina_best(result.stdout + result.stderr)
    if affinity is None:
        raise RuntimeError(f"Could not parse Vina affinity:\n{result.stdout}")

    poses = 0
    if os.path.exists(out_path):
        with open(out_path) as f:
            poses = f.read().count("MODEL")
    return affinity, poses


def run_docking(
    smiles: str,
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    num_modes: int = 9,
) -> dict:
    """Run dual-state AutoDock Vina and return affinity + ΔΔScore.

    Returns:
        {
          "affinity_kcal_mol": float,     # active-state best pose (kcal/mol)
          "inactive_affinity_kcal_mol": float,
          "delta_delta_score": float,      # active − inactive (kcal/mol)
          "num_modes": int,
          "active_receptor": str,
          "inactive_receptor": str,
          "warning": str | None,
        }

    ΔΔScore < 0 → active-state preference → agonist signal
    ΔΔScore > 0 → inactive-state preference → antagonist/nonbinder signal
    """
    for path, label in [
        (ACTIVE_RECEPTOR_PDBQT, "active"),
        (INACTIVE_RECEPTOR_PDBQT, "inactive"),
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Receptor PDBQT ({label}) not found: {path}"
            )

    vina = _vina_binary()

    with tempfile.TemporaryDirectory() as tmp:
        ligand_path = _smiles_to_pdbqt(smiles, tmp)

        active_affinity, active_poses = _run_vina(
            ACTIVE_RECEPTOR_PDBQT, ligand_path,
            ACTIVE_GRID_CENTER, GRID_SIZE,
            os.path.join(tmp, "out_active.pdbqt"),
            exhaustiveness, num_modes, vina,
        )
        inactive_affinity, _ = _run_vina(
            INACTIVE_RECEPTOR_PDBQT, ligand_path,
            INACTIVE_GRID_CENTER, GRID_SIZE,
            os.path.join(tmp, "out_inactive.pdbqt"),
            exhaustiveness, num_modes, vina,
        )

    delta_delta = round(active_affinity - inactive_affinity, 3)

    return {
        "affinity_kcal_mol": active_affinity,
        "inactive_affinity_kcal_mol": inactive_affinity,
        "delta_delta_score": delta_delta,
        "num_modes": active_poses,
        "active_receptor": "7VDH chain R (Gq-coupled active state)",
        "inactive_receptor": "AlphaFold2 Q96LB1 (inactive-like)",
        "warning": None,
    }
