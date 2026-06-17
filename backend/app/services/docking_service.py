"""AutoDock Vina docking service for MRGPRX2 ligand affinity estimation.

Workflow:
  SMILES → 3D conformer (RDKit ETKDGv3) → PDBQT (Meeko)
      → Vina subprocess (receptor: 7VDH chain R, orthosteric pocket)
      → best affinity (kcal/mol)

Grid box is centred on the 6IB ligand in 7VDH (compound 48/80 fragment),
covering the MRGPRX2 orthosteric binding pocket (25 × 25 × 25 Å).
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

RECEPTOR_PDBQT = _BASE / "data" / "receptor" / "mrgprx2_receptor.pdbqt"

# Grid centred on compound-48/80 fragment (6IB) in PDB 7VDH chain R
GRID_CENTER = (99.35, 65.30, 82.77)
GRID_SIZE = (25.0, 25.0, 25.0)

# Docking accuracy: 8 is publication-quality; 4 is fast/preview
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
    """Convert SMILES to a PDBQT string via RDKit + Meeko."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) == -1:
        # ETKDG failed – try random coords
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())
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


def _parse_vina_output(stdout: str) -> float | None:
    """Extract best affinity (mode 1) from Vina stdout table."""
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                return float(parts[1])
            except ValueError:
                continue
    return None


def run_docking(
    smiles: str,
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    num_modes: int = 9,
) -> dict:
    """Run AutoDock Vina and return the best docking affinity.

    Returns:
        {
          "affinity_kcal_mol": float,   # best pose energy (negative = better)
          "num_modes": int,             # poses found
          "warning": str | None,        # non-fatal issues
        }
    Raises:
        RuntimeError on Vina failure.
        ValueError on bad SMILES or missing receptor.
    """
    if not RECEPTOR_PDBQT.exists():
        raise FileNotFoundError(
            f"Receptor PDBQT not found: {RECEPTOR_PDBQT}. "
            "Run scripts/prep_receptor.py to generate it."
        )

    vina = _vina_binary()
    warning = None

    with tempfile.TemporaryDirectory() as tmp:
        ligand_path = _smiles_to_pdbqt(smiles, tmp)
        out_path = os.path.join(tmp, "out.pdbqt")

        cx, cy, cz = GRID_CENTER
        sx, sy, sz = GRID_SIZE

        cmd = [
            str(vina),
            "--receptor", str(RECEPTOR_PDBQT),
            "--ligand", ligand_path,
            "--center_x", str(cx),
            "--center_y", str(cy),
            "--center_z", str(cz),
            "--size_x", str(sx),
            "--size_y", str(sy),
            "--size_z", str(sz),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--out", out_path,
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            raise RuntimeError("Vina timed out after 120 s")

        if result.returncode != 0:
            raise RuntimeError(
                f"Vina exited with code {result.returncode}:\n{result.stderr}"
            )

        output = result.stdout + result.stderr
        affinity = _parse_vina_output(output)
        if affinity is None:
            raise RuntimeError(
                f"Could not parse Vina affinity from output:\n{output}"
            )

        # Count poses in output file
        poses = 0
        if os.path.exists(out_path):
            with open(out_path) as f:
                poses = f.read().count("MODEL")

        return {
            "affinity_kcal_mol": affinity,
            "num_modes": poses,
            "warning": warning,
        }
