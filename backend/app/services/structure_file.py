from rdkit import Chem

from app.schemas import StructureFileParseResponse


def _mol_to_smiles(molecule: Chem.Mol) -> str | None:
    try:
        Chem.SanitizeMol(molecule, catchErrors=True)
        stripped = Chem.RemoveHs(molecule)
        return Chem.MolToSmiles(stripped)
    except Exception:
        return None


def _parse_sdf_or_mol(text: str) -> tuple[Chem.Mol | None, str | None]:
    """Returns (molecule, name) for the first record in an SDF/MOL block.

    SDF "title" lines (the first line of a molfile) are commonly the
    compound name/ID - we surface that as a hint, but the structure itself
    comes from the parsed atom/bond block, not the title text.
    """
    first_block = text.split("$$$$")[0]
    molecule = Chem.MolFromMolBlock(first_block, sanitize=False)
    if molecule is None:
        return None, None

    title = first_block.splitlines()[0].strip() if first_block.splitlines() else ""
    return molecule, (title or None)


def _parse_mol2(text: str) -> tuple[Chem.Mol | None, str | None]:
    molecule = Chem.MolFromMol2Block(text, sanitize=False)
    if molecule is None:
        return None, None

    name = None
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.strip().upper() == "@<TRIPOS>MOLECULE" and index + 1 < len(lines):
            name = lines[index + 1].strip() or None
            break
    return molecule, name


def parse_structure_file(filename: str, content: bytes) -> StructureFileParseResponse:
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        text = content.decode("latin-1", errors="ignore")

    lowered = filename.lower()
    if lowered.endswith(".mol2"):
        molecule, name = _parse_mol2(text)
    else:
        molecule, name = _parse_sdf_or_mol(text)

    if molecule is None:
        return StructureFileParseResponse(
            name=None,
            smiles=None,
            message=(
                f"Couldn't read a molecular structure from \"{filename}\". "
                "Please upload a valid SDF (.sdf) or MOL2 (.mol2) file, or "
                "paste the SMILES string directly."
            ),
        )

    smiles = _mol_to_smiles(molecule)
    if not smiles:
        return StructureFileParseResponse(
            name=name,
            smiles=None,
            message=(
                f"Read \"{filename}\" but couldn't convert its structure into "
                "a valid molecule. Please paste the SMILES string directly."
            ),
        )

    return StructureFileParseResponse(
        name=name,
        smiles=smiles,
        message=f"Loaded the structure from \"{filename}\".",
    )
