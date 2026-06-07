import base64
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem.Draw import rdMolDraw2D


def _clean_svg(svg: str) -> str:
    cleaned = svg.replace("<?xml version='1.0' encoding='iso-8859-1'?>", "")
    cleaned = cleaned.replace("svg:", "")
    cleaned = cleaned.replace("fill:#FFFFFF;stroke:none", "fill:transparent;stroke:none")
    cleaned = cleaned.replace("width='520px'", "width='520'")
    cleaned = cleaned.replace("height='380px'", "height='380'")
    if "preserveAspectRatio=" not in cleaned:
        cleaned = cleaned.replace(
            "viewBox='0 0 520 380'>",
            "viewBox='0 0 520 380' preserveAspectRatio='xMidYMid meet'>",
        )
    return cleaned.strip()


def _build_molecule(smiles: str) -> Chem.Mol | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return molecule

    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    if molecule is None:
        return None

    Chem.SanitizeMol(molecule, catchErrors=True)
    return molecule


def smiles_to_svg(smiles: str) -> str | None:
    molecule = _build_molecule(smiles)
    if molecule is None:
        return None

    AllChem.Compute2DCoords(molecule)
    prepared = rdMolDraw2D.PrepareMolForDrawing(molecule, kekulize=False)
    drawer = rdMolDraw2D.MolDraw2DSVG(520, 380)
    draw_options = drawer.drawOptions()
    draw_options.addAtomIndices = False
    draw_options.bondLineWidth = 2.6
    draw_options.padding = 0.12
    draw_options.multipleBondOffset = 0.18
    draw_options.fixedBondLength = 42
    drawer.DrawMolecule(prepared)
    drawer.FinishDrawing()
    return _clean_svg(drawer.GetDrawingText())


def smiles_to_png_markup(smiles: str) -> str | None:
    png_bytes = smiles_to_png_bytes(smiles)
    if png_bytes is None:
        return None

    encoded = base64.b64encode(png_bytes).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def smiles_to_png_bytes(smiles: str) -> bytes | None:
    molecule = _build_molecule(smiles)
    if molecule is None:
        return None

    AllChem.Compute2DCoords(molecule)
    image = MolToImage(molecule, size=(520, 380), kekulize=False)
    buffer = BytesIO()
    image.save(buffer, format="PNG")
    return buffer.getvalue()


def smiles_to_descriptors(smiles: str) -> dict[str, float] | None:
    molecule = _build_molecule(smiles)
    if molecule is None:
        return None

    return {
        "molecular_weight": float(Descriptors.MolWt(molecule)),
        "logp": float(Crippen.MolLogP(molecule)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(molecule)),
    }
