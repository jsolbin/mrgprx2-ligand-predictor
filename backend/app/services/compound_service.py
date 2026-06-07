import httpx
from rdkit import Chem

from app.schemas import CompoundSearchResponse
from app.services.curated_activity import display_name_for_cid, resolve_curated_cid
from app.services.mock_data import MOCK_COMPOUNDS, get_demo_svg
from app.services.pubchem_service import (
    fetch_pubchem_compound_by_cid,
    fetch_pubchem_png_bytes_by_name,
    fetch_pubchem_png_bytes_by_smiles,
    search_pubchem_by_name,
    search_pubchem_by_smiles,
)
from app.services.rdkit_renderer import (
    smiles_to_descriptors,
    smiles_to_png_markup,
    smiles_to_png_bytes,
    smiles_to_svg,
)


def _rdkit_mol(smiles: str) -> Chem.Mol | None:
    """Lenient SMILES parse matching predict_service/training_store/etc.

    Strict `Chem.MolFromSmiles` rejects some valid-looking structures over
    issues like kekulization (e.g. certain fused-ring lactone SMILES) that
    RDKit can still build and sanitize when given the chance to recover via
    `sanitize=False` + `catchErrors=True`. Using the same lenient parse here
    keeps "is this a valid structure?" consistent across search/predict/label
    - otherwise the same SMILES can be "valid" for labeling but "Unknown
    compound" for search/predict, which is exactly the kind of inconsistency
    that left a user unable to get any result for a structure they could label.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return molecule
    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    if molecule is None:
        return None
    Chem.SanitizeMol(molecule, catchErrors=True)
    return molecule


def _looks_like_smiles(query: str) -> bool:
    normalized = query.strip()
    if not normalized:
        return False

    allowed_characters = set(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789[]()=#@+-\\/."
    )
    has_smiles_tokens = any(character in normalized for character in "[]()=#@\\/0123456789")
    return has_smiles_tokens and all(character in allowed_characters for character in normalized)


def _build_manual_smiles_response(smiles: str) -> CompoundSearchResponse:
    molecule = _rdkit_mol(smiles.strip())
    normalized_smiles = Chem.MolToSmiles(molecule) if molecule is not None else smiles.strip()
    structure_svg = smiles_to_svg(normalized_smiles) or smiles_to_png_markup(
        normalized_smiles
    )
    descriptors = smiles_to_descriptors(normalized_smiles) or {
        "molecular_weight": 0.0,
        "logp": 0.0,
        "tpsa": 0.0,
    }
    return CompoundSearchResponse(
        name="Manual SMILES Input",
        cid=None,
        smiles=normalized_smiles,
        inchi_key=None,
        molecular_weight=descriptors["molecular_weight"],
        logp=descriptors["logp"],
        tpsa=descriptors["tpsa"],
        structure_svg=structure_svg,
        source="input",
    )


def unknown_compound_response(query: str) -> CompoundSearchResponse:
    """Graceful fallback for input that can't be resolved or parsed.

    Rather than surfacing a raw lookup/parse error to the user, every
    structure-input path (typed name, pasted SMILES, file/photo upload)
    should still "work" - it just reports the compound as unknown so the
    rest of the UI can render normally.
    """
    normalized = query.strip()
    is_smiles = _looks_like_smiles(normalized)
    molecule = _rdkit_mol(normalized) if is_smiles else None
    smiles = Chem.MolToSmiles(molecule) if molecule is not None else ""
    structure_svg = (
        smiles_to_svg(smiles) or smiles_to_png_markup(smiles) if smiles else None
    )
    descriptors = (smiles_to_descriptors(smiles) if smiles else None) or {
        "molecular_weight": 0.0,
        "logp": 0.0,
        "tpsa": 0.0,
    }
    return CompoundSearchResponse(
        name="Unknown compound",
        cid=None,
        smiles=smiles,
        inchi_key=None,
        molecular_weight=descriptors["molecular_weight"],
        logp=descriptors["logp"],
        tpsa=descriptors["tpsa"],
        structure_svg=structure_svg,
        source="unknown",
    )


def search_compound(query: str) -> CompoundSearchResponse:
    normalized_query = query.strip()
    is_smiles = _looks_like_smiles(normalized_query)

    pubchem_result: CompoundSearchResponse | None = None

    if not is_smiles:
        curated_cid = resolve_curated_cid(normalized_query)
        if curated_cid is not None:
            try:
                pubchem_result = fetch_pubchem_compound_by_cid(curated_cid, normalized_query)
            except httpx.HTTPError:
                pubchem_result = None

    if pubchem_result is None:
        try:
            pubchem_result = (
                search_pubchem_by_smiles(normalized_query)
                if is_smiles
                else search_pubchem_by_name(normalized_query)
            )
        except httpx.HTTPError:
            pubchem_result = None

    if pubchem_result is not None:
        display_name = display_name_for_cid(pubchem_result.cid, pubchem_result.name)
        if display_name != pubchem_result.name:
            return pubchem_result.model_copy(update={"name": display_name})
        return pubchem_result

    if is_smiles:
        if _rdkit_mol(normalized_query) is None:
            raise LookupError(f"Compound '{query}' could not be parsed as a valid structure.")
        return _build_manual_smiles_response(normalized_query)

    normalized = normalized_query.lower()
    compound = MOCK_COMPOUNDS.get(normalized)

    if compound is None:
        raise LookupError(f"Compound '{query}' was not found in PubChem.")

    return CompoundSearchResponse(
        **compound,
        structure_svg=get_demo_svg(compound["name"]),
        source="mock",
    )


def render_compound_png(query: str) -> bytes:
    normalized_query = query.strip()
    is_smiles = _looks_like_smiles(normalized_query)

    if is_smiles:
        smiles_png = smiles_to_png_bytes(normalized_query)
        if smiles_png is not None:
            return smiles_png

    try:
        if is_smiles:
            pubchem_result = search_pubchem_by_smiles(normalized_query)
            if pubchem_result is not None:
                canonical_png = smiles_to_png_bytes(pubchem_result.smiles)
                if canonical_png is not None:
                    return canonical_png
            png_bytes = fetch_pubchem_png_bytes_by_smiles(normalized_query)
        else:
            pubchem_result = search_pubchem_by_name(normalized_query)
            if pubchem_result is not None:
                canonical_png = smiles_to_png_bytes(pubchem_result.smiles)
                if canonical_png is not None:
                    return canonical_png
            png_bytes = fetch_pubchem_png_bytes_by_name(normalized_query)
    except httpx.HTTPError:
        png_bytes = None

    if png_bytes is not None:
        return png_bytes

    compound = MOCK_COMPOUNDS.get(normalized_query.lower())
    if compound is not None:
        mock_png = smiles_to_png_bytes(compound["smiles"])
        if mock_png is not None:
            return mock_png

    if is_smiles:
        raise LookupError(f"Compound '{query}' could not be rendered.")

    raise LookupError(f"Compound '{query}' was not found in PubChem.")
