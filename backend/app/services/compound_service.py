import httpx

from app.schemas import CompoundSearchResponse
from app.services.mock_data import MOCK_COMPOUNDS, get_demo_svg
from app.services.pubchem_service import search_pubchem_by_name


def search_compound_by_name(name: str) -> CompoundSearchResponse:
    try:
        pubchem_result = search_pubchem_by_name(name)
    except httpx.HTTPError:
        pubchem_result = None

    if pubchem_result is not None:
        return pubchem_result

    normalized = name.strip().lower()
    compound = MOCK_COMPOUNDS.get(normalized)

    if compound is None:
        raise LookupError(f"Compound '{name}' was not found in PubChem.")

    return CompoundSearchResponse(
        **compound,
        structure_svg=get_demo_svg(compound["name"]),
        source="mock",
    )
