from app.schemas import CompoundSearchResponse
from app.services.mock_data import MOCK_COMPOUNDS, get_demo_svg


def search_compound_by_name(name: str) -> CompoundSearchResponse:
    normalized = name.strip().lower()
    compound = MOCK_COMPOUNDS.get(normalized)

    if compound is None:
        compound = {
            "name": name.strip(),
            "cid": None,
            "smiles": "C",
            "inchi_key": None,
            "molecular_weight": 16.04,
            "logp": 0.64,
            "tpsa": 0.00,
        }

    return CompoundSearchResponse(
        **compound,
        structure_svg=get_demo_svg(compound["name"]),
        source="mock",
    )
