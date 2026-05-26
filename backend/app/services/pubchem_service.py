import base64
from urllib.parse import quote

import httpx

from app.schemas import CompoundSearchResponse

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_TIMEOUT = 15.0
PROPERTY_FIELDS = [
    "Title",
    "CanonicalSMILES",
    "InChIKey",
    "MolecularWeight",
    "XLogP",
    "TPSA",
]


def _fetch_json(url: str) -> dict:
    response = httpx.get(url, timeout=PUBCHEM_TIMEOUT)
    response.raise_for_status()
    return response.json()


def _fetch_png_data_uri(cid: int) -> str:
    response = httpx.get(
        f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/PNG",
        params={"image_size": "300x300"},
        timeout=PUBCHEM_TIMEOUT,
    )
    response.raise_for_status()
    encoded = base64.b64encode(response.content).decode("ascii")
    return f'<img src="data:image/png;base64,{encoded}" alt="PubChem structure" />'


def search_pubchem_by_name(name: str) -> CompoundSearchResponse | None:
    encoded_name = quote(name.strip(), safe="")
    cid_payload = _fetch_json(
        f"{PUBCHEM_BASE_URL}/compound/name/{encoded_name}/cids/JSON"
    )
    cids = cid_payload.get("IdentifierList", {}).get("CID", [])

    if not cids:
        return None

    cid = cids[0]
    property_payload = _fetch_json(
        f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/{','.join(PROPERTY_FIELDS)}/JSON"
    )
    properties = property_payload.get("PropertyTable", {}).get("Properties", [])

    if not properties:
        return None

    record = properties[0]
    return CompoundSearchResponse(
        name=record.get("Title") or name.strip(),
        cid=cid,
        smiles=record.get("CanonicalSMILES", ""),
        inchi_key=record.get("InChIKey"),
        molecular_weight=float(record.get("MolecularWeight") or 0.0),
        logp=float(record.get("XLogP") or 0.0),
        tpsa=float(record.get("TPSA") or 0.0),
        structure_svg=_fetch_png_data_uri(cid),
        source="pubchem",
    )
