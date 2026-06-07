import base64
from urllib.parse import quote

import httpx

from app.schemas import CompoundSearchResponse
from app.services.molfile_renderer import molfile_to_highlight_svg
from app.services.rdkit_renderer import smiles_to_png_markup, smiles_to_svg

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_TIMEOUT = 15.0
PROPERTY_FIELDS = [
    "Title",
    "SMILES",
    "ConnectivitySMILES",
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
    return f"data:image/png;base64,{base64.b64encode(_fetch_png_bytes(cid)).decode('ascii')}"


def _fetch_png_bytes(cid: int) -> bytes:
    response = httpx.get(
        f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/PNG",
        params={"image_size": "300x300"},
        timeout=PUBCHEM_TIMEOUT,
    )
    response.raise_for_status()
    return response.content


def _fetch_sdf(cid: int) -> str:
    response = httpx.get(
        f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/SDF",
        timeout=PUBCHEM_TIMEOUT,
    )
    response.raise_for_status()
    return response.text


def _fetch_structure_svg(cid: int) -> str:
    try:
        molfile = _fetch_sdf(cid)
        return molfile_to_highlight_svg(molfile)
    except (httpx.HTTPError, ValueError):
        return _fetch_png_data_uri(cid)


def render_smiles_to_svg(smiles: str) -> str | None:
    encoded_smiles = quote(smiles.strip(), safe="")
    response = httpx.get(
        f"{PUBCHEM_BASE_URL}/standardize/smiles/{encoded_smiles}/SDF",
        timeout=PUBCHEM_TIMEOUT,
    )
    if response.status_code >= 400:
        return None

    try:
        return molfile_to_highlight_svg(response.text)
    except ValueError:
        return None


def _build_compound_response(cid: int, fallback_name: str) -> CompoundSearchResponse | None:
    property_payload = _fetch_json(
        f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/{','.join(PROPERTY_FIELDS)}/JSON"
    )
    properties = property_payload.get("PropertyTable", {}).get("Properties", [])

    if not properties:
        return None

    record = properties[0]
    canonical_smiles = (
        record.get("SMILES")
        or record.get("CanonicalSMILES")
        or record.get("ConnectivitySMILES")
        or ""
    )
    structure_svg = (
        smiles_to_svg(canonical_smiles)
        or smiles_to_png_markup(canonical_smiles)
        or _fetch_structure_svg(cid)
    )
    return CompoundSearchResponse(
        name=record.get("Title") or fallback_name.strip(),
        cid=cid,
        smiles=canonical_smiles,
        inchi_key=record.get("InChIKey"),
        molecular_weight=float(record.get("MolecularWeight") or 0.0),
        logp=float(record.get("XLogP") or 0.0),
        tpsa=float(record.get("TPSA") or 0.0),
        structure_svg=structure_svg,
        source="pubchem",
    )


def fetch_pubchem_compound_by_cid(cid: int, fallback_name: str) -> CompoundSearchResponse | None:
    return _build_compound_response(cid, fallback_name)


def search_pubchem_by_name(name: str) -> CompoundSearchResponse | None:
    encoded_name = quote(name.strip(), safe="")
    cid_payload = _fetch_json(
        f"{PUBCHEM_BASE_URL}/compound/name/{encoded_name}/cids/JSON"
    )
    cids = cid_payload.get("IdentifierList", {}).get("CID", [])

    if not cids:
        return None

    cid = cids[0]
    return _build_compound_response(cid, name)


def search_pubchem_by_smiles(smiles: str) -> CompoundSearchResponse | None:
    encoded_smiles = quote(smiles.strip(), safe="")
    cid_payload = _fetch_json(
        f"{PUBCHEM_BASE_URL}/compound/smiles/{encoded_smiles}/cids/JSON"
    )
    cids = cid_payload.get("IdentifierList", {}).get("CID", [])

    if not cids:
        return None

    cid = cids[0]
    return _build_compound_response(cid, smiles)


def fetch_pubchem_png_bytes_by_name(name: str) -> bytes | None:
    encoded_name = quote(name.strip(), safe="")
    cid_payload = _fetch_json(
        f"{PUBCHEM_BASE_URL}/compound/name/{encoded_name}/cids/JSON"
    )
    cids = cid_payload.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return None
    return _fetch_png_bytes(cids[0])


def fetch_pubchem_png_bytes_by_smiles(smiles: str) -> bytes | None:
    encoded_smiles = quote(smiles.strip(), safe="")
    cid_payload = _fetch_json(
        f"{PUBCHEM_BASE_URL}/compound/smiles/{encoded_smiles}/cids/JSON"
    )
    cids = cid_payload.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return None
    return _fetch_png_bytes(cids[0])
