import os

from fastapi import APIRouter, File, HTTPException, Query, UploadFile
from fastapi.responses import Response

from app.schemas import (
    CompoundSearchResponse,
    DockingRequest,
    DockingResponse,
    StructureFileParseResponse,
)
from app.services.compound_service import (
    render_compound_png,
    search_compound as search_compound_service,
    unknown_compound_response,
)
from app.services.structure_file import parse_structure_file

_DOCKING_ENABLED = os.getenv("ENABLE_DOCKING", "true").lower() == "true"

router = APIRouter(prefix="/compound", tags=["compound"])


@router.get("/search", response_model=CompoundSearchResponse)
def search_compound(name: str = Query(..., min_length=1)) -> CompoundSearchResponse:
    try:
        return search_compound_service(name)
    except LookupError:
        return unknown_compound_response(name)


@router.post("/parse-structure-file", response_model=StructureFileParseResponse)
async def parse_structure_file_upload(
    file: UploadFile = File(...),
) -> StructureFileParseResponse:
    file_bytes = await file.read()
    if not file_bytes:
        raise HTTPException(status_code=400, detail="Uploaded file was empty.")
    return parse_structure_file(file.filename or "structure", file_bytes)


@router.post("/dock", response_model=DockingResponse)
def dock_compound(body: DockingRequest) -> DockingResponse:
    """Run AutoDock Vina against the MRGPRX2 orthosteric pocket (7VDH).

    Returns the best binding affinity (kcal/mol; more negative = stronger).
    Typical range for binders: −5 to −10 kcal/mol.
    Disabled when ENABLE_DOCKING=false (e.g. memory-constrained deployments).
    """
    if not _DOCKING_ENABLED:
        raise HTTPException(
            status_code=503,
            detail=(
                "DOCKING_UNAVAILABLE: AutoDock Vina is disabled on this server "
                "due to memory constraints (free-tier deployment). "
                "Please enter the docking score manually, or run the backend locally."
            ),
        )
    from app.services.docking_service import run_docking  # lazy import — saves ~150 MB on startup
    try:
        result = run_docking(body.smiles, exhaustiveness=body.exhaustiveness)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except RuntimeError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    return DockingResponse(**result)


@router.get("/render")
def render_compound(name: str = Query(..., min_length=1)) -> Response:
    try:
        image_bytes = render_compound_png(name)
        return Response(content=image_bytes, media_type="image/png")
    except LookupError as error:
        raise HTTPException(status_code=404, detail=str(error)) from error
