from fastapi import APIRouter, File, HTTPException, Query, UploadFile
from fastapi.responses import Response

from app.schemas import (
    CompoundSearchResponse,
    StructureFileParseResponse,
)
from app.services.compound_service import (
    render_compound_png,
    search_compound as search_compound_service,
    unknown_compound_response,
)
from app.services.structure_file import parse_structure_file

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


@router.get("/render")
def render_compound(name: str = Query(..., min_length=1)) -> Response:
    try:
        image_bytes = render_compound_png(name)
        return Response(content=image_bytes, media_type="image/png")
    except LookupError as error:
        raise HTTPException(status_code=404, detail=str(error)) from error
