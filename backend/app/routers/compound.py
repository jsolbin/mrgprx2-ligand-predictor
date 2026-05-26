from fastapi import APIRouter, HTTPException, Query

from app.schemas import CompoundSearchResponse
from app.services.compound_service import search_compound_by_name

router = APIRouter(prefix="/compound", tags=["compound"])


@router.get("/search", response_model=CompoundSearchResponse)
def search_compound(name: str = Query(..., min_length=1)) -> CompoundSearchResponse:
    try:
        return search_compound_by_name(name)
    except LookupError as error:
        raise HTTPException(status_code=404, detail=str(error)) from error
