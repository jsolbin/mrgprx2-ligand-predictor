from fastapi import APIRouter, HTTPException

from app.schemas import (
    LabelSubmissionRequest,
    LabelSubmissionResponse,
    LabeledCompoundSummary,
)
from app.services.training_store import add_label, list_labels

router = APIRouter(prefix="/training", tags=["training"])


@router.post("/label", response_model=LabelSubmissionResponse)
def label_compound(payload: LabelSubmissionRequest) -> LabelSubmissionResponse:
    try:
        entry = add_label(
            smiles=payload.smiles,
            label=payload.label,
            name=payload.name,
            note=payload.note,
        )
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error

    return LabelSubmissionResponse(
        **entry,
        total_labeled_compounds=len(list_labels()),
    )


@router.get("/labels", response_model=list[LabeledCompoundSummary])
def get_labels() -> list[LabeledCompoundSummary]:
    return [LabeledCompoundSummary(**entry) for entry in list_labels()]
