from fastapi import APIRouter, HTTPException

from app.schemas import PredictRequest, PredictResponse
from app.services.predict_service import predict_activity

router = APIRouter(tags=["predict"])


@router.post("/predict", response_model=PredictResponse)
def predict(payload: PredictRequest) -> PredictResponse:
    try:
        return predict_activity(payload)
    except LookupError as error:
        raise HTTPException(status_code=404, detail=str(error)) from error
