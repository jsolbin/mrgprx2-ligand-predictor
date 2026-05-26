from fastapi import APIRouter

from app.schemas import PredictRequest, PredictResponse
from app.services.predict_service import predict_activity

router = APIRouter(tags=["predict"])


@router.post("/predict", response_model=PredictResponse)
def predict(payload: PredictRequest) -> PredictResponse:
    return predict_activity(payload)
