from app.schemas import (
    DescriptorSnapshot,
    PredictRequest,
    PredictResponse,
    PredictionProbabilities,
)
from app.services.compound_service import search_compound_by_name


def predict_activity(payload: PredictRequest) -> PredictResponse:
    compound = search_compound_by_name(payload.compound_name)
    name = compound.name.lower()

    if "qwf" in name:
        probabilities = PredictionProbabilities(
            agonist=0.12, antagonist=0.79, nonbinder=0.09
        )
        prediction = "antagonist-like"
        docking_score = -8.3
        interpretation = (
            "Demo baseline assigns higher antagonist probability based on the seed list."
        )
    elif "substance p" in name or "ciprofloxacin" in name:
        probabilities = PredictionProbabilities(
            agonist=0.78, antagonist=0.14, nonbinder=0.08
        )
        prediction = "agonist-like"
        docking_score = -7.4
        interpretation = (
            "Demo baseline assigns higher agonist probability based on the seed list."
        )
    else:
        probabilities = PredictionProbabilities(
            agonist=0.33, antagonist=0.27, nonbinder=0.40
        )
        prediction = "nonbinder-like"
        docking_score = None
        interpretation = (
            "No trained model is attached yet, so this is a placeholder uncertainty output."
        )

    return PredictResponse(
        compound_name=compound.name,
        receptor=payload.receptor,
        smiles=compound.smiles,
        prediction=prediction,
        probabilities=probabilities,
        descriptors=DescriptorSnapshot(
            molecular_weight=compound.molecular_weight,
            logp=compound.logp,
            tpsa=compound.tpsa,
        ),
        docking_score=docking_score,
        interpretation=interpretation,
        model_version="baseline-mock-v1",
    )
