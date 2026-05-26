from pydantic import BaseModel, Field


class CompoundSearchResponse(BaseModel):
    name: str
    cid: int | None
    smiles: str
    inchi_key: str | None
    molecular_weight: float
    logp: float
    tpsa: float
    structure_svg: str
    source: str


class PredictRequest(BaseModel):
    compound_name: str = Field(..., min_length=1)
    receptor: str = Field(default="MRGPRX2", min_length=1)


class PredictionProbabilities(BaseModel):
    agonist: float
    antagonist: float
    nonbinder: float


class DescriptorSnapshot(BaseModel):
    molecular_weight: float
    logp: float
    tpsa: float


class PredictResponse(BaseModel):
    compound_name: str
    receptor: str
    smiles: str
    prediction: str
    probabilities: PredictionProbabilities
    descriptors: DescriptorSnapshot
    docking_score: float | None
    interpretation: str
    model_version: str
