from pydantic import BaseModel, Field


class DockingRequest(BaseModel):
    smiles: str = Field(..., min_length=1)
    exhaustiveness: int = Field(default=8, ge=1, le=32)


class DockingResponse(BaseModel):
    affinity_kcal_mol: float
    num_modes: int
    warning: str | None = None


class CompoundSearchResponse(BaseModel):
    name: str
    cid: int | None
    smiles: str
    inchi_key: str | None
    molecular_weight: float
    logp: float
    tpsa: float
    structure_svg: str | None
    source: str


class ExperimentalDataInput(BaseModel):
    docking_score: float | None = None
    mrna_fold_change: float | None = None
    mrna_method: str | None = None
    protein_fold_change: float | None = None
    protein_method: str | None = None
    cell_line: str | None = None
    concentration: float | None = None
    time_hours: float | None = None


class PredictRequest(BaseModel):
    compound_name: str = Field(..., min_length=1)
    receptor: str = Field(default="MRGPRX2", min_length=1)
    experimental_data: ExperimentalDataInput | None = None
    binding_mode: str | None = None
    selected_residues: list[str] = Field(default_factory=list)


class PredictionProbabilities(BaseModel):
    agonist: float
    antagonist: float
    nonbinder: float


class DescriptorSnapshot(BaseModel):
    molecular_weight: float
    logp: float
    tpsa: float
    h_bond_donors: int
    h_bond_acceptors: int
    rotatable_bonds: int
    ring_count: int


class FactorScore(BaseModel):
    label: str
    weight: float
    score: float


class FeatureGroup(BaseModel):
    title: str
    items: list[str]


class SimilarLigand(BaseModel):
    name: str
    label: str
    similarity: float
    rationale: str
    source: str = "curated"


class ReceptorResidueSet(BaseModel):
    mrgprx1: list[str]
    mrgprx2: list[str]
    note: str


class ReceptorPocketFeature(BaseModel):
    residues: list[str]
    feature: str
    structural_signature: str


class Mrgprx1Comparison(BaseModel):
    mrgprx2_selectivity_score: float
    closest_mrgprx1_ligand: str | None
    mrgprx1_similarity: float
    shared_residues: ReceptorResidueSet
    mrgprx1_pocket: ReceptorPocketFeature
    mrgprx2_pocket: ReceptorPocketFeature
    signals: list[str]
    summary: str


class WeightedFactor(BaseModel):
    label: str
    weight: float
    score: float
    contribution_pct: float


class LabelSubmissionRequest(BaseModel):
    smiles: str = Field(..., min_length=1)
    label: str = Field(..., min_length=1)
    name: str | None = None
    note: str | None = None


class LabelSubmissionResponse(BaseModel):
    name: str
    smiles: str
    label: str
    note: str | None
    source: str
    submitted_at: str
    total_labeled_compounds: int


class LabeledCompoundSummary(BaseModel):
    name: str
    smiles: str
    label: str
    note: str | None
    source: str
    submitted_at: str


class DrugLikenessCheck(BaseModel):
    label: str
    passed: bool


class BindingAnalysis(BaseModel):
    agonist_likelihood: float
    antagonist_likelihood: float


class ModelNeighbor(BaseModel):
    name: str
    label: str
    similarity: float
    source: str


class ModelFeatureContribution(BaseModel):
    description: str
    importance: float


class ExperimentalAdjustmentComponent(BaseModel):
    description: str
    delta_pct: float


class ExperimentalAdjustment(BaseModel):
    target_label: str
    structure_probability: float
    adjusted_probability: float
    applied_delta_pct: float
    components: list[ExperimentalAdjustmentComponent]
    specificity_note: str | None = None


class ModelEvidence(BaseModel):
    model_type: str
    trained_on: int
    nonbinder_probability: float
    agonist_probability: float
    antagonist_probability: float
    nearest_neighbors: list[ModelNeighbor]
    top_features: list[ModelFeatureContribution]
    summary: str
    experimental_adjustment: ExperimentalAdjustment | None = None


class ApplicabilityDomain(BaseModel):
    in_domain: bool
    reason: str

class AssayBasis(BaseModel):
    readout: str
    note: str

class ReceptorRegulation(BaseModel):
    mrna_note: str | None = None
    protein_note: str | None = None
    warning: str


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
    confidence: float
    factor_analysis: list[FactorScore]
    weighted_factors: list[WeightedFactor]
    composite_score: float
    binding_analysis: BindingAnalysis
    feature_groups: list[FeatureGroup]
    similar_ligand: SimilarLigand
    mrgprx1_comparison: Mrgprx1Comparison
    drug_likeness: list[DrugLikenessCheck]
    model_evidence: ModelEvidence | None = None
    analyzed_at: str
    applicability_domain: ApplicabilityDomain | None = None
    assay_basis: AssayBasis | None = None
    receptor_regulation: ReceptorRegulation | None = None


class StructureFileParseResponse(BaseModel):
    name: str | None
    smiles: str | None
    message: str


