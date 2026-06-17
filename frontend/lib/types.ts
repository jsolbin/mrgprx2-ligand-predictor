export type CompoundSearchResponse = {
  name: string;
  cid: number | null;
  smiles: string;
  inchi_key: string | null;
  molecular_weight: number;
  logp: number;
  tpsa: number;
  structure_svg: string | null;
  source: "mock" | "pubchem" | "input";
};

export type PredictionResponse = {
  compound_name: string;
  receptor: string;
  smiles: string;
  prediction: "agonist-like" | "antagonist-like" | "nonbinder-like" | "indeterminate";
  probabilities: {
    agonist: number;
    antagonist: number;
    nonbinder: number;
  };
  descriptors: {
    molecular_weight: number;
    logp: number;
    tpsa: number;
    h_bond_donors: number;
    h_bond_acceptors: number;
    rotatable_bonds: number;
    ring_count: number;
  };
  docking_score: number | null;
  interpretation: string;
  model_version: string;
  confidence: number;
  factor_analysis: Array<{
    label: string;
    weight: number;
    score: number;
  }>;
  weighted_factors: Array<{
    label: string;
    weight: number;
    score: number;
    contribution_pct: number;
  }>;
  composite_score: number;
  binding_analysis: {
    agonist_likelihood: number;
    antagonist_likelihood: number;
  };
  feature_groups: Array<{
    title: string;
    items: string[];
  }>;
  similar_ligand: {
    name: string;
    label: string;
    similarity: number;
    rationale: string;
    source: "curated" | "user";
  };
  mrgprx1_comparison: {
    mrgprx2_selectivity_score: number;
    closest_mrgprx1_ligand: string | null;
    mrgprx1_similarity: number;
    shared_residues: {
      mrgprx1: string[];
      mrgprx2: string[];
      note: string;
    };
    mrgprx1_pocket: {
      residues: string[];
      feature: string;
      structural_signature: string;
    };
    mrgprx2_pocket: {
      residues: string[];
      feature: string;
      structural_signature: string;
    };
    signals: string[];
    summary: string;
  };
  drug_likeness: Array<{
    label: string;
    passed: boolean;
  }>;
  model_evidence: {
    model_type: string;
    trained_on: number;
    nonbinder_probability: number;
    agonist_probability: number;
    antagonist_probability: number;
    nearest_neighbors: Array<{
      name: string;
      label: string;
      similarity: number;
      source: string;
    }>;
    top_features: Array<{
      description: string;
      importance: number;
    }>;
    summary: string;
    experimental_adjustment: {
      target_label: string;
      structure_probability: number;
      adjusted_probability: number;
      applied_delta_pct: number;
      components: Array<{ description: string; delta_pct: number }>;
      specificity_note: string | null;
    } | null;
  } | null;
  applicability_domain: {
    in_domain: boolean;
    reason: string;
  } | null;
  assay_basis: {
    readout: string;
    note: string;
  } | null;
  receptor_regulation: {
    mrna_note: string | null;
    protein_note: string | null;
    warning: string;
  } | null;
  analyzed_at: string;
};

export type LabelSubmissionResponse = {
  name: string;
  smiles: string;
  label: "agonist" | "antagonist" | "nonbinder";
  note: string | null;
  source: string;
  submitted_at: string;
  total_labeled_compounds: number;
};

export type StructureFileParseResponse = {
  name: string | null;
  smiles: string | null;
  message: string;
};

