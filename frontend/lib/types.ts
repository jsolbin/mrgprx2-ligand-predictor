export type CompoundSearchResponse = {
  name: string;
  cid: number | null;
  smiles: string;
  inchi_key: string | null;
  molecular_weight: number;
  logp: number;
  tpsa: number;
  structure_svg: string;
  source: "mock" | "pubchem";
};

export type PredictionResponse = {
  compound_name: string;
  receptor: string;
  smiles: string;
  prediction: "agonist-like" | "antagonist-like" | "nonbinder-like";
  probabilities: {
    agonist: number;
    antagonist: number;
    nonbinder: number;
  };
  descriptors: {
    molecular_weight: number;
    logp: number;
    tpsa: number;
  };
  docking_score: number | null;
  interpretation: string;
  model_version: string;
};
