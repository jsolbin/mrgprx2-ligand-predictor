import {
  CompoundSearchResponse,
  LabelSubmissionResponse,
  PredictionResponse,
  StructureFileParseResponse,
} from "@/lib/types";

const API_BASE_URL =
  process.env.NEXT_PUBLIC_API_BASE_URL ?? "http://localhost:8000";

async function handleResponse<T>(response: Response): Promise<T> {
  if (!response.ok) {
    const text = await response.text();
    throw new Error(text || "Request failed");
  }

  return response.json() as Promise<T>;
}

export async function searchCompound(query: string) {
  const url = new URL("/compound/search", API_BASE_URL);
  url.searchParams.set("name", query);

  const response = await fetch(url.toString(), {
    method: "GET",
    cache: "no-store",
  });

  return handleResponse<CompoundSearchResponse>(response);
}

type PredictPayload = {
  compound_name: string;
  receptor: string;
  experimental_data: {
    docking_score: number | null;
    mrna_fold_change: number | null;
    mrna_method: string | null;
    protein_fold_change: number | null;
    protein_method: string | null;
    cell_line: string | null;
    concentration: number | null;
    time_hours: number | null;
  };
  binding_mode: string;
  selected_residues: string[];
};

export async function predictCompound(payload: PredictPayload) {
  const response = await fetch(new URL("/predict", API_BASE_URL), {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify(payload),
  });

  return handleResponse<PredictionResponse>(response);
}

type LabelPayload = {
  smiles: string;
  label: "agonist" | "antagonist" | "nonbinder";
  name?: string | null;
  note?: string | null;
};

export async function labelCompound(payload: LabelPayload) {
  const response = await fetch(new URL("/training/label", API_BASE_URL), {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify(payload),
  });

  return handleResponse<LabelSubmissionResponse>(response);
}

export async function parseStructureFile(file: File) {
  const body = new FormData();
  body.append("file", file);

  const response = await fetch(new URL("/compound/parse-structure-file", API_BASE_URL), {
    method: "POST",
    body,
  });

  return handleResponse<StructureFileParseResponse>(response);
}

export function getCompoundRenderUrl(query: string) {
  const url = new URL("/compound/render", API_BASE_URL);
  url.searchParams.set("name", query);
  return url.toString();
}
