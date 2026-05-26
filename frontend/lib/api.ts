import { CompoundSearchResponse, PredictionResponse } from "@/lib/types";

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

export async function predictCompound(query: string) {
  const response = await fetch(new URL("/predict", API_BASE_URL), {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      compound_name: query,
      receptor: "MRGPRX2",
    }),
  });

  return handleResponse<PredictionResponse>(response);
}
