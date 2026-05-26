"use client";

import { useState } from "react";
import { searchCompound, predictCompound } from "@/lib/api";
import type {
  CompoundSearchResponse,
  PredictionResponse,
} from "@/lib/types";
import { ProbabilityChart } from "@/components/probability-chart";

const EXAMPLES = ["Ciprofloxacin", "Morphine", "Vancomycin", "Substance P"];

export function SearchWorkspace() {
  const [query, setQuery] = useState("Ciprofloxacin");
  const [compound, setCompound] = useState<CompoundSearchResponse | null>(null);
  const [prediction, setPrediction] = useState<PredictionResponse | null>(null);
  const [loadingSearch, setLoadingSearch] = useState(false);
  const [loadingPredict, setLoadingPredict] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function onSearch() {
    setLoadingSearch(true);
    setError(null);
    setPrediction(null);

    try {
      const result = await searchCompound(query);
      setCompound(result);
    } catch (searchError) {
      setError(
        searchError instanceof Error
          ? searchError.message
          : "Compound search failed",
      );
    } finally {
      setLoadingSearch(false);
    }
  }

  async function onPredict() {
    setLoadingPredict(true);
    setError(null);

    try {
      const result = await predictCompound(query);
      setPrediction(result);
    } catch (predictError) {
      setError(
        predictError instanceof Error
          ? predictError.message
          : "Prediction failed",
      );
    } finally {
      setLoadingPredict(false);
    }
  }

  return (
    <main className="mx-auto flex min-h-screen w-full max-w-6xl flex-col gap-8 px-6 py-10 md:px-10">
      <section className="rounded-[32px] border border-black/10 bg-panel/90 p-8 shadow-[0_24px_80px_rgba(16,35,28,0.08)]">
        <p className="mb-3 text-sm font-medium uppercase tracking-[0.28em] text-accent">
          Structure-first workspace
        </p>
        <div className="flex flex-col gap-4 md:flex-row md:items-end md:justify-between">
          <div className="max-w-2xl">
            <h1 className="text-4xl font-semibold tracking-tight">
              MRGPRX2 Ligand Predictor
            </h1>
            <p className="mt-3 text-base leading-7 text-ink/75">
              Chemical name or SMILES input, PubChem lookup, molecular summary,
              and a first prediction panel for agonist-like, antagonist-like,
              or nonbinder-like activity.
            </p>
          </div>
          <div className="rounded-full border border-accent/20 bg-accent/5 px-4 py-2 text-sm text-accent">
            Frontend MVP
          </div>
        </div>

        <div className="mt-8 grid gap-3 md:grid-cols-[1fr_auto_auto]">
          <input
            className="h-14 rounded-2xl border border-black/10 bg-white px-4 text-base outline-none transition focus:border-accent"
            value={query}
            onChange={(event) => setQuery(event.target.value)}
            placeholder="Enter chemical name or SMILES"
          />
          <button
            className="h-14 rounded-2xl bg-accent px-6 text-base font-medium text-white transition hover:bg-accent/90 disabled:cursor-not-allowed disabled:opacity-60"
            onClick={onSearch}
            disabled={loadingSearch}
          >
            {loadingSearch ? "Searching..." : "Search"}
          </button>
          <button
            className="h-14 rounded-2xl border border-black/10 bg-white px-6 text-base font-medium text-ink transition hover:border-signal hover:text-signal disabled:cursor-not-allowed disabled:opacity-60"
            onClick={onPredict}
            disabled={loadingPredict}
          >
            {loadingPredict ? "Predicting..." : "Predict"}
          </button>
        </div>

        <div className="mt-4 flex flex-wrap gap-2 text-sm text-ink/70">
          {EXAMPLES.map((example) => (
            <button
              key={example}
              className="rounded-full border border-black/10 px-3 py-1.5 transition hover:border-accent hover:text-accent"
              onClick={() => setQuery(example)}
            >
              {example}
            </button>
          ))}
        </div>

        {error ? (
          <div className="mt-6 rounded-2xl border border-red-200 bg-red-50 px-4 py-3 text-sm text-red-700">
            {error}
          </div>
        ) : null}
      </section>

      <section className="grid gap-6 lg:grid-cols-[1.15fr_0.85fr]">
        <article className="rounded-[28px] border border-black/10 bg-white/85 p-6 shadow-[0_20px_60px_rgba(16,35,28,0.06)]">
          <div className="flex items-center justify-between">
            <h2 className="text-xl font-semibold">Compound Summary</h2>
            <span className="text-sm text-ink/60">
              {compound?.source === "pubchem" ? "PubChem" : "Demo data"}
            </span>
          </div>

          {compound ? (
            <div className="mt-5 grid gap-6 md:grid-cols-[220px_1fr]">
              <div
                className="flex min-h-[220px] items-center justify-center rounded-3xl border border-black/10 bg-canvas p-4"
                dangerouslySetInnerHTML={{ __html: compound.structure_svg }}
              />
              <div className="grid gap-3 text-sm">
                <InfoRow label="Name" value={compound.name} />
                <InfoRow label="PubChem CID" value={compound.cid?.toString() ?? "-"} />
                <InfoRow label="SMILES" value={compound.smiles} mono />
                <InfoRow label="InChIKey" value={compound.inchi_key ?? "-"} mono />
                <InfoRow
                  label="Molecular Weight"
                  value={compound.molecular_weight.toFixed(2)}
                />
                <InfoRow label="LogP" value={compound.logp.toFixed(2)} />
                <InfoRow label="TPSA" value={compound.tpsa.toFixed(2)} />
              </div>
            </div>
          ) : (
            <EmptyState message="Search a compound to load its basic structure and descriptors." />
          )}
        </article>

        <article className="rounded-[28px] border border-black/10 bg-white/85 p-6 shadow-[0_20px_60px_rgba(16,35,28,0.06)]">
          <div className="flex items-center justify-between">
            <h2 className="text-xl font-semibold">Prediction Panel</h2>
            <span className="text-sm text-ink/60">MRGPRX2</span>
          </div>

          {prediction ? (
            <div className="mt-5 flex flex-col gap-5">
              <div className="rounded-3xl bg-ink px-5 py-4 text-white">
                <p className="text-xs uppercase tracking-[0.28em] text-white/60">
                  Final call
                </p>
                <p className="mt-2 text-3xl font-semibold">
                  {prediction.prediction}
                </p>
                <p className="mt-2 text-sm text-white/75">
                  {prediction.interpretation}
                </p>
              </div>

              <ProbabilityChart probabilities={prediction.probabilities} />

              <div className="grid gap-3 text-sm">
                <InfoRow
                  label="Docking Score"
                  value={
                    prediction.docking_score !== null
                      ? `${prediction.docking_score.toFixed(2)} kcal/mol`
                      : "Not run"
                  }
                />
                <InfoRow
                  label="Model Version"
                  value={prediction.model_version}
                />
                <InfoRow
                  label="Descriptor Snapshot"
                  value={`MW ${prediction.descriptors.molecular_weight.toFixed(2)} / LogP ${prediction.descriptors.logp.toFixed(2)} / TPSA ${prediction.descriptors.tpsa.toFixed(2)}`}
                />
              </div>
            </div>
          ) : (
            <EmptyState message="Run prediction after entering a compound name. Current backend returns a mock baseline result." />
          )}
        </article>
      </section>
    </main>
  );
}

function InfoRow({
  label,
  value,
  mono = false,
}: {
  label: string;
  value: string;
  mono?: boolean;
}) {
  return (
    <div className="grid gap-1 rounded-2xl border border-black/5 bg-canvas/60 px-4 py-3">
      <span className="text-xs uppercase tracking-[0.18em] text-ink/45">
        {label}
      </span>
      <span className={mono ? "font-mono text-[13px]" : "text-[15px]"}>
        {value}
      </span>
    </div>
  );
}

function EmptyState({ message }: { message: string }) {
  return (
    <div className="mt-5 rounded-3xl border border-dashed border-black/15 bg-canvas/40 px-5 py-10 text-sm text-ink/60">
      {message}
    </div>
  );
}
