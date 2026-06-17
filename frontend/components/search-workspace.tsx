"use client";

import type { ChangeEvent, ReactNode } from "react";
import { useEffect, useRef, useState } from "react";
import { ProbabilityChart } from "@/components/probability-chart";
import {
  dockCompound,
  getCompoundRenderUrl,
  labelCompound,
  parseStructureFile,
  predictCompound,
  searchCompound,
} from "@/lib/api";
import type { DockingResponse } from "@/lib/api";
import type {
  CompoundSearchResponse,
  LabelSubmissionResponse,
  PredictionResponse,
} from "@/lib/types";

const QUICK_EXAMPLES = [
  {
    name: "Compound 48/80",
    smiles: "CNCCC1=CC(=C(C=C1)OC)CC2=CC(=CC(=C2OC)CC3=C(C=CC(=C3)CCNC)OC)CCNC",
  },
  {
    name: "Ciprofloxacin",
    smiles: "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
  },
  {
    name: "Levofloxacin",
    smiles: "C[C@H]1COC2=C3N1C=C(C(=O)C3=CC(=C2N4CCN(CC4)C)F)C(=O)O",
  },
  {
    name: "Atracurium",
    smiles:
      "C[N+]1(CCC2=CC(=C(C=C2C1CC3=CC(=C(C=C3)OC)OC)OC)OC)CCC(=O)OCCCCCOC(=O)CC[N+]4(CCC5=CC(=C(C=C5C4CC6=CC(=C(C=C6)OC)OC)OC)OC)C",
  },
  {
    name: "Rocuronium",
    smiles:
      "CC(=O)O[C@H]1[C@H](C[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(C[C@@H]([C@H](C4)O)N5CCOCC5)C)C)[N+]6(CCCC6)CC=C",
  },
  {
    name: "Morphine",
    smiles: "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=C[C@H]3[C@H]1C5)O",
  },
  {
    name: "QWF",
    smiles:
      "CC(C)(C)OC(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1ccccc1)C(N)=O",
  },
  {
    name: "Quercetin",
    smiles: "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
  },
  {
    name: "Luteolin",
    smiles: "C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O",
  },
  {
    name: "Osthole",
    smiles: "CC(=CCC1=C(C=CC2=C1OC(=O)C=C2)OC)C",
  },
  {
    name: "Substance P fragment",
    smiles: "CC(=O)NC(CCCCN)C(=O)NC(CC1=CC=CC=C1)C(=O)N",
  },
  {
    name: "Aspirin (nonbinder)",
    smiles: "CC(=O)Oc1ccccc1C(=O)O",
  },
  {
    name: "Cholesterol (nonbinder)",
    smiles: "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
  },
  {
    name: "(R)-Salbutamol [stress test ①]",
    smiles: "CC(C)(C)NC[C@@H](O)c1ccc(O)c(CO)c1",
  },
  {
    name: "(S)-Salbutamol [stress test ②]",
    smiles: "CC(C)(C)NC[C@H](O)c1ccc(O)c(CO)c1",
  },
] as const;

const BINDING_MODES = ["Unknown", "Allosteric", "Orthosteric"] as const;

// A muted, low-contrast palette in the app's existing blue/slate family - keeps
// the stacked breakdown bar readable without the clashing primary colors a
// default chart palette would introduce.
const FACTOR_COLORS = ["#4a84f6", "#F5C049", "#647FB5", "#A08E65"] as const;

const INITIAL_EXPERIMENTAL_DATA = {
  dockingScore: "",
  mrnaFoldChange: "",
  mrnaMethod: "",
  proteinFoldChange: "",
  proteinMethod: "",
  cellLine: "",
  concentration: "",
  timeHours: "",
};

const RESIDUE_OPTIONS = [
  { id: "asp184", name: "Asp184", strength: "critical", site: "orthosteric" },
  { id: "glu164", name: "Glu164", strength: "critical", site: "orthosteric" },
  { id: "tyr279", name: "Tyr279", strength: "important", site: "orthosteric" },
  { id: "phe170", name: "Phe170", strength: "important", site: "orthosteric" },
  { id: "his259", name: "His259", strength: "moderate", site: "allosteric" },
] as const;

type HistoryEntry = {
  label: string;
  smiles: string;
  query: string;
};

export function SearchWorkspace() {
  const [query, setQuery] = useState("");
  const [compound, setCompound] = useState<CompoundSearchResponse | null>(null);
  const [structureImageUrl, setStructureImageUrl] = useState<string | null>(
    null,
  );
  const [prediction, setPrediction] = useState<PredictionResponse | null>(null);
  const [loadingPredict, setLoadingPredict] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [experimentalDataOpen, setExperimentalDataOpen] = useState(false);
  const [bindingResiduesOpen, setBindingResiduesOpen] = useState(false);
  const [experimentalData, setExperimentalData] = useState(
    INITIAL_EXPERIMENTAL_DATA,
  );
  const [bindingMode, setBindingMode] =
    useState<(typeof BINDING_MODES)[number]>("Unknown");
  const [selectedResidues, setSelectedResidues] = useState<string[]>([]);
  const [analysisHistory, setAnalysisHistory] = useState<HistoryEntry[]>([]);
  const [labelSubmitting, setLabelSubmitting] = useState<
    "agonist" | "antagonist" | "nonbinder" | null
  >(null);
  const [labelResult, setLabelResult] =
    useState<LabelSubmissionResponse | null>(null);
  const [labelError, setLabelError] = useState<string | null>(null);
  const [editingLabel, setEditingLabel] = useState(false);
  // The backend re-canonicalizes SMILES with RDKit before storing/returning a
  // label, so `labelResult.smiles` won't text-match the PubChem-canonical
  // `compound.smiles` we display. Track the exact structure string the label
  // was submitted *for* so "already labeled" detection stays accurate.
  const [labeledForSmiles, setLabeledForSmiles] = useState<string | null>(null);
  const [fileProcessing, setFileProcessing] = useState(false);
  const [fileProcessingMessage, setFileProcessingMessage] = useState<
    string | null
  >(null);
  const fileInputRef = useRef<HTMLInputElement | null>(null);
  const [dockingResult, setDockingResult] = useState<DockingResponse | null>(null);
  const [dockingLoading, setDockingLoading] = useState(false);
  const [dockingError, setDockingError] = useState<string | null>(null);

  const canPredict = query.trim().length > 0;
  const activeSmiles = (compound?.smiles ?? query).trim();

  // The "Saved as agonist" confirmation (and any in-progress label edit) is
  // only meaningful for the structure it was recorded against - once the
  // user moves on to a different structure, stale confirmation text from the
  // previous one shouldn't keep showing.
  useEffect(() => {
    setLabelResult((current) =>
      labeledForSmiles === activeSmiles ? current : null,
    );
    setLabeledForSmiles((current) =>
      current === activeSmiles ? current : null,
    );
    setLabelError(null);
    setEditingLabel(false);
    // labeledForSmiles intentionally excluded - it's the value we compare
    // against, not a trigger; including it would clear itself on every set.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeSmiles]);

  const isCurrentStructureLabeled =
    labelResult !== null && labeledForSmiles === activeSmiles;

  const displayedProbabilities = prediction?.probabilities ?? {
    agonist: 1,
    antagonist: 0,
    nonbinder: 0,
  };
  const displayedLabel = prediction?.prediction ?? "agonist-like";
  const displayedConfidence = prediction?.confidence ?? 0;
  const displayedCompoundName = compound?.name ?? "";
  const compositeScore = prediction?.composite_score ?? 0;
  const factorAnalysis = prediction?.factor_analysis ?? [];
  const weightedFactors = prediction?.weighted_factors ?? [];
  const featureGroups = prediction?.feature_groups ?? [];
  const mrgprx1Comparison = prediction?.mrgprx1_comparison ?? null;
  const modelEvidence = prediction?.model_evidence ?? null;
  const molecularProperties = prediction
    ? [
        {
          label: "Molecular Weight",
          value: `${prediction.descriptors.molecular_weight.toFixed(1)} g/mol`,
        },
        { label: "LogP", value: prediction.descriptors.logp.toFixed(2) },
        {
          label: "H-Bond Donors",
          value: String(prediction.descriptors.h_bond_donors),
        },
        {
          label: "H-Bond Acceptors",
          value: String(prediction.descriptors.h_bond_acceptors),
        },
        {
          label: "Rotatable Bonds",
          value: String(prediction.descriptors.rotatable_bonds),
        },
        {
          label: "TPSA",
          value: `${prediction.descriptors.tpsa.toFixed(1)} A²`,
        },
        {
          label: "Ring Count",
          value: String(prediction.descriptors.ring_count),
        },
      ]
    : [];

  async function onPredict(
    overrideQuery?: string,
    overrideDockingScore?: number | null,
    overrideDeltaDeltaScore?: number | null,
    overrideInactiveDockingScore?: number | null,
  ) {
    const targetQuery = (overrideQuery ?? query).trim();
    if (!targetQuery) {
      return;
    }

    if (overrideQuery !== undefined) {
      setQuery(overrideQuery);
      // Switching to a different compound — clear the previous docking cache
      // so stale scores from the old structure don't carry over.
      setDockingResult(null);
      setDockingError(null);
      updateExperimentalData("dockingScore", "");
    }

    setLoadingPredict(true);
    setError(null);
    setPrediction(null);
    setCompound(null);
    setStructureImageUrl(null);

    try {
      const [resultCompound, result] = await Promise.all([
        searchCompound(targetQuery),
        predictCompound({
          compound_name: targetQuery,
          receptor: "MRGPRX2",
          experimental_data: {
            docking_score: overrideDockingScore !== undefined
              ? overrideDockingScore
              : parseNullableNumber(experimentalData.dockingScore),
            delta_delta_score: overrideDeltaDeltaScore !== undefined
              ? overrideDeltaDeltaScore
              : null,
            inactive_docking_score: overrideInactiveDockingScore !== undefined
              ? overrideInactiveDockingScore
              : null,
            mrna_fold_change: parseNullableNumber(
              experimentalData.mrnaFoldChange,
            ),
            mrna_method: experimentalData.mrnaMethod || null,
            protein_fold_change: parseNullableNumber(
              experimentalData.proteinFoldChange,
            ),
            protein_method: experimentalData.proteinMethod || null,
            cell_line: experimentalData.cellLine || null,
            concentration: parseNullableNumber(experimentalData.concentration),
            time_hours: parseNullableNumber(experimentalData.timeHours),
          },
          binding_mode: bindingMode,
          selected_residues: selectedResidues,
        }),
      ]);

      // Reveal the structure and the full result set in one frame, so the
      // diagram never appears before its analysis is ready.
      setCompound(resultCompound);
      setStructureImageUrl(getCompoundRenderUrl(targetQuery));
      setPrediction(result);
      setAnalysisHistory((current) => [
        {
          label: `${formatPredictionLabel(result.prediction)} (${result.confidence.toFixed(1)}%)`,
          smiles: result.smiles,
          query: targetQuery,
        },
        ...current.filter((item) => item.smiles !== result.smiles).slice(0, 4),
      ]);
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

  async function onLabel(label: "agonist" | "antagonist" | "nonbinder") {
    const smiles = (compound?.smiles ?? query).trim();
    if (!smiles) {
      return;
    }

    setLabelSubmitting(label);
    setLabelError(null);
    setLabelResult(null);

    try {
      const result = await labelCompound({
        smiles,
        label,
        name: compound?.name ?? null,
      });
      setLabelResult(result);
      setLabeledForSmiles(smiles);
      setEditingLabel(false);
    } catch (labelSubmitError) {
      setLabelError(
        labelSubmitError instanceof Error
          ? labelSubmitError.message
          : "Saving the label failed",
      );
    } finally {
      setLabelSubmitting(null);
    }
  }

  function updateExperimentalData(
    key: keyof typeof INITIAL_EXPERIMENTAL_DATA,
    value: string,
  ) {
    setExperimentalData((current) => ({
      ...current,
      [key]: value,
    }));
  }

  function toggleResidue(id: string) {
    setSelectedResidues((current) =>
      current.includes(id)
        ? current.filter((residueId) => residueId !== id)
        : [...current, id],
    );
  }

  function onUploadClick() {
    fileInputRef.current?.click();
  }

  async function onFileChange(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    event.target.value = "";
    if (!file) {
      return;
    }

    setFileProcessing(true);
    setFileProcessingMessage(null);
    try {
      const result = await parseStructureFile(file);
      setFileProcessingMessage(result.message);
      if (result.smiles) {
        setQuery(result.smiles);
      }
    } catch (caughtError) {
      setFileProcessingMessage(
        caughtError instanceof Error
          ? caughtError.message
          : "Couldn't read that file. Please upload a valid SDF/MOL2 file or paste the SMILES directly.",
      );
    } finally {
      setFileProcessing(false);
    }
  }

  function onSelectExample(example: (typeof QUICK_EXAMPLES)[number]) {
    setQuery(example.smiles);
    setCompound(null);
    setStructureImageUrl(null);
    setPrediction(null);
    setError(null);
    setExperimentalData(INITIAL_EXPERIMENTAL_DATA);
    setBindingMode("Unknown");
    setSelectedResidues([]);
    setDockingResult(null);
    setDockingError(null);
  }

  return (
    <main className="h-screen overflow-hidden bg-[#f7f8fc] text-[#171a22]">
      <header className="flex items-center gap-5 border-b border-[#dfe3ef] bg-white px-6 py-3.5">
        <div className="flex items-center gap-3">
          <div className="grid h-6 w-6 place-items-center rounded-md text-[#4481f3]">
            <AtomGlyph />
          </div>
          <div className="text-[14px] font-semibold">MRGPRX2 Classifier</div>
        </div>
        <div className="text-[14px] text-[#7d8393]">
          Multi-factor agonist / antagonist prediction
        </div>
      </header>

      <div className="grid h-[calc(100vh-57px)] grid-cols-[430px_minmax(0,1fr)_410px] overflow-hidden">
        <aside className="min-h-0 min-w-0 overflow-y-auto overflow-x-hidden border-r border-[#dfe3ef] bg-[#fafbff] px-5 py-6">
          <PanelTitle>Structure Input</PanelTitle>

          <textarea
            className="mt-4 min-h-[136px] w-full rounded-[20px] border border-[#dfe3ef] bg-white px-4 py-4 font-mono text-[14px] leading-6 text-[#22262f] shadow-[0_1px_2px_rgba(17,24,39,0.02)] outline-none transition placeholder:font-sans placeholder:text-[15px] placeholder:text-[#8a90a0] focus:border-[#7ea5f8]"
            value={query}
            onChange={(event) => {
              setQuery(event.target.value);
              setCompound(null);
              setPrediction(null);
              setStructureImageUrl(null);
              setError(null);
              setDockingResult(null);
              setDockingError(null);
              updateExperimentalData("dockingScore", "");
            }}
            placeholder="Paste a SMILES string here. e.g. C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O"
          />

          {compound ? (
            <div className="mt-3 flex w-full items-center justify-between rounded-[18px] border border-[#dfe3ef] bg-white px-4 py-3.5 shadow-[0_1px_2px_rgba(17,24,39,0.02)]">
              <span className="flex items-center gap-3 text-[14px] font-medium text-[#262b36]">
                <span className="text-[#7b8090]">
                  <FlaskGlyph />
                </span>
                {displayedCompoundName}
              </span>
            </div>
          ) : null}

          {query.trim() || compound ? (
            <div className="mt-3 rounded-[18px] border border-[#dfe3ef] bg-white px-4 py-3.5 shadow-[0_1px_2px_rgba(17,24,39,0.02)]">
              <div className="flex items-center gap-3 text-[14px] font-medium text-[#232833]">
                <span className="text-[#747b8d]">
                  <SparkGlyph />
                </span>
                Teach the Classifier
              </div>
              <p className="mt-2 text-[13px] leading-5 text-[#7b8090]">
                Know this compound&apos;s MRGPRX2 activity? Label it and its
                structure is folded into future similarity-based predictions.
              </p>

              {isCurrentStructureLabeled && labelResult && !editingLabel ? (
                <div className="mt-3 rounded-[14px] border border-[#9ebeff] bg-[#f4f8ff] px-3.5 py-3">
                  <p className="text-[13px] leading-5 text-[#4a84f6]">
                    Already saved as <strong>{labelResult.label}</strong>. The
                    structure is in the reference set (
                    {labelResult.total_labeled_compounds} user-labeled compound
                    {labelResult.total_labeled_compounds === 1 ? "" : "s"}{" "}
                    learned so far).
                  </p>
                  <button
                    type="button"
                    className="mt-2.5 text-[13px] font-medium text-[#4a84f6] underline-offset-2 hover:underline"
                    onClick={() => {
                      setEditingLabel(true);
                      setLabelError(null);
                    }}
                  >
                    Mislabeled? Change the label
                  </button>
                </div>
              ) : (
                <>
                  {editingLabel ? (
                    <p className="mt-2 text-[12.5px] leading-5 text-[#9aa0b0]">
                      Choosing a new label below will overwrite the saved one
                      for this structure.
                    </p>
                  ) : null}
                  <div className="mt-3 grid grid-cols-3 gap-2.5">
                    <button
                      type="button"
                      className="flex h-11 items-center justify-center rounded-[14px] border border-[#9bbcff] bg-[#f4f8ff] px-2 text-[13px] font-medium text-[#4a84f6] transition hover:bg-[#e9f1ff] disabled:cursor-not-allowed disabled:opacity-60"
                      onClick={() => onLabel("agonist")}
                      disabled={labelSubmitting !== null}
                    >
                      {labelSubmitting === "agonist"
                        ? "Saving..."
                        : "Label as Agonist"}
                    </button>
                    <button
                      type="button"
                      className="flex h-11 items-center justify-center rounded-[14px] border border-[#ff9fb0] bg-white px-2 text-[13px] font-medium text-[#ff4158] transition hover:bg-[#fff3f5] disabled:cursor-not-allowed disabled:opacity-60"
                      onClick={() => onLabel("antagonist")}
                      disabled={labelSubmitting !== null}
                    >
                      {labelSubmitting === "antagonist"
                        ? "Saving..."
                        : "Label as Antagonist"}
                    </button>
                    <button
                      type="button"
                      className="flex h-11 items-center justify-center rounded-[14px] border border-[#d4d9e7] bg-white px-2 text-[13px] font-medium text-[#6f7584] transition hover:bg-[#f4f5f9] disabled:cursor-not-allowed disabled:opacity-60"
                      onClick={() => onLabel("nonbinder")}
                      disabled={labelSubmitting !== null}
                    >
                      {labelSubmitting === "nonbinder"
                        ? "Saving..."
                        : "Label as Nonbinder"}
                    </button>
                  </div>
                  {editingLabel ? (
                    <button
                      type="button"
                      className="mt-2.5 text-[12.5px] font-medium text-[#9aa0b0] hover:text-[#7b8090]"
                      onClick={() => {
                        setEditingLabel(false);
                        setLabelError(null);
                      }}
                    >
                      Cancel
                    </button>
                  ) : null}
                </>
              )}
              {labelError ? (
                <p className="mt-3 text-[13px] leading-5 text-[#ff4158]">
                  {labelError}
                </p>
              ) : null}
            </div>
          ) : null}

          <div className="mt-4">
            <button
              type="button"
              className="flex h-14 w-full cursor-pointer items-center justify-between rounded-[18px] border border-[#dfe3ef] bg-white px-4 shadow-[0_1px_2px_rgba(17,24,39,0.02)]"
              onClick={() => setExperimentalDataOpen((v) => !v)}
            >
              <span className="flex items-center gap-3 text-[14px] font-medium text-[#232833]">
                <span className="text-[#747b8d]">
                  <BeakerGlyph />
                </span>
                Experimental Data
              </span>
              <span
                className={`text-[#7b8191] transition-transform ${experimentalDataOpen ? "rotate-180" : ""}`}
              >
                <ChevronGlyph />
              </span>
            </button>
            {experimentalDataOpen && <div className="mt-3 grid gap-3">
              <SectionLabel icon={<FlaskGlyph />}>
                AutoDock Vina
              </SectionLabel>
              <FormCard>
                <div className="flex flex-col gap-2">
                  <button
                    type="button"
                    disabled={!activeSmiles || dockingLoading}
                    onClick={async () => {
                      if (!activeSmiles) return;
                      setDockingLoading(true);
                      setDockingError(null);
                      setDockingResult(null);
                      try {
                        const result = await dockCompound(activeSmiles);
                        setDockingResult(result);
                        updateExperimentalData(
                          "dockingScore",
                          String(result.affinity_kcal_mol)
                        );
                        // Auto-predict with fresh docking scores (state hasn't
                        // flushed yet, so pass all values directly as overrides)
                        if (query.trim()) {
                          onPredict(
                            undefined,
                            result.affinity_kcal_mol,
                            result.delta_delta_score,
                            result.inactive_affinity_kcal_mol,
                          );
                        }
                      } catch (e) {
                        setDockingError(
                          e instanceof Error ? e.message : "Docking failed"
                        );
                      } finally {
                        setDockingLoading(false);
                      }
                    }}
                    className="flex items-center justify-center gap-2 rounded-xl border border-[#4a84f6] bg-[#f0f5ff] px-4 py-2.5 text-[13px] font-medium text-[#2563eb] transition hover:bg-[#dbeafe] disabled:cursor-not-allowed disabled:opacity-50"
                  >
                    {dockingLoading ? (
                      <>
                        <span className="inline-block h-4 w-4 animate-spin rounded-full border-2 border-[#2563eb] border-t-transparent" />
                        Running AutoDock Vina…
                      </>
                    ) : (
                      "Run AutoDock Vina (MRGPRX2 pocket)"
                    )}
                  </button>
                  {dockingLoading && (
                    <p className="text-[12px] text-[#6b7280] italic">
                      AutoDock Vina is running against the MRGPRX2 orthosteric pocket (PDB 7VDH).
                      This typically takes <span className="font-medium">10–30 seconds</span> depending on molecular size.
                      The prediction panel will update automatically when done.
                    </p>
                  )}
                  {dockingResult && !dockingLoading && (
                    <div className="rounded-lg bg-[#f0fdf4] border border-[#86efac] px-3 py-2 text-[13px] text-[#166534] space-y-1">
                      <div>
                        <span className="font-semibold">Active state (7VDH):</span>{" "}
                        {dockingResult.affinity_kcal_mol.toFixed(2)} kcal/mol
                        {" "}({dockingResult.num_modes} pose{dockingResult.num_modes !== 1 ? "s" : ""})
                      </div>
                      {dockingResult.inactive_affinity_kcal_mol != null && (
                        <div>
                          <span className="font-semibold">Inactive state (AlphaFold2):</span>{" "}
                          {dockingResult.inactive_affinity_kcal_mol.toFixed(2)} kcal/mol
                        </div>
                      )}
                      {dockingResult.delta_delta_score != null && (
                        <div>
                          <span className="font-semibold">ΔΔScore:</span>{" "}
                          {dockingResult.delta_delta_score > 0 ? "+" : ""}
                          {dockingResult.delta_delta_score.toFixed(2)} kcal/mol
                          {" — "}
                          {dockingResult.delta_delta_score < -1.0
                            ? "active-state preference → agonist signal"
                            : dockingResult.delta_delta_score > 1.0
                            ? "inactive-state preference → antagonist signal"
                            : "state-indifferent (no directional signal)"}
                        </div>
                      )}
                      <div className="text-[#4ade80] font-medium">Prediction updated automatically.</div>
                    </div>
                  )}
                  {dockingError && (
                    <div className={`rounded-lg px-3 py-2 text-[13px] space-y-1 ${
                      dockingError.includes("manually")
                        ? "bg-[#fffbeb] border border-[#fcd34d] text-[#92400e]"
                        : "bg-[#fef2f2] border border-[#fca5a5] text-[#991b1b]"
                    }`}>
                      {dockingError.includes("manually") ? (
                        <>
                          <div className="font-medium">AutoDock Vina could not prepare this molecule for 3D docking.</div>
                          <div className="text-[12px]">
                            This SMILES uses aromatic lactone/chromone notation that RDKit cannot embed into 3D coordinates.
                            Please enter the docking score manually in the field below (from an external docking tool), or leave it blank to use the structure-only prediction.
                          </div>
                        </>
                      ) : (
                        dockingError
                      )}
                    </div>
                  )}
                </div>
                <Field
                  label={dockingError?.includes("manually") ? "Docking Score (kcal/mol) — enter manually" : "Docking Score (kcal/mol)"}
                  value={experimentalData.dockingScore}
                  onChange={(value) =>
                    updateExperimentalData("dockingScore", value)
                  }
                  placeholder="e.g. -6.5 (or run Vina above)"
                />
              </FormCard>

              <div className="flex items-center gap-2">
                <SectionLabel icon={<SparkGlyph />}>mRNA Expression</SectionLabel>
                <span className="rounded px-1.5 py-0.5 text-[10px] font-semibold uppercase tracking-wide bg-amber-100 text-amber-700 border border-amber-300">
                  Reference Only
                </span>
              </div>
              <FormCard className="border-amber-200 bg-amber-50/40">
                <p className="text-[11px] text-amber-700 mb-2">
                  Not used in the agonist/antagonist classifier. Displayed separately as receptor regulation data.
                </p>
                <Field
                  label="Fold Change (vs control)"
                  value={experimentalData.mrnaFoldChange}
                  onChange={(value) =>
                    updateExperimentalData("mrnaFoldChange", value)
                  }
                  placeholder="e.g. 2.5"
                />
                <SelectField
                  label="Method"
                  value={experimentalData.mrnaMethod}
                  onChange={(value) =>
                    updateExperimentalData("mrnaMethod", value)
                  }
                  options={["Select method", "qPCR", "RNA-seq", "Microarray"]}
                />
              </FormCard>

              <div className="flex items-center gap-2">
                <SectionLabel icon={<TubeGlyph />}>Protein Expression</SectionLabel>
                <span className="rounded px-1.5 py-0.5 text-[10px] font-semibold uppercase tracking-wide bg-amber-100 text-amber-700 border border-amber-300">
                  Reference Only
                </span>
              </div>
              <FormCard className="border-amber-200 bg-amber-50/40">
                <p className="text-[11px] text-amber-700 mb-2">
                  Not used in the agonist/antagonist classifier. Displayed separately as receptor regulation data.
                </p>
                <Field
                  label="Fold Change (vs control)"
                  value={experimentalData.proteinFoldChange}
                  onChange={(value) =>
                    updateExperimentalData("proteinFoldChange", value)
                  }
                  placeholder="e.g. 1.8"
                />
                <SelectField
                  label="Method"
                  value={experimentalData.proteinMethod}
                  onChange={(value) =>
                    updateExperimentalData("proteinMethod", value)
                  }
                  options={[
                    "Select method",
                    "Western blot",
                    "ELISA",
                    "Flow Cytometry",
                    "Mass Spec",
                  ]}
                />
              </FormCard>

              <SectionLabel>Conditions</SectionLabel>
              <FormCard>
                <Field
                  label="Cell Line"
                  value={experimentalData.cellLine}
                  onChange={(value) =>
                    updateExperimentalData("cellLine", value)
                  }
                  placeholder="e.g. HMC-1, LAD2, RBL-2H3"
                />
                <div className="grid grid-cols-[minmax(0,1fr)_minmax(0,1fr)] gap-3">
                  <Field
                    label="Conc. (μM)"
                    value={experimentalData.concentration}
                    onChange={(value) =>
                      updateExperimentalData("concentration", value)
                    }
                    placeholder="e.g. 10"
                  />
                  <Field
                    label="Time (hrs)"
                    value={experimentalData.timeHours}
                    onChange={(value) =>
                      updateExperimentalData("timeHours", value)
                    }
                    placeholder="e.g. 24"
                  />
                </div>
              </FormCard>
            </div>}
          </div>

          <div className="mt-4">
            <button
              type="button"
              className="flex h-14 w-full cursor-pointer items-center justify-between rounded-[18px] border border-[#dfe3ef] bg-white px-4 shadow-[0_1px_2px_rgba(17,24,39,0.02)]"
              onClick={() => setBindingResiduesOpen((v) => !v)}
            >
              <span className="flex items-center gap-3 text-[14px] font-medium text-[#232833]">
                <span className="text-[#747b8d]">
                  <TargetGlyph />
                </span>
                Binding Site Residues
              </span>
              <span
                className={`text-[#7b8191] transition-transform ${bindingResiduesOpen ? "rotate-180" : ""}`}
              >
                <ChevronGlyph />
              </span>
            </button>
            {bindingResiduesOpen && <div className="mt-3 rounded-[22px] border border-[#dfe3ef] bg-white px-4 py-4 shadow-[0_1px_2px_rgba(17,24,39,0.02)]">
              <div className="text-[13px] text-[#6f7584]">Binding Mode</div>
              <div className="mt-2">
                <select
                  className="h-11 w-full rounded-xl border border-[#d7dceb] bg-white px-4 text-[14px] text-[#262b36] outline-none"
                  value={bindingMode}
                  onChange={(event) =>
                    setBindingMode(
                      event.target.value as (typeof BINDING_MODES)[number],
                    )
                  }
                >
                  {BINDING_MODES.map((option) => (
                    <option key={option} value={option}>
                      {option}
                    </option>
                  ))}
                </select>
              </div>

              <div className="mt-4 text-[13px] font-medium text-[#6f7584]">
                Key Residue Interactions
              </div>

              <div className="mt-3 grid gap-3">
                {RESIDUE_OPTIONS.map((residue) => {
                  const isSelected = selectedResidues.includes(residue.id);

                  return (
                    <div
                      key={residue.id}
                      className="flex items-center justify-between gap-3"
                    >
                      <div className="min-w-0">
                        <div className="text-[14px] font-medium text-[#171a22]">
                          {residue.name}
                        </div>
                        <div className="mt-1 flex flex-wrap gap-2">
                          <Pill
                            tone={
                              residue.strength === "critical"
                                ? "red"
                                : residue.strength === "important"
                                  ? "blue"
                                  : "neutral"
                            }
                          >
                            {residue.strength}
                          </Pill>
                          <Pill tone="neutral">{residue.site}</Pill>
                        </div>
                      </div>
                      <Toggle
                        checked={isSelected}
                        onChange={() => toggleResidue(residue.id)}
                      />
                    </div>
                  );
                })}
              </div>
            </div>}
          </div>

          {error ? (
            <div className="mt-4 rounded-2xl border border-red-200 bg-red-50 px-4 py-3 text-sm text-red-700">
              {error}
            </div>
          ) : null}

          <button
            className="mt-6 flex h-14 w-full items-center justify-center gap-3 rounded-2xl bg-[#4a84f6] px-4 text-[15px] font-medium text-white shadow-[0_8px_20px_rgba(74,132,246,0.18)] transition hover:bg-[#3f79ec] disabled:cursor-not-allowed disabled:bg-[#aac2f5]"
            onClick={() => onPredict()}
            disabled={!canPredict || loadingPredict}
          >
            <SparkGlyph />
            {loadingPredict
              ? "Predicting..."
              : "Predict (Structure + Pharmacophore)"}
          </button>

          <div className="mt-4">
            <button
              className="flex h-16 w-full items-center justify-center gap-2.5 rounded-2xl border border-dashed border-[#d4d9e7] bg-white px-3 text-[14px] text-[#707789] transition hover:border-[#8faef9] hover:text-[#4a84f6] disabled:cursor-not-allowed disabled:opacity-60"
              onClick={onUploadClick}
              disabled={fileProcessing}
            >
              <UploadGlyph />
              {fileProcessing ? "Reading file..." : "Upload SDF / MOL2"}
            </button>
          </div>
          <input
            ref={fileInputRef}
            type="file"
            accept=".sdf,.mol2"
            className="hidden"
            onChange={onFileChange}
          />
          {fileProcessingMessage ? (
            <p className="mt-2.5 rounded-[14px] border border-[#e3e7f0] bg-[#fbfcff] px-3.5 py-2.5 text-[13px] leading-5 text-[#6f7584]">
              {fileProcessingMessage}
            </p>
          ) : null}

          <div className="mt-5">
            <Subheading>Quick Examples</Subheading>
            <div className="mt-3 grid gap-3">
              {QUICK_EXAMPLES.map((example) => (
                <button
                  key={example.name}
                  className="block w-full min-w-0 text-left"
                  onClick={() => onSelectExample(example)}
                >
                  <div className="text-[14px] font-medium text-[#151922]">
                    {example.name}
                  </div>
                  <div className="min-w-0 overflow-hidden text-ellipsis whitespace-nowrap font-mono text-[12px] text-[#7a8090]">
                    {example.smiles}
                  </div>
                </button>
              ))}
            </div>
          </div>

          {analysisHistory.length > 0 ? (
            <div className="mt-8 border-t border-[#dfe3ef] pt-8">
              <Subheading icon={<ClockGlyph />}>Analysis History</Subheading>
              <div className="mt-4 grid gap-3">
                {analysisHistory.map((entry) => (
                  <button
                    key={`${entry.label}-${entry.smiles}`}
                    type="button"
                    onClick={() => onPredict(entry.query)}
                    disabled={loadingPredict}
                    className="w-full rounded-[20px] border border-[#dfe3ef] bg-white px-5 py-4 text-left shadow-[0_1px_2px_rgba(17,24,39,0.02)] transition hover:border-[#9bbcff] hover:bg-[#f8faff] disabled:cursor-not-allowed disabled:opacity-60"
                  >
                    <div
                      className="text-[16px] font-semibold"
                      style={{ color: historyLabelColor(entry.label) }}
                    >
                      {entry.label}
                    </div>
                    <div className="mt-2 break-all font-mono text-[12px] leading-5 text-[#6f7584]">
                      {entry.smiles}
                    </div>
                  </button>
                ))}
              </div>
            </div>
          ) : null}
        </aside>

        <section className="flex min-w-0 min-h-0 flex-col overflow-hidden bg-white">
          {prediction ? (
            <div className="border-b border-[#dfe3ef] px-8 py-7">
              <div className="flex items-center gap-4">
                <span
                  className="text-[24px]"
                  style={{ color: predictionColor(displayedLabel) }}
                >
                  •
                </span>
                <span
                  className="text-[22px] font-semibold"
                  style={{ color: predictionColor(displayedLabel) }}
                >
                  {formatPredictionLabel(displayedLabel)}
                </span>
                <span className="text-[22px] font-semibold text-[#151922]">
                  {displayedConfidence.toFixed(1)}%
                </span>
                <span className="text-[16px] text-[#767d8d]">confidence</span>
              </div>
              {displayedLabel === "indeterminate" ? (
                <p className="mt-3 max-w-[640px] text-[14px] leading-6 text-[#a9711f]">
                  {prediction?.interpretation}
                </p>
              ) : null}
            </div>
          ) : null}

          <div className="flex flex-1 flex-col items-center justify-center overflow-hidden px-8 py-10">
            <div className="flex h-full w-full max-w-[720px] flex-col items-center justify-center">
              <div className="flex min-h-[620px] w-full items-center justify-center rounded-[28px] border border-[#dfe3ef] bg-[#fdfdff] px-8 py-10">
                {loadingPredict ? (
                  <AnalyzingState />
                ) : structureImageUrl ? (
                  <StructureGraphic markup={structureImageUrl} />
                ) : (
                  <EmptyStructureState />
                )}
              </div>
              <div className="mt-4 font-mono text-[13px] text-[#6f7584]">
                {loadingPredict ? "" : (compound?.smiles ?? "")}
              </div>
            </div>
          </div>
        </section>

        <aside className="min-h-0 min-w-0 overflow-y-auto overflow-x-hidden border-l border-[#dfe3ef] bg-[#fafbff] px-5 py-6">
          {loadingPredict ? (
            <div className="flex h-full min-h-[70vh] items-center justify-center px-6 text-center">
              <div className="max-w-[260px] text-[#7b8090]">
                <div className="mb-4 flex justify-center">
                  <Spinner />
                </div>
                <p className="text-[15px] leading-7">
                  Running the multi-factor analysis - results and structure will
                  appear together once everything is ready.
                </p>
              </div>
            </div>
          ) : prediction ? (
            <>
              {prediction.applicability_domain && !prediction.applicability_domain.in_domain && (
                <div className="mb-4 rounded-2xl border border-amber-200 bg-amber-50 px-4 py-3">
                  <div className="text-[12px] font-semibold uppercase tracking-[0.1em] text-amber-700">
                    Outside Applicability Domain
                  </div>
                  <p className="mt-1.5 text-[13px] leading-5 text-amber-800">
                    {prediction.applicability_domain.reason}
                  </p>
                </div>
              )}

              {prediction.assay_basis && (
                <div className="mb-4 rounded-2xl border border-[#e3e7f1] bg-[#f7f8fc] px-4 py-3">
                  <div className="text-[12px] font-semibold uppercase tracking-[0.1em] text-[#71788a]">
                    Prediction Basis
                  </div>
                  <p className="mt-1.5 text-[12px] leading-5 text-[#6f7584]">
                    <span className="font-medium">Readout:</span> {prediction.assay_basis.readout}
                  </p>
                  <p className="mt-1 text-[12px] leading-5 text-[#6f7584]">
                    {prediction.assay_basis.note}
                  </p>
                </div>
              )}

              {prediction.receptor_regulation && (
                <div className="mb-4 rounded-2xl border border-purple-100 bg-purple-50 px-4 py-3">
                  <div className="text-[12px] font-semibold uppercase tracking-[0.1em] text-purple-700">
                    Receptor Regulation (independent axis)
                  </div>
                  {prediction.receptor_regulation.mrna_note && (
                    <p className="mt-1.5 text-[13px] text-purple-800">{prediction.receptor_regulation.mrna_note}</p>
                  )}
                  {prediction.receptor_regulation.protein_note && (
                    <p className="mt-1 text-[13px] text-purple-800">{prediction.receptor_regulation.protein_note}</p>
                  )}
                  <p className="mt-2 text-[12px] leading-5 text-purple-700 italic">
                    {prediction.receptor_regulation.warning}
                  </p>
                </div>
              )}

              <SideSection title="Factor Analysis">
                <PanelCard>
                  <div className="grid gap-4">
                    {factorAnalysis.map((factor) => (
                      <MetricBar
                        key={factor.label}
                        label={factor.label}
                        weight={`w=${factor.weight.toFixed(1)}`}
                        value={factor.score}
                      />
                    ))}
                  </div>
                  <div className="my-5 border-t border-[#eceff7]" />
                  <MetricBar
                    label="Composite Score"
                    weight=""
                    value={compositeScore}
                    emphasis
                  />

                  <div className="mt-5">
                    <div className="text-[13px] font-medium uppercase tracking-[0.1em] text-[#71788a]">
                      Weighted Contribution Breakdown
                    </div>
                    <div className="mt-3 flex h-4 w-full overflow-hidden rounded-full bg-[#edf0f7]">
                      {weightedFactors.map((factor, index) => (
                        <div
                          key={factor.label}
                          className="h-full"
                          style={{
                            width: `${factor.contribution_pct}%`,
                            backgroundColor:
                              FACTOR_COLORS[index % FACTOR_COLORS.length],
                          }}
                          title={`${factor.label}: ${factor.contribution_pct.toFixed(1)}% of composite`}
                        />
                      ))}
                    </div>
                    <div className="mt-3 grid gap-2">
                      {weightedFactors.map((factor, index) => (
                        <div
                          key={factor.label}
                          className="flex items-center justify-between gap-3 text-[13px]"
                        >
                          <span className="flex min-w-0 items-center gap-2 text-[#6f7584]">
                            <span
                              className="h-2.5 w-2.5 shrink-0 rounded-full"
                              style={{
                                backgroundColor:
                                  FACTOR_COLORS[index % FACTOR_COLORS.length],
                              }}
                            />
                            <span className="truncate">{factor.label}</span>
                            <span className="shrink-0 text-[#a0a6b6]">
                              w={factor.weight.toFixed(1)}
                            </span>
                          </span>
                          <span className="shrink-0 font-medium text-[#171a22]">
                            {factor.contribution_pct.toFixed(1)}%
                          </span>
                        </div>
                      ))}
                    </div>
                  </div>

                  <p className="mt-5 text-[15px] italic leading-6 text-[#7b8090]">
                    Add experimental data &amp; binding site interactions for
                    higher-confidence multi-factor prediction.
                  </p>
                </PanelCard>
              </SideSection>

              <SideSection title="Structural & Binding Analysis">
                <PanelCard>
                  <LikelihoodRow
                    label="Agonist Likelihood"
                    value={prediction.binding_analysis.agonist_likelihood}
                    color="#4a84f6"
                  />
                  <LikelihoodRow
                    label="Antagonist Likelihood"
                    value={prediction.binding_analysis.antagonist_likelihood}
                    color="#ff4158"
                  />

                  <div className="mt-8 grid gap-7">
                    {featureGroups.map((group) => (
                      <div key={group.title}>
                        <div className="text-[13px] font-medium uppercase tracking-[0.12em] text-[#71788a]">
                          {group.title}
                        </div>
                        <div className="mt-3 flex flex-wrap gap-2">
                          {group.items.map((item) => (
                            <FeatureTag key={item}>{item}</FeatureTag>
                          ))}
                        </div>
                      </div>
                    ))}
                  </div>

                  <div className="mt-8 rounded-[18px] bg-[#f7f8fc] p-5">
                    <div className="flex items-start justify-between gap-4">
                      <div className="text-[12px] uppercase tracking-[0.08em] text-[#707789]">
                        Most Similar Known Ligand
                      </div>
                      <div className="text-right text-[14px] text-[#7c8393]">
                        <div>
                          {prediction.similar_ligand.similarity.toFixed(1)}%
                        </div>
                        <div>match</div>
                      </div>
                    </div>
                    <div className="mt-3 text-[16px] font-semibold text-[#171a22]">
                      {prediction.similar_ligand.name}{" "}
                      <span className="ml-2 rounded-full border border-[#9ebeff] bg-[#f4f8ff] px-2.5 py-0.5 text-[12px] font-medium text-[#4a84f6]">
                        {prediction.similar_ligand.label}
                      </span>
                    </div>
                    <p className="mt-3 text-[14px] leading-6 text-[#6f7584]">
                      {prediction.similar_ligand.rationale}
                    </p>
                  </div>

                  <div className="mt-6 border-t border-[#eceff7] pt-6">
                    <ProbabilityChart probabilities={displayedProbabilities} />
                  </div>
                </PanelCard>
              </SideSection>

              {mrgprx1Comparison ? (
                <SideSection title="MRGPRX2 vs MRGPRX1 Selectivity">
                  <PanelCard>
                    <div className="mb-1 flex items-center justify-between gap-3">
                      <span className="text-[15px] text-[#6f7584]">
                        Selectivity estimate
                      </span>
                      <SelectivityBadge
                        score={mrgprx1Comparison.mrgprx2_selectivity_score}
                      />
                    </div>
                    <p className="mb-4 text-[13px] leading-6 text-[#9aa0b0]">
                      A heuristic read of structural divergence vs. the MRGPRX1
                      paralog pocket - not a measured binding selectivity. See
                      the signals below for the structural evidence behind this
                      call.
                    </p>

                    <div className="mt-2 grid gap-3 sm:grid-cols-2">
                      <PocketCard
                        title="MRGPRX2 pocket"
                        feature={mrgprx1Comparison.mrgprx2_pocket.feature}
                        residues={mrgprx1Comparison.mrgprx2_pocket.residues}
                        signature={
                          mrgprx1Comparison.mrgprx2_pocket.structural_signature
                        }
                        tone="blue"
                      />
                      <PocketCard
                        title="MRGPRX1 pocket"
                        feature={mrgprx1Comparison.mrgprx1_pocket.feature}
                        residues={mrgprx1Comparison.mrgprx1_pocket.residues}
                        signature={
                          mrgprx1Comparison.mrgprx1_pocket.structural_signature
                        }
                        tone="neutral"
                      />
                    </div>

                    <div className="mt-4 rounded-[16px] bg-[#f7f8fc] p-4 text-[13px] leading-6 text-[#6f7584]">
                      <span className="font-medium text-[#171a22]">
                        Conserved orthosteric anchor:{" "}
                      </span>
                      MRGPRX2{" "}
                      {mrgprx1Comparison.shared_residues.mrgprx2.join(", ")} ≈
                      MRGPRX1{" "}
                      {mrgprx1Comparison.shared_residues.mrgprx1.join(", ")}.{" "}
                      {mrgprx1Comparison.shared_residues.note}
                    </div>

                    {mrgprx1Comparison.closest_mrgprx1_ligand ? (
                      <div className="mt-4 flex items-center justify-between gap-3 text-[14px]">
                        <span className="text-[#71788a]">
                          Closest MRGPRX1 reference
                        </span>
                        <span className="font-medium text-[#171a22]">
                          {mrgprx1Comparison.closest_mrgprx1_ligand} (
                          {mrgprx1Comparison.mrgprx1_similarity.toFixed(1)}%)
                        </span>
                      </div>
                    ) : null}

                    <div className="mt-4 grid gap-2">
                      {mrgprx1Comparison.signals.map((signal) => (
                        <div
                          key={signal}
                          className="flex items-start gap-2 text-[13px] leading-6 text-[#6f7584]"
                        >
                          <span className="mt-2 h-1.5 w-1.5 shrink-0 rounded-full bg-[#9bbcff]" />
                          {signal}
                        </div>
                      ))}
                    </div>
                  </PanelCard>
                </SideSection>
              ) : null}

              {modelEvidence ? (
                <SideSection title="Model Evidence">
                  <PanelCard>
                    <div className="flex items-start justify-between gap-4">
                      <div>
                        <div className="text-[16px] font-semibold text-[#171a22]">
                          {modelEvidence.model_type}
                        </div>
                        <div className="mt-1 text-[13px] text-[#9aa0b0]">
                          Trained on {modelEvidence.trained_on} labeled MRGPRX2
                          ligands
                        </div>
                      </div>
                      <span className="rounded-full bg-[#edf4ff] px-3 py-1 text-[12px] font-medium text-[#4a84f6]">
                        Learned
                      </span>
                    </div>

                    <p className="mt-4 text-[14px] leading-6 text-[#6f7584]">
                      {modelEvidence.summary}
                    </p>

                    <LikelihoodRow
                      label="Nonbinder Probability"
                      value={modelEvidence.nonbinder_probability}
                      color="#8a90a0"
                    />
                    <LikelihoodRow
                      label="Agonist Probability"
                      value={modelEvidence.agonist_probability}
                      color="#4a84f6"
                    />
                    <LikelihoodRow
                      label="Antagonist Probability"
                      value={modelEvidence.antagonist_probability}
                      color="#ff4158"
                    />

                    {modelEvidence.experimental_adjustment ? (
                      <div className="mt-6 rounded-[14px] border border-[#e3e7f1] bg-[#f7f8fc] px-4 py-3.5">
                        <div className="text-[13px] font-medium uppercase tracking-[0.12em] text-[#71788a]">
                          Calibrated by your experimental data
                        </div>
                        <div className="mt-2 flex items-center gap-2 text-[14px] text-[#171a22]">
                          <span className="font-medium capitalize">
                            {modelEvidence.experimental_adjustment.target_label}
                          </span>
                          <span className="text-[#9aa0b0]">
                            {modelEvidence.experimental_adjustment.structure_probability.toFixed(1)}%
                            (structure-only)
                          </span>
                          <span className="text-[#9aa0b0]">→</span>
                          <span className="font-semibold">
                            {modelEvidence.experimental_adjustment.adjusted_probability.toFixed(1)}%
                            (adjusted)
                          </span>
                        </div>
                        <div className="mt-2 grid gap-1">
                          {modelEvidence.experimental_adjustment.components.map((component) => (
                            <div
                              key={component.description}
                              className="flex items-center justify-between text-[13px] text-[#6f7584]"
                            >
                              <span>{component.description}</span>
                              <span className="font-medium text-[#171a22]">
                                +{component.delta_pct.toFixed(1)}pp
                              </span>
                            </div>
                          ))}
                        </div>
                        <p className="mt-2 text-[12px] leading-5 text-[#9aa0b0]">
                          Experimental evidence shows MRGPRX2 engagement strength, not
                          direction - it calibrates the structure-grounded probability
                          within a bounded range (
                          {modelEvidence.experimental_adjustment.applied_delta_pct.toFixed(1)}
                          pp applied here). The agonist-vs-antagonist call itself stays
                          structure-based.
                        </p>
                        {modelEvidence.experimental_adjustment.specificity_note ? (
                          <p className="mt-2 rounded-[10px] bg-[#fff6e9] px-3 py-2 text-[12px] leading-5 text-[#9a6b1f]">
                            ⚠ {modelEvidence.experimental_adjustment.specificity_note}
                          </p>
                        ) : null}
                      </div>
                    ) : null}

                    {modelEvidence.nearest_neighbors.length > 0 ? (
                      <div className="mt-6">
                        <div className="text-[13px] font-medium uppercase tracking-[0.12em] text-[#71788a]">
                          Nearest labeled neighbors
                        </div>
                        <div className="mt-3 grid gap-2">
                          {modelEvidence.nearest_neighbors.map(
                            (neighbor, index) => (
                              <div
                                key={`${neighbor.name}-${index}`}
                                className="flex items-center justify-between gap-3 rounded-[14px] bg-[#f7f8fc] px-4 py-3 text-[14px]"
                              >
                                <div className="min-w-0">
                                  <div className="truncate font-medium text-[#171a22]">
                                    {neighbor.name}
                                  </div>
                                  <div className="text-[12px] text-[#9aa0b0]">
                                    {neighbor.source === "user"
                                      ? "from your labels"
                                      : "curated reference"}
                                  </div>
                                </div>
                                <div className="flex items-center gap-2 text-right">
                                  <span
                                    className="rounded-full border px-2.5 py-0.5 text-[12px] font-medium"
                                    style={{
                                      color: historyLabelColor(
                                        neighbor.label.charAt(0).toUpperCase() +
                                          neighbor.label.slice(1),
                                      ),
                                      borderColor: "currentColor",
                                    }}
                                  >
                                    {neighbor.label}
                                  </span>
                                  <span className="text-[#7c8393]">
                                    {neighbor.similarity.toFixed(1)}%
                                  </span>
                                </div>
                              </div>
                            ),
                          )}
                        </div>
                      </div>
                    ) : null}

                    {modelEvidence.top_features.length > 0 ? (
                      <div className="mt-6">
                        <div className="text-[13px] font-medium uppercase tracking-[0.12em] text-[#71788a]">
                          Features driving this call
                        </div>
                        <div className="mt-3 grid gap-3">
                          {modelEvidence.top_features.map((feature) => (
                            <div key={feature.description}>
                              <div className="flex items-center justify-between text-[13px] text-[#6f7584]">
                                <span>{feature.description}</span>
                                <span className="font-medium text-[#171a22]">
                                  {feature.importance.toFixed(1)}%
                                </span>
                              </div>
                              <div className="mt-1.5 h-2 overflow-hidden rounded-full bg-[#eceff7]">
                                <div
                                  className="h-full rounded-full bg-[#9bbcff]"
                                  style={{
                                    width: `${Math.min(100, feature.importance)}%`,
                                  }}
                                />
                              </div>
                            </div>
                          ))}
                        </div>
                      </div>
                    ) : null}
                  </PanelCard>
                </SideSection>
              ) : null}

              <SideSection title="Molecular Properties">
                <div className="mb-4 flex justify-end">
                  <span className="rounded-full bg-[#edf4ff] px-4 py-2 text-[14px] font-medium text-[#4a84f6]">
                    ✓ PubChem
                  </span>
                </div>
                <PanelCard>
                  <div className="grid gap-0">
                    {molecularProperties.map((item, index) => (
                      <div
                        key={item.label}
                        className={`grid grid-cols-[1fr_auto] items-center gap-4 py-4 ${
                          index === 0 ? "" : "border-t border-[#eceff7]"
                        }`}
                      >
                        <span className="text-[16px] text-[#71788a]">
                          {item.label}
                        </span>
                        <span className="text-[16px] font-medium text-[#151922]">
                          {item.value}
                        </span>
                      </div>
                    ))}
                  </div>
                </PanelCard>
              </SideSection>

              <SideSection title="Drug-likeness">
                <PanelCard>
                  <div className="grid gap-4 py-2">
                    {prediction.drug_likeness.map((item) => (
                      <div
                        key={item.label}
                        className="flex items-center justify-between gap-4 text-[16px]"
                      >
                        <span className="text-[#71788a]">{item.label}</span>
                        <span
                          className={`font-medium ${
                            item.passed ? "text-[#4a84f6]" : "text-[#a0a6b6]"
                          }`}
                        >
                          {item.passed ? "✓ Pass" : "Fail"}
                        </span>
                      </div>
                    ))}
                  </div>
                </PanelCard>
                <div className="mt-10 text-center text-[16px] text-[#767d8d]">
                  {formatAnalyzedAt(prediction.analyzed_at)}
                </div>
              </SideSection>
            </>
          ) : (
            <div className="flex h-full min-h-[70vh] items-center justify-center px-6 text-center">
              <div className="max-w-[240px] text-[#7b8090]">
                <div className="mb-4 flex justify-center text-[#707789]">
                  <InfoGlyph />
                </div>
                <p className="text-[15px] leading-7">
                  Properties, pharmacophore analysis, and binding site
                  evaluation will appear here after prediction.
                </p>
              </div>
            </div>
          )}
        </aside>
      </div>
    </main>
  );
}

function predictionColor(label: string) {
  if (label === "agonist-like") return "#4a84f6";
  if (label === "antagonist-like") return "#ff4158";
  if (label === "indeterminate") return "#c98a1f";
  return "#8a90a0";
}

function historyLabelColor(label: string) {
  if (label.startsWith("Antagonist")) return "#ff4158";
  if (label.startsWith("Nonbinder")) return "#8a90a0";
  if (label.startsWith("Indeterminate")) return "#c98a1f";
  if (label.startsWith("Agonist")) return "#4a84f6";
  return "#4a84f6";
}

function formatPredictionLabel(label: string) {
  if (label === "agonist-like") return "Agonist";
  if (label === "antagonist-like") return "Antagonist";
  if (label === "nonbinder-like") return "Nonbinder";
  if (label === "indeterminate") return "Indeterminate";
  return label;
}

function parseNullableNumber(value: string) {
  const normalized = value.trim();
  if (!normalized) {
    return null;
  }

  const parsed = Number(normalized);
  return Number.isFinite(parsed) ? parsed : null;
}

function formatAnalyzedAt(value: string) {
  const parsed = new Date(value);
  if (Number.isNaN(parsed.getTime())) {
    return "Analyzed";
  }

  return `Analyzed ${parsed.toLocaleTimeString([], {
    hour: "numeric",
    minute: "2-digit",
  })}`;
}

function PanelTitle({ children }: { children: ReactNode }) {
  return (
    <div className="text-[17px] font-medium uppercase tracking-[0.06em] text-[#6f7584]">
      {children}
    </div>
  );
}

function Subheading({
  children,
  icon,
}: {
  children: ReactNode;
  icon?: ReactNode;
}) {
  return (
    <div className="flex items-center gap-3 text-[17px] font-medium text-[#6f7584]">
      {icon ? <span>{icon}</span> : null}
      {children}
    </div>
  );
}

function SideSection({
  title,
  children,
}: {
  title: string;
  children: ReactNode;
}) {
  return (
    <section className="mb-10">
      <div className="mb-4 text-[17px] font-medium uppercase tracking-[0.06em] text-[#6f7584]">
        {title}
      </div>
      {children}
    </section>
  );
}

function PanelCard({ children }: { children: ReactNode }) {
  return (
    <div className="rounded-[24px] border border-[#dfe3ef] bg-white p-5 shadow-[0_1px_2px_rgba(17,24,39,0.03)]">
      {children}
    </div>
  );
}

function FormCard({ children, className }: { children: ReactNode; className?: string }) {
  return (
    <div className={`rounded-[22px] border border-[#e5e8f1] bg-white px-4 py-4 shadow-[0_1px_2px_rgba(17,24,39,0.02)] ${className ?? ""}`}>
      <div className="grid gap-3">{children}</div>
    </div>
  );
}

function SectionLabel({
  children,
  icon,
}: {
  children: ReactNode;
  icon?: ReactNode;
}) {
  return (
    <div className="flex items-center gap-3 text-[16px] text-[#6f7584]">
      {icon ? <span>{icon}</span> : null}
      {children}
    </div>
  );
}

function Field({
  label,
  value,
  onChange,
  placeholder,
}: {
  label: string;
  value: string;
  onChange: (value: string) => void;
  placeholder?: string;
}) {
  return (
    <label className="grid min-w-0 gap-1.5">
      <span className="text-[15px] text-[#6f7584]">{label}</span>
      <input
        className="h-12 w-full min-w-0 rounded-[16px] border border-[#e3e7f0] bg-[#fbfcff] px-4 text-[15px] text-[#202530] outline-none placeholder:text-[#aab0bf]"
        value={value}
        placeholder={placeholder}
        onChange={(event) => onChange(event.target.value)}
      />
    </label>
  );
}

function SelectField({
  label,
  value,
  onChange,
  options,
}: {
  label: string;
  value: string;
  onChange: (value: string) => void;
  options: readonly string[];
}) {
  return (
    <label className="grid min-w-0 gap-1.5">
      <span className="text-[15px] text-[#6f7584]">{label}</span>
      <div className="relative">
        <select
          className="h-12 w-full min-w-0 appearance-none rounded-[16px] border border-[#e3e7f0] bg-white px-4 text-[15px] text-[#202530] outline-none"
          value={value}
          onChange={(event) => onChange(event.target.value)}
        >
          {options.map((option) => (
            <option
              key={option}
              value={option === "Select method" ? "" : option}
            >
              {option}
            </option>
          ))}
        </select>
        <span className="pointer-events-none absolute right-5 top-1/2 -translate-y-1/2 text-[#8a90a0]">
          <ChevronGlyph />
        </span>
      </div>
    </label>
  );
}

function Pill({
  children,
  tone,
}: {
  children: ReactNode;
  tone: "blue" | "neutral" | "red";
}) {
  const classes =
    tone === "red"
      ? "border-[#ff9fb0] bg-white text-[#ff4158]"
      : tone === "blue"
        ? "border-[#9bbcff] bg-[#f4f8ff] text-[#4a84f6]"
        : "border-[#d7dceb] bg-white text-[#6f7584]";

  return (
    <span
      className={`rounded-full border px-2 py-0.5 text-[12px] leading-5 ${classes}`}
    >
      {children}
    </span>
  );
}

function selectivityTier(score: number): {
  label: string;
  tone: "blue" | "neutral" | "red";
} {
  if (score >= 65)
    return { label: "Predicted MRGPRX2-selective", tone: "blue" };
  if (score >= 40)
    return { label: "Moderate / mixed selectivity", tone: "neutral" };
  return { label: "Predicted MRGPRX1-leaning", tone: "red" };
}

function SelectivityBadge({ score }: { score: number }) {
  const tier = selectivityTier(score);
  return <Pill tone={tier.tone}>{tier.label}</Pill>;
}

function Toggle({
  checked,
  onChange,
}: {
  checked: boolean;
  onChange: () => void;
}) {
  return (
    <label className="relative inline-flex cursor-pointer items-center">
      <input
        type="checkbox"
        className="peer sr-only"
        checked={checked}
        onChange={() => onChange()}
      />
      <span className="h-8 w-14 rounded-full bg-[#dde0e8] transition peer-checked:bg-[#4a84f6]" />
      <span className="absolute left-1 h-6 w-6 rounded-full bg-white shadow-sm transition peer-checked:left-7" />
    </label>
  );
}

function MetricBar({
  label,
  weight,
  value,
  emphasis = false,
}: {
  label: string;
  weight: string;
  value: number;
  emphasis?: boolean;
}) {
  return (
    <div>
      <div className="mb-2 flex items-center justify-between gap-3">
        <div
          className={`text-[14px] ${emphasis ? "font-semibold text-[#171a22]" : "text-[#6f7584]"}`}
        >
          {label}
        </div>
        <div className="flex items-center gap-3">
          {weight ? (
            <span className="text-[14px] text-[#6f7584]">{weight}</span>
          ) : null}
          <span className="text-[14px] font-semibold text-[#171a22]">
            {value.toFixed(1)}
          </span>
        </div>
      </div>
      <div className="h-2 rounded-full bg-[#edf0f7]">
        <div
          className="h-2 rounded-full bg-[#4a84f6]"
          style={{ width: `${Math.min(100, value)}%` }}
        />
      </div>
    </div>
  );
}

function LikelihoodRow({
  label,
  value,
  color,
}: {
  label: string;
  value: number;
  color: string;
}) {
  return (
    <div className="mb-4">
      <div className="mb-2 flex items-center justify-between gap-3">
        <span className="text-[15px] text-[#6f7584]">{label}</span>
        <span className="text-[15px] font-medium" style={{ color }}>
          {value}%
        </span>
      </div>
      <div className="h-2 rounded-full bg-[#eef0f5]">
        <div
          className="h-2 rounded-full"
          style={{
            width: `${value}%`,
            backgroundColor: value === 0 ? "#eef0f5" : color,
          }}
        />
      </div>
    </div>
  );
}

function PocketCard({
  title,
  feature,
  residues,
  signature,
  tone,
}: {
  title: string;
  feature: string;
  residues: string[];
  signature: string;
  tone: "blue" | "neutral";
}) {
  const borderClass =
    tone === "blue"
      ? "border-[#9bbcff] bg-[#f4f8ff]"
      : "border-[#dfe3ef] bg-white";

  return (
    <div className={`rounded-[16px] border ${borderClass} p-4`}>
      <div className="text-[12px] font-medium uppercase tracking-[0.08em] text-[#707789]">
        {title}
      </div>
      <div className="mt-1 text-[14px] font-semibold text-[#171a22]">
        {feature}
      </div>
      <div className="mt-2 flex flex-wrap gap-1.5">
        {residues.map((residue) => (
          <Pill key={residue} tone="neutral">
            {residue}
          </Pill>
        ))}
      </div>
      <p className="mt-2 text-[12.5px] leading-5 text-[#6f7584]">{signature}</p>
    </div>
  );
}

function FeatureTag({ children }: { children: ReactNode }) {
  return (
    <span className="rounded-full border border-[#9ebeff] bg-[#f4f8ff] px-3 py-1.5 text-[13px] font-medium text-[#4a84f6]">
      ✓ {children}
    </span>
  );
}

function StructurePlaceholder() {
  return (
    <div className="flex h-full w-full flex-col items-center justify-center gap-2 px-8 text-center">
      <FlaskGlyph />
      <p className="text-[14px] text-[#9aa0b0]">
        Couldn&apos;t render a structure for this input.
      </p>
    </div>
  );
}

function StructureImage({ src }: { src: string }) {
  const [failed, setFailed] = useState(false);

  if (failed) {
    return <StructurePlaceholder />;
  }

  return (
    <div className="flex h-full w-full items-center justify-center">
      <img
        src={src}
        alt="Predicted structure"
        className="h-auto max-h-[360px] w-full max-w-[520px] object-contain"
        onError={() => setFailed(true)}
      />
    </div>
  );
}

function StructureGraphic({ markup }: { markup: string }) {
  const normalized = markup.trim();

  if (
    normalized.startsWith("data:image/") ||
    normalized.startsWith("http://") ||
    normalized.startsWith("https://") ||
    normalized.startsWith("/")
  ) {
    return <StructureImage src={normalized} />;
  }

  if (normalized.startsWith("<svg")) {
    const encodedSvg = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(
      normalized,
    )}`;

    return (
      <div className="flex h-full w-full items-center justify-center">
        <img
          src={encodedSvg}
          alt="Predicted structure"
          className="h-auto max-h-[360px] w-full max-w-[520px] object-contain"
        />
      </div>
    );
  }

  return (
    <div
      className="flex h-full w-full items-center justify-center"
      dangerouslySetInnerHTML={{ __html: normalized }}
    />
  );
}

function AnalyzingState() {
  return (
    <div className="flex flex-col items-center justify-center text-center">
      <Spinner size={28} />
      <p className="mt-6 text-[18px] font-semibold text-[#171a22]">
        Analyzing structure
      </p>
      <p className="mt-2 max-w-[360px] text-[15px] leading-7 text-[#7b8090]">
        Resolving the compound, rendering its 2D structure, and running the
        multi-factor prediction together.
      </p>
    </div>
  );
}

function Spinner({ size = 22 }: { size?: number }) {
  return (
    <svg
      width={size}
      height={size}
      viewBox="0 0 24 24"
      fill="none"
      className="animate-spin text-[#4a84f6]"
      aria-hidden="true"
    >
      <circle
        cx="12"
        cy="12"
        r="9"
        stroke="currentColor"
        strokeOpacity="0.18"
        strokeWidth="3"
      />
      <path
        d="M21 12a9 9 0 0 0-9-9"
        stroke="currentColor"
        strokeWidth="3"
        strokeLinecap="round"
      />
    </svg>
  );
}

function EmptyStructureState() {
  return (
    <div className="flex flex-col items-center justify-center text-center">
      <div className="grid h-28 w-28 place-items-center rounded-[24px] bg-[#f7f8fc] text-[#707789]">
        <FlaskGlyph />
      </div>
      <p className="mt-6 text-[18px] font-semibold text-[#171a22]">
        No structure loaded
      </p>
      <p className="mt-2 max-w-[360px] text-[15px] leading-7 text-[#7b8090]">
        Paste a SMILES string or upload a file to start analysis.
      </p>
    </div>
  );
}

function AtomGlyph() {
  return (
    <svg
      width="22"
      height="22"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <path
        d="M12 3c3.5 2.2 6 5.3 6 9s-2.5 6.8-6 9c-3.5-2.2-6-5.3-6-9s2.5-6.8 6-9Z"
        stroke="currentColor"
        strokeWidth="1.8"
      />
      <path
        d="M4 7c4 .1 7.4 1.2 9.5 4.3C15.6 14.4 16.3 18 16 22c-4-.1-7.4-1.2-9.5-4.3C4.4 14.6 3.7 11 4 7Z"
        stroke="currentColor"
        strokeWidth="1.8"
      />
      <path
        d="M20 7c-4 .1-7.4 1.2-9.5 4.3C8.4 14.4 7.7 18 8 22c4-.1 7.4-1.2 9.5-4.3 2.1-3.1 2.8-6.7 2.5-10.7Z"
        stroke="currentColor"
        strokeWidth="1.8"
      />
      <circle cx="12" cy="12" r="1.8" fill="currentColor" />
    </svg>
  );
}

function FlaskGlyph() {
  return (
    <svg
      width="18"
      height="18"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <path
        d="M9 3h6M10 3v6l-4 8a2 2 0 0 0 1.8 3h8.4A2 2 0 0 0 18 17l-4-8V3"
        stroke="currentColor"
        strokeWidth="1.8"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
}

function BeakerGlyph() {
  return (
    <svg
      width="18"
      height="18"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <path
        d="M8 3h8M9 3v5l-4 9a2 2 0 0 0 1.8 3h10.4A2 2 0 0 0 19 17l-4-9V3"
        stroke="currentColor"
        strokeWidth="1.8"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
      <path
        d="M8 13h8"
        stroke="currentColor"
        strokeWidth="1.8"
        strokeLinecap="round"
      />
    </svg>
  );
}

function TubeGlyph() {
  return (
    <svg
      width="18"
      height="18"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <path
        d="M6 3h4v10a3 3 0 1 1-4 0V3Zm8 0h4v10a3 3 0 1 1-4 0V3Z"
        stroke="currentColor"
        strokeWidth="1.8"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
}

function TargetGlyph() {
  return (
    <svg
      width="18"
      height="18"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <circle cx="12" cy="12" r="8" stroke="currentColor" strokeWidth="1.8" />
      <circle cx="12" cy="12" r="4" stroke="currentColor" strokeWidth="1.8" />
      <circle cx="12" cy="12" r="1.4" fill="currentColor" />
    </svg>
  );
}

function ClockGlyph() {
  return (
    <svg
      width="18"
      height="18"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <circle cx="12" cy="12" r="9" stroke="currentColor" strokeWidth="1.8" />
      <path
        d="M12 7v5l3 2"
        stroke="currentColor"
        strokeWidth="1.8"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
}

function ChevronGlyph() {
  return (
    <svg
      width="18"
      height="18"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <path
        d="m6 9 6 6 6-6"
        stroke="currentColor"
        strokeWidth="1.8"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
}

function SparkGlyph() {
  return (
    <svg
      width="18"
      height="18"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <path
        d="M12 3 13.8 8.2 19 10l-5.2 1.8L12 17l-1.8-5.2L5 10l5.2-1.8L12 3ZM5 18l.8 2.2L8 21l-2.2.8L5 24l-.8-2.2L2 21l2.2-.8L5 18Zm14-2 .8 2.2L22 19l-2.2.8L19 22l-.8-2.2L16 19l2.2-.8L19 16Z"
        fill="currentColor"
      />
    </svg>
  );
}

function UploadGlyph() {
  return (
    <svg
      width="18"
      height="18"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <path
        d="M12 16V4m0 0 4 4m-4-4-4 4M4 16v3a1 1 0 0 0 1 1h14a1 1 0 0 0 1-1v-3"
        stroke="currentColor"
        strokeWidth="1.8"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
}

function InfoGlyph() {
  return (
    <svg
      width="34"
      height="34"
      viewBox="0 0 24 24"
      fill="none"
      aria-hidden="true"
    >
      <circle
        cx="12"
        cy="12"
        r="9.25"
        stroke="currentColor"
        strokeWidth="1.5"
      />
      <circle cx="12" cy="8" r="1.1" fill="currentColor" />
      <path
        d="M12 11v6"
        stroke="currentColor"
        strokeWidth="1.5"
        strokeLinecap="round"
      />
    </svg>
  );
}
