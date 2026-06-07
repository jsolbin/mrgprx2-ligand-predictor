"use client";

type ProbabilityKey = "agonist" | "antagonist" | "nonbinder";

const LABELS: Record<ProbabilityKey, string> = {
  agonist: "Agonist",
  antagonist: "Antagonist",
  nonbinder: "Nonbinder",
};

const COLORS: Record<ProbabilityKey, string> = {
  agonist: "#4a84f6",
  antagonist: "#ff4158",
  nonbinder: "#a0a6b6",
};

export function ProbabilityChart({
  probabilities,
}: {
  probabilities: {
    agonist: number;
    antagonist: number;
    nonbinder: number;
  };
}) {
  const data = (Object.keys(probabilities) as ProbabilityKey[]).map((key) => ({
    key,
    label: LABELS[key],
    percent: Math.max(0, Math.min(1, probabilities[key])) * 100,
    color: COLORS[key],
  }));

  return (
    <div>
      <div className="mb-4 flex items-center justify-between">
        <span className="text-[13px] font-medium uppercase tracking-[0.12em] text-[#71788a]">
          Class Probabilities
        </span>
        <span className="text-[13px] text-[#a0a6b6]">Normalized 0-100%</span>
      </div>

      <div className="grid gap-4">
        {data.map((item) => (
          <div key={item.key}>
            <div className="mb-2 flex items-center justify-between gap-3">
              <span className="flex items-center gap-2.5 text-[15px] text-[#6f7584]">
                <span
                  className="h-2.5 w-2.5 shrink-0 rounded-full"
                  style={{ backgroundColor: item.color }}
                />
                {item.label}
              </span>
              <span className="text-[15px] font-medium" style={{ color: item.color }}>
                {item.percent.toFixed(1)}%
              </span>
            </div>
            <div className="h-2 rounded-full bg-[#eef0f5]">
              <div
                className="h-2 rounded-full transition-[width] duration-500"
                style={{
                  width: `${item.percent}%`,
                  backgroundColor: item.percent === 0 ? "#eef0f5" : item.color,
                }}
              />
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}
