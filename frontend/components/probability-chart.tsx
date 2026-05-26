"use client";

import {
  Bar,
  BarChart,
  CartesianGrid,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";

export function ProbabilityChart({
  probabilities,
}: {
  probabilities: {
    agonist: number;
    antagonist: number;
    nonbinder: number;
  };
}) {
  const data = [
    { name: "Agonist", probability: probabilities.agonist },
    { name: "Antagonist", probability: probabilities.antagonist },
    { name: "Nonbinder", probability: probabilities.nonbinder },
  ];

  return (
    <div className="h-64 rounded-3xl border border-black/10 bg-canvas/70 p-4">
      <ResponsiveContainer width="100%" height="100%">
        <BarChart data={data}>
          <CartesianGrid strokeDasharray="3 3" stroke="rgba(16,35,28,0.1)" />
          <XAxis dataKey="name" stroke="#3c544b" fontSize={12} />
          <YAxis
            stroke="#3c544b"
            fontSize={12}
            tickFormatter={(value: number) => `${Math.round(value * 100)}%`}
            domain={[0, 1]}
          />
          <Tooltip
            formatter={(value: number) => `${(value * 100).toFixed(1)}%`}
          />
          <Bar dataKey="probability" fill="#1d6b52" radius={[12, 12, 0, 0]} />
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
}
