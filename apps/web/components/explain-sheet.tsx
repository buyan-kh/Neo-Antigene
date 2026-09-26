"use client";

import { useEffect, useState } from "react";
import {
  Bar,
  BarChart,
  Cell,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip as ChartTooltip,
  XAxis,
  YAxis,
} from "recharts";

import { MutantPeptide } from "@/components/peptide";
import { Badge } from "@/components/ui/badge";
import {
  Sheet,
  SheetContent,
  SheetDescription,
  SheetHeader,
  SheetTitle,
} from "@/components/ui/sheet";
import { Separator } from "@/components/ui/separator";
import { Skeleton } from "@/components/ui/skeleton";
import { api, type Explanation } from "@/lib/api";

/** Two muted tones from the same family, not a red/green traffic light. */
const POSITIVE = "var(--chart-2)";
const NEGATIVE = "var(--chart-4)";

type Row = Record<string, unknown>;

export function ExplainSheet({
  runId,
  row,
  onClose,
}: {
  runId: string;
  row: Row | null;
  onClose: () => void;
}) {
  const [explanation, setExplanation] = useState<Explanation | null>(null);
  const [error, setError] = useState<string | null>(null);

  const peptide = row ? String(row.mutant_peptide ?? "") : "";
  const allele = row ? String(row.allele ?? "") : "";

  useEffect(() => {
    if (!peptide || !allele) return;
    setExplanation(null);
    setError(null);
    api
      .explain(runId, peptide, allele)
      .then(setExplanation)
      .catch((cause: Error) => setError(cause.message));
  }, [runId, peptide, allele]);

  const data = explanation
    ? [
        ...explanation.contributions.map((item) => ({
          feature: item.feature,
          contribution: item.contribution,
          weight: item.weight,
          value: item.value,
        })),
        {
          feature: "bias",
          contribution: explanation.bias,
          weight: explanation.bias,
          value: 1,
        },
      ].sort((a, b) => Math.abs(b.contribution) - Math.abs(a.contribution))
    : [];

  return (
    <Sheet open={row !== null} onOpenChange={(open) => !open && onClose()}>
      <SheetContent className="w-full overflow-y-auto sm:max-w-xl">
        <SheetHeader>
          <SheetTitle className="flex items-center gap-2">
            {row?.gene ? String(row.gene) : "Candidate"}
            <span className="font-mono text-sm font-normal text-muted-foreground">
              {row?.hgvsp ? String(row.hgvsp) : ""}
            </span>
          </SheetTitle>
          <SheetDescription>
            Every term that produced this score, with its sign.
          </SheetDescription>
        </SheetHeader>

        <div className="space-y-6 px-4 pb-8">
          <div className="flex items-start justify-between gap-4">
            <MutantPeptide
              mutant={peptide}
              wildtype={
                row?.wildtype_peptide ? String(row.wildtype_peptide) : null
              }
              position={Number(row?.mutation_position ?? 0)}
            />
            <Badge variant="secondary" className="font-mono">
              {allele}
            </Badge>
          </div>

          {error && <p className="text-sm text-destructive">{error}</p>}

          {!explanation && !error && <Skeleton className="h-72 w-full" />}

          {explanation && (
            <>
              <div className="grid grid-cols-3 gap-4 text-sm">
                <Stat label="Score" value={explanation.score.toFixed(4)} />
                <Stat label="Logit" value={explanation.logit.toFixed(3)} />
                <Stat label="Bias" value={explanation.bias.toFixed(2)} />
              </div>

              {explanation.gate_failures.length > 0 && (
                <div className="rounded-md border border-amber-500/40 bg-amber-500/5 p-3">
                  <p className="text-sm font-medium">Excluded by gates</p>
                  <div className="mt-2 flex flex-wrap gap-1.5">
                    {explanation.gate_failures.map((gate) => (
                      <Badge key={gate} variant="secondary" className="font-mono">
                        {gate}
                      </Badge>
                    ))}
                  </div>
                </div>
              )}

              <Separator />

              <div>
                <h3 className="text-sm font-medium">
                  Signed contribution to the logit
                </h3>
                <p className="mt-1 text-xs leading-relaxed text-muted-foreground">
                  Each bar is weight &times; normalized feature value. They sum,
                  with the bias, to the logit above. Nothing else enters the
                  score.
                </p>

                <div className="mt-4 h-[340px]">
                  <ResponsiveContainer width="100%" height="100%">
                    <BarChart
                      data={data}
                      layout="vertical"
                      margin={{ top: 4, right: 16, bottom: 4, left: 8 }}
                    >
                      <XAxis
                        type="number"
                        tick={{ fontSize: 11 }}
                        stroke="var(--border)"
                        tickLine={false}
                        axisLine={{ stroke: "var(--border)" }}
                      />
                      <YAxis
                        type="category"
                        dataKey="feature"
                        width={124}
                        tick={{ fontSize: 11, fontFamily: "var(--font-mono)" }}
                        stroke="var(--border)"
                        tickLine={false}
                        axisLine={false}
                      />
                      <ReferenceLine x={0} stroke="var(--border)" />
                      <ChartTooltip
                        cursor={{ fill: "var(--accent)", opacity: 0.4 }}
                        contentStyle={{
                          background: "var(--popover)",
                          border: "1px solid var(--border)",
                          borderRadius: 8,
                          fontSize: 12,
                        }}
                        formatter={(value, _name, entry) => {
                          const item = (entry as { payload: (typeof data)[number] })
                            .payload;
                          return [
                            `${Number(value).toFixed(3)}  (w=${item.weight.toFixed(2)} × ${item.value.toFixed(3)})`,
                            item.feature,
                          ];
                        }}
                      />
                      <Bar dataKey="contribution" radius={2} barSize={14}>
                        {data.map((item) => (
                          <Cell
                            key={item.feature}
                            fill={item.contribution >= 0 ? POSITIVE : NEGATIVE}
                          />
                        ))}
                      </Bar>
                    </BarChart>
                  </ResponsiveContainer>
                </div>
              </div>
            </>
          )}
        </div>
      </SheetContent>
    </Sheet>
  );
}

function Stat({ label, value }: { label: string; value: string }) {
  return (
    <div>
      <div className="text-xs text-muted-foreground">{label}</div>
      <div className="mt-0.5 font-mono tabular-nums">{value}</div>
    </div>
  );
}
