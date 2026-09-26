"use client";

import Link from "next/link";
import { useState } from "react";
import { Upload } from "lucide-react";
import { toast } from "sonner";
import {
  Bar,
  BarChart,
  Cell,
  ResponsiveContainer,
  Tooltip as ChartTooltip,
  XAxis,
  YAxis,
} from "recharts";

import { EmptyState } from "@/components/empty-state";
import { Skeleton } from "@/components/ui/skeleton";
import { useRun } from "@/components/run-context";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { Slider } from "@/components/ui/slider";
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table";
import { api, type BenchmarkReport } from "@/lib/api";

export default function BenchmarkScreen() {
  const { runId, restoring } = useRun();
  const [report, setReport] = useState<BenchmarkReport | null>(null);
  const [k, setK] = useState(20);
  const [busy, setBusy] = useState(false);

  async function load(file: File | null, nextK = k) {
    if (!file || !runId) return;
    setBusy(true);
    try {
      const body = new FormData();
      body.append("labels", file);
      setReport(await api.benchmark(runId, body));
      setK(nextK);
    } catch (error) {
      toast.error((error as Error).message);
    } finally {
      setBusy(false);
    }
  }

  if (restoring) {
    return (
      <div className="space-y-3">
        <Skeleton className="h-9 w-72" />
        <Skeleton className="h-[480px] w-full" />
      </div>
    );
  }

  if (!runId) {
    return (
      <EmptyState
        title="No run to benchmark"
        body="Benchmarking compares a finished ranking against baselines over the same peptides. Start a run first."
        action={
          <Button asChild>
            <Link href="/">Start a run</Link>
          </Button>
        }
      />
    );
  }

  if (!report) {
    return (
      <div className="mx-auto max-w-2xl space-y-6">
        <div>
          <h1 className="text-2xl font-semibold tracking-tight">Benchmark</h1>
          <p className="mt-1 text-sm leading-relaxed text-muted-foreground">
            Compare this ranking against the built-in baselines over the same
            peptides.
          </p>
        </div>

        <Card className="border-dashed">
          <CardHeader>
            <CardTitle className="text-base">
              No experimentally validated labels are loaded
            </CardTitle>
            <CardDescription className="leading-relaxed">
              This screen stays empty until a real assay-label file is loaded.
              Generating owned T-cell assay labels is the current work, and
              until those exist there is no honest number to put here. Nothing
              on this screen is ever sample, placeholder, or illustrative data.
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <p className="text-sm leading-relaxed text-muted-foreground">
              A label file is a TSV of peptide, allele and assay call. The
              repository ships one built from the published ELISPOT table of
              Ott et al. 2017 at{" "}
              <code className="font-mono text-[13px]">
                data/benchmarks/ott2017/validated.tsv
              </code>
              , with provenance in{" "}
              <code className="font-mono text-[13px]">SOURCES.md</code>.
            </p>
            <Button asChild variant="secondary" disabled={busy}>
              <label className="cursor-pointer">
                <Upload className="size-4" />
                Load a label file
                <input
                  type="file"
                  accept=".tsv,.csv,.txt"
                  className="sr-only"
                  onChange={(event) =>
                    load(event.target.files?.[0] ?? null)
                  }
                />
              </label>
            </Button>
          </CardContent>
        </Card>
      </div>
    );
  }

  const ranks = report.positive_ranks;
  const methods = Object.keys(ranks);
  const distribution = methods.map((method) => ({
    method,
    ranks: ranks[method] ?? [],
  }));

  return (
    <div className="space-y-6">
      <div className="flex flex-wrap items-end justify-between gap-4">
        <div>
          <h1 className="text-2xl font-semibold tracking-tight">Benchmark</h1>
          <p className="mt-1 text-sm text-muted-foreground">
            {report.source.filename} · {report.source.assayed} assayed ·{" "}
            {report.source.positives} positive · {report.source.matched_to_run}{" "}
            matched to this run
          </p>
        </div>
        <div className="flex flex-wrap gap-1.5">
          {Object.entries(report.source.assays ?? {}).map(([assay, count]) => (
            <Badge key={assay} variant="secondary" className="font-mono">
              {assay} × {count}
            </Badge>
          ))}
        </div>
      </div>

      <Card>
        <CardHeader>
          <CardTitle className="text-base">
            Retrieval over the assayed intersection
          </CardTitle>
          <CardDescription>
            Restricted to the {report.pool_size} peptides that were both scored
            by this run and assayed, so no method wins by nominating an easier
            set.
          </CardDescription>
        </CardHeader>
        <CardContent className="px-0">
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead>method</TableHead>
                <TableHead className="text-right">recall@10</TableHead>
                <TableHead className="text-right">recall@20</TableHead>
                <TableHead className="text-right">precision@10</TableHead>
                <TableHead className="text-right">best rank</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {report.retrieval.map((summary) => (
                <TableRow key={summary.method}>
                  <TableCell className="font-mono">{summary.method}</TableCell>
                  <TableCell className="text-right font-mono tabular-nums">
                    {percent(summary.recall_at_10)}
                  </TableCell>
                  <TableCell className="text-right font-mono tabular-nums">
                    {percent(summary.recall_at_20)}
                  </TableCell>
                  <TableCell className="text-right font-mono tabular-nums">
                    {percent(summary.precision_at_10)}
                  </TableCell>
                  <TableCell className="text-right font-mono tabular-nums">
                    {summary.best_positive_rank ?? "—"}
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </CardContent>
      </Card>

      <Card>
        <CardHeader>
          <CardTitle className="text-base">Validation rate at k</CardTitle>
          <CardDescription>
            Share of the top k that came back positive in the assay.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-6">
          <div className="flex items-center gap-4">
            <span className="w-16 font-mono text-sm tabular-nums">k = {k}</span>
            <Slider
              value={[k]}
              min={1}
              max={Math.max(report.pool_size, 1)}
              step={1}
              onValueChange={([value]) => setK(value)}
              className="max-w-sm"
            />
          </div>

          <div className="h-[260px]">
            <ResponsiveContainer width="100%" height="100%">
              <BarChart
                data={report.validation.map((item) => ({
                  method: item.method,
                  rate: item.validation_rate,
                }))}
                margin={{ top: 4, right: 8, bottom: 4, left: 0 }}
              >
                <XAxis
                  dataKey="method"
                  tick={{ fontSize: 11, fontFamily: "var(--font-mono)" }}
                  tickLine={false}
                  axisLine={{ stroke: "var(--border)" }}
                />
                <YAxis
                  tickFormatter={(value: number) => `${Math.round(value * 100)}%`}
                  tick={{ fontSize: 11 }}
                  tickLine={false}
                  axisLine={false}
                />
                <ChartTooltip
                  cursor={{ fill: "var(--accent)", opacity: 0.4 }}
                  contentStyle={{
                    background: "var(--popover)",
                    border: "1px solid var(--border)",
                    borderRadius: 8,
                    fontSize: 12,
                  }}
                  formatter={(value) => percent(Number(value))}
                />
                <Bar dataKey="rate" radius={2} barSize={36}>
                  {report.validation.map((item) => (
                    <Cell
                      key={item.method}
                      fill={
                        item.method === "neoantigene"
                          ? "var(--chart-1)"
                          : "var(--muted-foreground)"
                      }
                    />
                  ))}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </div>
        </CardContent>
      </Card>

      <Card>
        <CardHeader>
          <CardTitle className="text-base">
            Where the validated peptides landed
          </CardTitle>
          <CardDescription>
            Each mark is one experimentally positive peptide, at its rank.
            Further left is better.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          {distribution.map((item) => (
            <div key={item.method} className="space-y-1.5">
              <div className="flex items-baseline justify-between">
                <span className="font-mono text-sm">{item.method}</span>
                <span className="text-xs text-muted-foreground">
                  {item.ranks.length} positives
                </span>
              </div>
              <div className="relative h-7 rounded-md border bg-muted/40">
                {item.ranks.map((rank, index) => (
                  <span
                    key={`${rank}-${index}`}
                    title={`rank ${rank}`}
                    className="absolute top-1 h-5 w-[3px] rounded-full bg-primary"
                    style={{
                      left: `calc(${(rank / Math.max(report.pool_size, 1)) * 100}% - 1.5px)`,
                    }}
                  />
                ))}
              </div>
            </div>
          ))}
          <p className="text-xs text-muted-foreground">
            Scale runs from rank 1 to {report.pool_size}.
          </p>
        </CardContent>
      </Card>
    </div>
  );
}

function percent(value: number) {
  return `${(value * 100).toFixed(1)}%`;
}
