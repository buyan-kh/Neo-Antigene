"use client";

import Link from "next/link";
import { useEffect, useState } from "react";
import { AlertTriangle, Check, Loader2 } from "lucide-react";

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
import { Separator } from "@/components/ui/separator";
import { Skeleton } from "@/components/ui/skeleton";
import { EmptyState } from "@/components/empty-state";
import type { RunReport, RunState } from "@/lib/api";

const POLL_MS = 700;

export default function RunPage() {
  const { runId, run, refresh, restoring } = useRun();
  const [failedToPoll, setFailedToPoll] = useState<string | null>(null);

  useEffect(() => {
    if (!runId) return;
    let live = true;

    const tick = async () => {
      try {
        const state = await refresh();
        setFailedToPoll(null);
        if (!live) return;
        if (state && (state.status === "running" || state.status === "queued")) {
          setTimeout(tick, POLL_MS);
        }
      } catch (error) {
        if (live) setFailedToPoll((error as Error).message);
      }
    };
    void tick();
    return () => {
      live = false;
    };
  }, [runId, refresh]);

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
        title="No run started"
        body="Nothing is running yet. Start one from the Input screen and its progress will appear here as the pipeline reports it."
        action={
          <Button asChild>
            <Link href="/">Go to Input</Link>
          </Button>
        }
      />
    );
  }

  if (!run) {
    return (
      <div className="mx-auto max-w-3xl space-y-3">
        <Skeleton className="h-8 w-56" />
        <Skeleton className="h-40 w-full" />
      </div>
    );
  }

  const events = run.events ?? [];
  const done = run.status === "succeeded";
  const failed = run.status === "failed";

  return (
    <div className="mx-auto max-w-3xl space-y-6">
      <div className="flex items-start justify-between gap-6">
        <div>
          <h1 className="text-2xl font-semibold tracking-tight">
            {run.sample_id ?? "Run"}
          </h1>
          <p className="mt-1 font-mono text-xs text-muted-foreground">
            {run.run_id}
          </p>
        </div>
        <StatusBadge status={run.status} />
      </div>

      {failedToPoll && (
        <Card className="border-destructive/40">
          <CardContent className="pt-6 text-sm text-destructive">
            {failedToPoll}
          </CardContent>
        </Card>
      )}

      {failed && (
        <Card className="border-destructive/40">
          <CardHeader>
            <CardTitle className="text-base text-destructive">
              {run.error_type ?? "Run failed"}
            </CardTitle>
            <CardDescription>
              The pipeline raised this. It is the real message, not a summary.
            </CardDescription>
          </CardHeader>
          <CardContent>
            <pre className="overflow-x-auto whitespace-pre-wrap rounded-md bg-muted p-3 font-mono text-[13px] leading-relaxed">
              {run.error}
            </pre>
          </CardContent>
        </Card>
      )}

      <Warnings run={run} />

      <Card>
        <CardHeader>
          <CardTitle className="text-base">Progress</CardTitle>
          <CardDescription>
            Reported by the pipeline itself, as each stage completes.
          </CardDescription>
        </CardHeader>
        <CardContent>
          {events.length === 0 ? (
            <p className="text-sm text-muted-foreground">
              Waiting for the first stage to report…
            </p>
          ) : (
            <ol className="space-y-2.5">
              {events.map((event, index) => (
                <li key={index} className="flex items-start gap-3 text-sm">
                  <span className="mt-0.5">
                    {event.level === "WARNING" ? (
                      <AlertTriangle className="size-4 text-amber-600 dark:text-amber-500" />
                    ) : (
                      <Check className="size-4 text-muted-foreground" />
                    )}
                  </span>
                  <span
                    className={
                      event.level === "WARNING"
                        ? "leading-relaxed text-amber-700 dark:text-amber-400"
                        : "leading-relaxed"
                    }
                  >
                    {event.message}
                  </span>
                </li>
              ))}
              {(run.status === "running" || run.status === "queued") && (
                <li className="flex items-center gap-3 text-sm text-muted-foreground">
                  <Loader2 className="size-4 animate-spin" />
                  Working…
                </li>
              )}
            </ol>
          )}
        </CardContent>
      </Card>

      {run.report && <ReportCard report={run.report} />}

      {done && (
        <div className="flex justify-end">
          <Button asChild>
            <Link href="/results">View results</Link>
          </Button>
        </div>
      )}
    </div>
  );
}

function StatusBadge({ status }: { status: string }) {
  const variant =
    status === "failed"
      ? "destructive"
      : status === "succeeded"
        ? "default"
        : "secondary";
  return (
    <Badge variant={variant} className="capitalize">
      {status}
    </Badge>
  );
}

function Warnings({ run }: { run: RunState }) {
  const warnings = run.warnings;
  const items = [
    warnings?.development_detail,
    warnings?.proteome_detail,
    warnings?.expression_detail,
  ].filter(Boolean) as string[];

  if (items.length === 0) return null;

  return (
    <Card className="border-amber-500/40 bg-amber-500/5">
      <CardHeader>
        <CardTitle className="text-base">Caveats for this run</CardTitle>
      </CardHeader>
      <CardContent>
        <ul className="space-y-2">
          {items.map((item) => (
            <li key={item} className="flex gap-2.5 text-sm leading-relaxed">
              <AlertTriangle className="mt-0.5 size-4 shrink-0 text-amber-600 dark:text-amber-500" />
              {item}
            </li>
          ))}
        </ul>
      </CardContent>
    </Card>
  );
}

function ReportCard({ report }: { report: RunReport }) {
  const rows: [string, number | string][] = [
    ["Variants read", report.variants_read],
    ["Variants kept after filtering", report.variants_kept],
    ["Peptides generated", report.peptides_generated],
    ["Peptide-HLA pairs scored", report.pairs_scored],
    [
      "Pairs gated out",
      Object.values(report.pairs_gated ?? {}).reduce((sum, n) => sum + n, 0),
    ],
    ["Shortlist", report.shortlisted],
    ["Proteome proteins", report.proteome_proteins ?? 0],
  ];

  return (
    <Card>
      <CardHeader>
        <CardTitle className="text-base">Run report</CardTitle>
      </CardHeader>
      <CardContent>
        <dl className="grid grid-cols-2 gap-x-8">
          {rows.map(([label, value], index) => (
            <div key={label}>
              {index > 1 && <Separator className="my-0" />}
              <div className="flex items-baseline justify-between py-1.5 text-sm">
                <dt className="text-muted-foreground">{label}</dt>
                <dd className="font-mono tabular-nums">{value}</dd>
              </div>
            </div>
          ))}
        </dl>
      </CardContent>
    </Card>
  );
}
