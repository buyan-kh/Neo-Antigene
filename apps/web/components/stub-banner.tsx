"use client";

import { AlertTriangle } from "lucide-react";
import { useEffect, useState } from "react";

import { useRun } from "@/components/run-context";
import { api, type ReferenceInfo, type RunState } from "@/lib/api";

/**
 * Persistent warning whenever the bundled two-protein stub is in play.
 *
 * The self-peptide gate and the nearest-self feature are meaningless against
 * a stub, so a run made with it must never be presented as a real result.
 * This sits in the shell rather than on one screen so the caveat cannot be
 * scrolled past or navigated away from.
 */
export function StubProteomeBanner() {
  const { runId } = useRun();
  const [reference, setReference] = useState<ReferenceInfo | null>(null);
  const [run, setRun] = useState<RunState | null>(null);

  useEffect(() => {
    api.reference().then(setReference).catch(() => setReference(null));
  }, []);

  // A real proteome being installed is not the same as this run having used
  // one -- a manifest can point anywhere. Usage wins over installation.
  useEffect(() => {
    if (!runId) {
      setRun(null);
      return;
    }
    const poll = () => api.run(runId).then(setRun).catch(() => setRun(null));
    void poll();
    const timer = setInterval(poll, 3000);
    return () => clearInterval(timer);
  }, [runId]);

  if (!reference) return null;

  const installed =
    reference.cached_path !== null && reference.cached_is_stub === false;
  const runUsedStub = run?.warnings?.proteome_is_stub ?? false;
  if (installed && !runUsedStub) return null;

  if (runUsedStub) {
    return (
      <div className="border-b border-amber-500/30 bg-amber-500/10">
        <div className="mx-auto flex max-w-[1600px] items-start gap-3 px-6 py-2.5">
          <AlertTriangle className="mt-0.5 size-4 shrink-0 text-amber-600 dark:text-amber-500" />
          <p className="text-sm leading-relaxed text-amber-900 dark:text-amber-200">
            <span className="font-medium">
              This run used a stub reference proteome
            </span>{" "}
            ({run?.warnings?.proteome_proteins ?? 0} proteins). The self-peptide
            gate and the self-dissimilarity feature are meaningless here, so
            these results must not be read as a real ranking. Point the
            manifest at a full Ensembl release to fix it.
          </p>
        </div>
      </div>
    );
  }

  return (
    <div className="border-b border-amber-500/30 bg-amber-500/10">
      <div className="mx-auto flex max-w-[1600px] items-start gap-3 px-6 py-2.5">
        <AlertTriangle className="mt-0.5 size-4 shrink-0 text-amber-600 dark:text-amber-500" />
        <p className="text-sm leading-relaxed text-amber-900 dark:text-amber-200">
          <span className="font-medium">
            No full reference proteome is installed.
          </span>{" "}
          Runs fall back to the bundled{" "}
          <code className="font-mono text-[13px]">proteome.mini.fa</code> stub
          ({reference.bundled_stub_proteins} proteins), which makes the
          self-peptide gate and the self-dissimilarity feature meaningless —
          nearly every peptide will look foreign. Install the real one with{" "}
          <code className="font-mono text-[13px]">
            {reference.fetch_command}
          </code>
          .
        </p>
      </div>
    </div>
  );
}
