"use client";

import {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useState,
} from "react";

import { ApiError, api, type RunState } from "@/lib/api";

type RunContextValue = {
  runId: string | null;
  run: RunState | null;
  /** True until the stored run id has been checked against the API. */
  restoring: boolean;
  setRunId: (id: string | null) => void;
  refresh: () => Promise<RunState | null>;
};

const RunContext = createContext<RunContextValue | null>(null);

const STORAGE_KEY = "neoantigene.runId";

/**
 * The current run, shared across screens.
 *
 * Only the run *id* is remembered, in `sessionStorage`, and every screen
 * re-reads the run itself from the API. The run state lives in the API
 * process's memory, so this adds no persistence the backend does not already
 * have: close the tab or restart the API and the run is gone. Without this,
 * clicking from Results to Benchmark and back discarded the run, which is
 * exactly what happens mid-demo.
 */
export function RunProvider({ children }: { children: React.ReactNode }) {
  const [runId, setRunId] = useState<string | null>(null);
  const [run, setRun] = useState<RunState | null>(null);
  const [restoring, setRestoring] = useState(true);

  const refresh = useCallback(async () => {
    if (!runId) return null;
    const next = await api.run(runId);
    setRun(next);
    return next;
  }, [runId]);

  // Restore on first mount. A 404 means the API restarted or forgot the run,
  // so the stale id is dropped rather than left to fail on every screen.
  useEffect(() => {
    let cancelled = false;
    const stored = window.sessionStorage.getItem(STORAGE_KEY);
    // Every state update below happens in a promise callback rather than in
    // the effect body, including the no-stored-id case, so this never sets
    // state synchronously during the effect.
    const restore = stored ? api.run(stored) : Promise.resolve(null);
    restore
      .then((state) => {
        if (cancelled || !stored || !state) return;
        setRunId(stored);
        setRun(state);
      })
      .catch((error: unknown) => {
        if (cancelled) return;
        if (!(error instanceof ApiError) || error.status !== 404) {
          console.error("could not restore the previous run", error);
        }
        window.sessionStorage.removeItem(STORAGE_KEY);
      })
      .finally(() => {
        if (!cancelled) setRestoring(false);
      });
    return () => {
      cancelled = true;
    };
  }, []);

  const value = useMemo<RunContextValue>(
    () => ({
      runId,
      run,
      refresh,
      restoring,
      setRunId: (id) => {
        setRunId(id);
        setRun(null);
        if (id) window.sessionStorage.setItem(STORAGE_KEY, id);
        else window.sessionStorage.removeItem(STORAGE_KEY);
      },
    }),
    [runId, run, refresh, restoring],
  );

  return <RunContext.Provider value={value}>{children}</RunContext.Provider>;
}

export function useRun() {
  const context = useContext(RunContext);
  if (!context) throw new Error("useRun must be used inside RunProvider");
  return context;
}
