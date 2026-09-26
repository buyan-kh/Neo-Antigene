/**
 * Thin fetch wrapper over the FastAPI app.
 *
 * Every type here is re-exported from `api-types.ts`, which is generated from
 * the API's OpenAPI schema, which is in turn generated from the Python
 * pydantic models. Nothing about the response shape is written twice, so a
 * change in the library surfaces as a TypeScript error rather than as a
 * silently wrong table.
 */
import type { components } from "./api-types";

type Schemas = components["schemas"];

export type BackendInfo = Schemas["BackendInfo"];
export type ReferenceInfo = Schemas["ReferenceInfo"];
export type RunState = Schemas["RunState"];
export type RunCreated = Schemas["RunCreated"];
export type RunRequest = Schemas["RunRequest"];
export type SampleInput = Schemas["SampleInput"];
export type ResultsPage = Schemas["ResultsPage"];
export type Explanation = Schemas["Explanation"];
export type Contribution = Schemas["Contribution"];
export type BenchmarkReport = Schemas["BenchmarkReport"];
export type Workspace = Schemas["Workspace"];
/**
 * Not generated: the manifest-validation errors travel inside `detail`, which
 * FastAPI types as unknown, so the shape is declared once here to match
 * `neoantigene_api.schemas.FieldError`.
 */
export type FieldError = { field: string; message: string };
export type ProgressEvent = Schemas["ProgressEvent"];
export type RunWarnings = Schemas["RunWarnings"];
export type RunReport = Schemas["RunReport"];

/** `/api/example` returns the parsed manifest plus the file names it points at. */
export type ExampleSample = {
  path: string;
  sample_id: string;
  cancer_type: string | null;
  tumor_purity: number | null;
  hla: string[];
  variant_tsv: string | null;
  somatic_vcf: string | null;
  expression_tsv: string | null;
  normal_expression_tsv: string | null;
  proteome_fasta: string | null;
};

export const API_BASE =
  process.env.NEXT_PUBLIC_API_BASE ?? "http://127.0.0.1:8000";

/**
 * An error carrying whatever the API actually said.
 *
 * The API reports manifest validation failures as per-field errors so the UI
 * can attach them to the offending input. Collapsing that into a generic
 * message would throw away the only part the user can act on.
 */
export class ApiError extends Error {
  readonly status: number;
  readonly fields: FieldError[];

  constructor(status: number, message: string, fields: FieldError[] = []) {
    super(message);
    this.name = "ApiError";
    this.status = status;
    this.fields = fields;
  }
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  let response: Response;
  try {
    response = await fetch(`${API_BASE}${path}`, {
      ...init,
      headers: {
        ...(init?.body instanceof FormData
          ? {}
          : { "Content-Type": "application/json" }),
        ...init?.headers,
      },
    });
  } catch {
    throw new ApiError(
      0,
      `Cannot reach the API at ${API_BASE}. Is it running? Try \`make api\`.`,
    );
  }

  if (!response.ok) {
    throw await toApiError(response);
  }
  return (await response.json()) as T;
}

async function toApiError(response: Response): Promise<ApiError> {
  let detail: unknown;
  try {
    detail = (await response.json())?.detail;
  } catch {
    return new ApiError(response.status, await safeText(response));
  }

  // Manifest validation returns a list of {field, message}.
  if (Array.isArray(detail)) {
    const fields = detail.filter(
      (item): item is FieldError =>
        typeof item === "object" && item !== null && "field" in item,
    );
    if (fields.length > 0) {
      return new ApiError(response.status, "The sample manifest is invalid.", fields);
    }
    return new ApiError(response.status, JSON.stringify(detail));
  }
  return new ApiError(
    response.status,
    typeof detail === "string" ? detail : `Request failed (${response.status})`,
  );
}

async function safeText(response: Response): Promise<string> {
  try {
    return (await response.text()) || `Request failed (${response.status})`;
  } catch {
    return `Request failed (${response.status})`;
  }
}

export const api = {
  backends: () => request<BackendInfo>("/api/backends"),
  reference: () => request<ReferenceInfo>("/api/reference"),
  example: () => request<ExampleSample>("/api/example"),

  upload: (files: FormData) =>
    request<Workspace>("/api/uploads", { method: "POST", body: files }),

  startRun: (body: RunRequest) =>
    request<RunCreated>("/api/runs", {
      method: "POST",
      body: JSON.stringify(body),
    }),

  run: (runId: string) => request<RunState>(`/api/runs/${runId}`),

  results: (runId: string, params?: Record<string, string | undefined>) => {
    const query = new URLSearchParams(
      Object.entries(params ?? {}).filter(
        (entry): entry is [string, string] => Boolean(entry[1]),
      ),
    );
    const suffix = query.toString() ? `?${query}` : "";
    return request<ResultsPage>(`/api/runs/${runId}/results${suffix}`);
  },

  explain: (runId: string, peptide: string, allele: string) =>
    request<Explanation>(
      `/api/runs/${runId}/explain?peptide=${encodeURIComponent(
        peptide,
      )}&allele=${encodeURIComponent(allele)}`,
    ),

  benchmark: (runId: string, labels: FormData) =>
    request<BenchmarkReport>(`/api/runs/${runId}/benchmark`, {
      method: "POST",
      body: labels,
    }),

  downloadUrl: (runId: string, artifact: "ranked.tsv" | "features.json") =>
    `${API_BASE}/api/runs/${runId}/download/${artifact}`,
};
