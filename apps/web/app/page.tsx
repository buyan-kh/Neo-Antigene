"use client";

import { useRouter } from "next/navigation";
import { useEffect, useState } from "react";
import { toast } from "sonner";
import { Loader2, Upload, X } from "lucide-react";

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
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import { Separator } from "@/components/ui/separator";
import { Skeleton } from "@/components/ui/skeleton";
import {
  ApiError,
  api,
  type BackendInfo,
  type ExampleSample,
  type FieldError,
} from "@/lib/api";

type Mode = "example" | "upload";

export default function InputPage() {
  const router = useRouter();
  const { setRunId } = useRun();

  const [backends, setBackends] = useState<BackendInfo | null>(null);
  const [example, setExample] = useState<ExampleSample | null>(null);
  const [mode, setMode] = useState<Mode>("example");
  const [backend, setBackend] = useState<string>("");
  const [starting, setStarting] = useState(false);
  const [errors, setErrors] = useState<FieldError[]>([]);

  const [sampleId, setSampleId] = useState("");
  const [cancerType, setCancerType] = useState("");
  const [purity, setPurity] = useState("");
  const [hla, setHla] = useState<string[]>([]);
  const [hlaDraft, setHlaDraft] = useState("");
  const [workspace, setWorkspace] = useState<string | null>(null);
  const [files, setFiles] = useState<Record<string, string>>({});

  useEffect(() => {
    api
      .backends()
      .then((info) => {
        setBackends(info);
        setBackend(info.default);
      })
      .catch((error: Error) => toast.error(error.message));
    api.example().then(setExample).catch(() => setExample(null));
  }, []);

  const isDevelopmentBackend = backends?.development.includes(backend) ?? false;
  const fieldError = (name: string) =>
    errors.find((item) => item.field === name)?.message;

  async function onUpload(slot: string, file: File | null) {
    if (!file) return;
    const body = new FormData();
    body.append("files", file);
    try {
      const result = await api.upload(body);
      setWorkspace(result.id);
      setFiles((current) => ({ ...current, [slot]: file.name }));
    } catch (error) {
      toast.error((error as Error).message);
    }
  }

  async function start() {
    setStarting(true);
    setErrors([]);
    try {
      const created = await api.startRun(
        mode === "example"
          ? {
              source: "example",
              backend,
              acknowledge_development_backend: isDevelopmentBackend,
            }
          : {
              source: "upload",
              workspace,
              backend,
              acknowledge_development_backend: isDevelopmentBackend,
              sample: {
                sample_id: sampleId,
                cancer_type: cancerType || null,
                tumor_purity: purity ? Number(purity) : null,
                hla,
                variant_tsv: files.variant_tsv ?? null,
                somatic_vcf: files.somatic_vcf ?? null,
                expression_tsv: files.expression_tsv ?? null,
                normal_expression_tsv: files.normal_expression_tsv ?? null,
                proteome: "ensembl",
              },
            },
      );
      setRunId(created.run_id);
      router.push("/run");
    } catch (error) {
      if (error instanceof ApiError && error.fields.length > 0) {
        setErrors(error.fields);
        toast.error(error.message);
      } else {
        toast.error((error as Error).message);
      }
    } finally {
      setStarting(false);
    }
  }

  const bothSources = Boolean(files.variant_tsv) && Boolean(files.somatic_vcf);

  return (
    <div className="mx-auto max-w-3xl space-y-6">
      <div>
        <h1 className="text-2xl font-semibold tracking-tight">New run</h1>
        <p className="mt-1 text-sm leading-relaxed text-muted-foreground">
          Rank mutant peptides for one patient from somatic variants and a class I
          HLA type.
        </p>
      </div>

      <Card>
        <CardHeader>
          <CardTitle className="text-base">Input</CardTitle>
          <CardDescription>
            Start from the bundled example, or upload your own files.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-2 gap-3">
            <ModeCard
              active={mode === "example"}
              onClick={() => setMode("example")}
              title="Bundled example"
              detail={
                example
                  ? `${example.sample_id} · ${example.hla.length} alleles`
                  : "Loading…"
              }
            />
            <ModeCard
              active={mode === "upload"}
              onClick={() => setMode("upload")}
              title="Upload files"
              detail="Variants TSV or VEP VCF, plus HLA"
            />
          </div>

          {mode === "example" && example && (
            <div className="rounded-lg border bg-muted/40 p-4 text-sm">
              <dl className="grid grid-cols-[auto_1fr] gap-x-6 gap-y-1.5">
                <dt className="text-muted-foreground">Manifest</dt>
                <dd className="font-mono text-[13px]">{example.path}</dd>
                <dt className="text-muted-foreground">Sample</dt>
                <dd className="font-mono text-[13px]">{example.sample_id}</dd>
                <dt className="text-muted-foreground">Variants</dt>
                <dd className="font-mono text-[13px]">
                  {example.variant_tsv ?? example.somatic_vcf ?? "—"}
                </dd>
                <dt className="text-muted-foreground">HLA</dt>
                <dd className="flex flex-wrap gap-1">
                  {example.hla.map((allele) => (
                    <Badge key={allele} variant="secondary" className="font-mono">
                      {allele}
                    </Badge>
                  ))}
                </dd>
              </dl>
            </div>
          )}

          {mode === "upload" && (
            <div className="space-y-4">
              <Field
                label="Sample ID"
                error={fieldError("sample_id")}
                htmlFor="sample_id"
              >
                <Input
                  id="sample_id"
                  value={sampleId}
                  onChange={(event) => setSampleId(event.target.value)}
                  placeholder="PT-001"
                />
              </Field>

              <div className="grid grid-cols-2 gap-4">
                <Field label="Cancer type" htmlFor="cancer_type">
                  <Input
                    id="cancer_type"
                    value={cancerType}
                    onChange={(event) => setCancerType(event.target.value)}
                    placeholder="melanoma"
                  />
                </Field>
                <Field
                  label="Tumor purity"
                  error={fieldError("tumor_purity")}
                  htmlFor="tumor_purity"
                >
                  <Input
                    id="tumor_purity"
                    value={purity}
                    onChange={(event) => setPurity(event.target.value)}
                    placeholder="0.6"
                    inputMode="decimal"
                  />
                </Field>
              </div>

              <Field label="Class I HLA" error={fieldError("hla")} htmlFor="hla">
                <div className="space-y-2">
                  <div className="flex gap-2">
                    <Input
                      id="hla"
                      value={hlaDraft}
                      placeholder="HLA-A*02:01"
                      className="font-mono"
                      onChange={(event) => setHlaDraft(event.target.value)}
                      onKeyDown={(event) => {
                        if (event.key !== "Enter") return;
                        event.preventDefault();
                        const value = hlaDraft.trim();
                        if (value && !hla.includes(value)) setHla([...hla, value]);
                        setHlaDraft("");
                      }}
                    />
                    <Button
                      type="button"
                      variant="secondary"
                      onClick={() => {
                        const value = hlaDraft.trim();
                        if (value && !hla.includes(value)) setHla([...hla, value]);
                        setHlaDraft("");
                      }}
                    >
                      Add
                    </Button>
                  </div>
                  {hla.length === 0 ? (
                    <p className="text-xs text-muted-foreground">
                      No alleles yet. Press Enter to add each one.
                    </p>
                  ) : (
                    <div className="flex flex-wrap gap-1.5">
                      {hla.map((allele) => (
                        <Badge key={allele} variant="secondary" className="gap-1 font-mono">
                          {allele}
                          <button
                            type="button"
                            aria-label={`Remove ${allele}`}
                            onClick={() => setHla(hla.filter((a) => a !== allele))}
                            className="opacity-60 transition-opacity duration-150 hover:opacity-100"
                          >
                            <X className="size-3" />
                          </button>
                        </Badge>
                      ))}
                    </div>
                  )}
                </div>
              </Field>

              <Separator />

              <div className="grid gap-3">
                <FileSlot
                  label="Variants TSV"
                  slot="variant_tsv"
                  files={files}
                  onUpload={onUpload}
                />
                <FileSlot
                  label="VEP VCF"
                  slot="somatic_vcf"
                  files={files}
                  onUpload={onUpload}
                />
                <FileSlot
                  label="Expression TSV (optional)"
                  slot="expression_tsv"
                  files={files}
                  onUpload={onUpload}
                />
                <FileSlot
                  label="Normal expression TSV (optional)"
                  slot="normal_expression_tsv"
                  files={files}
                  onUpload={onUpload}
                />
              </div>

              {bothSources && (
                <p className="text-sm text-destructive">
                  Provide exactly one of variants TSV or VEP VCF, not both.
                </p>
              )}
            </div>
          )}
        </CardContent>
      </Card>

      <Card>
        <CardHeader>
          <CardTitle className="text-base">Presentation backend</CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          {!backends ? (
            <Skeleton className="h-9 w-full" />
          ) : (
            <>
              <Select value={backend} onValueChange={setBackend}>
                <SelectTrigger className="w-full">
                  <SelectValue />
                </SelectTrigger>
                <SelectContent>
                  {backends.available.map((name) => (
                    <SelectItem key={name} value={name}>
                      {name}
                      {backends.development.includes(name) && " (development)"}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>

              {!backends.mhcflurry_ready && backend === "mhcflurry" && (
                <p className="text-sm text-muted-foreground">
                  {backends.mhcflurry_detail}
                </p>
              )}

              {isDevelopmentBackend && (
                <div className="rounded-md border border-destructive/40 bg-destructive/10 p-3">
                  <p className="text-sm font-medium text-destructive">
                    {backend} is a development stand-in, not a predictor.
                  </p>
                  <p className="mt-1 text-sm leading-relaxed text-destructive/90">
                    Its scores are structural placeholders for exercising the
                    pipeline. Any ranking produced with it is meaningless and must
                    not be read as a result.
                  </p>
                </div>
              )}
            </>
          )}
        </CardContent>
      </Card>

      <div className="flex items-center justify-end gap-3">
        <Button onClick={start} disabled={starting || bothSources}>
          {starting && <Loader2 className="size-4 animate-spin" />}
          Start run
        </Button>
      </div>
    </div>
  );
}

function ModeCard({
  active,
  onClick,
  title,
  detail,
}: {
  active: boolean;
  onClick: () => void;
  title: string;
  detail: string;
}) {
  return (
    <button
      type="button"
      onClick={onClick}
      className={`rounded-lg border p-4 text-left transition-colors duration-150 ${
        active ? "border-primary bg-accent/50" : "hover:bg-accent/30"
      }`}
    >
      <div className="text-sm font-medium">{title}</div>
      <div className="mt-0.5 text-xs text-muted-foreground">{detail}</div>
    </button>
  );
}

function Field({
  label,
  error,
  htmlFor,
  children,
}: {
  label: string;
  error?: string;
  htmlFor: string;
  children: React.ReactNode;
}) {
  return (
    <div className="space-y-1.5">
      <Label htmlFor={htmlFor}>{label}</Label>
      {children}
      {error && <p className="text-xs text-destructive">{error}</p>}
    </div>
  );
}

function FileSlot({
  label,
  slot,
  files,
  onUpload,
}: {
  label: string;
  slot: string;
  files: Record<string, string>;
  onUpload: (slot: string, file: File | null) => void;
}) {
  return (
    <div className="flex items-center justify-between gap-4 rounded-md border px-3 py-2">
      <div className="min-w-0">
        <div className="text-sm">{label}</div>
        <div className="truncate font-mono text-xs text-muted-foreground">
          {files[slot] ?? "none"}
        </div>
      </div>
      <Button asChild variant="secondary" size="sm">
        <label className="cursor-pointer">
          <Upload className="size-3.5" />
          Choose
          <input
            type="file"
            className="sr-only"
            onChange={(event) => onUpload(slot, event.target.files?.[0] ?? null)}
          />
        </label>
      </Button>
    </div>
  );
}
