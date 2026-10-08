# Demo apps

A screen-shareable view of the ranking engine. **This is a demo, not a
product**: no auth, no accounts, no database, no multi-tenancy, no billing.
Run state lives in memory and disappears when the API restarts.

```
apps/api/   FastAPI wrapper around the installed `neoantigene` package
apps/web/   Next.js App Router UI that renders what the API returns
```

Both depend on `src/neoantigene/`. The package depends on neither, and nothing
in `src/neoantigene/` was restructured to accommodate them.

## One command

```bash
make demo
```

That runs [`scripts/dev.sh`](../scripts/dev.sh), which starts the API on
`:8000` and the UI on `:3000`, installs web dependencies on first run, and
shuts both down together on Ctrl-C.

From a clean checkout the full sequence is:

```bash
make setup   # deps + mhcflurry models + Ensembl proteome
make demo
```

`make setup` is the only manual data preparation, and it is the documented
`mhcflurry-downloads fetch` plus the proteome fetch. The bundled example at
`data/examples/sample.yaml` then runs with one click from cold start.

## What the four screens do

**Input** loads the bundled example in one click, or takes uploaded variants
TSV / VEP VCF / expression files with HLA entered as tags. Validation reuses
the library's own `Sample` model, and its errors are attached to the offending
field rather than dumped as a traceback. The backend selector shows a
persistent destructive banner for development backends, matching the CLI's
guard.

**Run** polls the run and prints the pipeline's own progress as it arrives:
variants read, somatic filtering summary, peptides generated, pairs scored,
pairs gated, shortlist size. Failures show the real exception message in a
readable card. Nothing is swallowed into "something went wrong".

**Results** is the main screen. A sortable table over `io.writers.COLUMNS`
with column visibility control, peptide sequences in monospace with the
mutated residue in the accent colour and the wild-type stacked underneath,
filters for gene, allele and gated status, and muted badges naming the gate
that excluded a row. Clicking a row opens a side sheet with the per-feature
breakdown from `scoring.rank.contributions` as a horizontal diverging bar
chart of signed contributions, bias included. Downloads are served from the
library's existing writers.

**Benchmark** renders only when a real validated-label file is loaded. With no
labels it shows a designed empty state stating plainly that no experimentally
validated labels are loaded and that generating owned labels is the current
work. It never displays sample, placeholder or illustrative numbers. If an
uploaded label file has no overlap with the run, the API refuses and says so
rather than reporting a metric over an empty intersection.

## Honesty constraints wired into the UI

- A **persistent banner** appears whenever the bundled two-protein
  `proteome.mini.fa` stub is in use, either because no full proteome is
  installed or because the current run's manifest pointed at the stub. The
  self-peptide gate and the self-dissimilarity feature are meaningless against
  it, and the UI says so rather than leaving the caveat in the README.
- Every number on screen comes from a real pipeline run. There is no fixture
  data, no seeded demo run, and no fallback values.
- Empty states say what is empty and why.

## Types

The UI does not hand-write response types. `lib/api-types.ts` is generated
from the API's OpenAPI schema, which FastAPI derives from the Python pydantic
models. Regenerate after changing an API model:

```bash
uv run --extra api python -c "import json,sys; sys.path.insert(0,'apps/api'); \
  from neoantigene_api.main import app; \
  json.dump(app.openapi(), open('apps/web/openapi.json','w'), indent=2)"
cd apps/web && npx openapi-typescript openapi.json -o lib/api-types.ts
```

A library change that alters a response shape then surfaces as a TypeScript
error instead of a silently wrong table.

## Architecture notes

**The API wraps, it does not reimplement.** Ranking, scoring, peptide
generation, writers and metrics are all called from `neoantigene`. No ranking
logic exists in TypeScript.

**Progress without modifying the library.** `apps/api/neoantigene_api/progress.py`
attaches a `logging.Handler` to the `neoantigene` logger and routes records to
the owning run using a `threading.local`. No callback parameter was added to
the library to make the UI work.

**Runs are async and in memory.** `POST /api/runs` returns a run id
immediately; `GET /api/runs/{id}` returns status, progress events and the
pipeline's `RunReport`. One run executes at a time.

**Python deps are optional.** The API's dependencies live in the `api` extra
in `pyproject.toml`, so a plain `uv sync` stays lean and CI is unaffected. The
API's tests live under `apps/api/tests/` and are outside the root `testpaths`
for the same reason.

## Known rough edges

- Run state is in memory, so restarting the API loses history. Deliberate.
- One concurrent run. A second request queues.
- Results filtering and sorting happen client-side over the returned page,
  which is fine at shortlist size and would not be at genome scale.
- Mobile is not a target; the layout is built for a laptop screen share.

## Deploy on Vercel

One Vercel project, two services, configured in the repository-root
`vercel.json`. Leave the project's Root Directory at the repository root.
Do not point it at `apps/web`.

| Service | Root | Public path |
| --- | --- | --- |
| `web` | `apps/web` | everything except `/api` |
| `api` | repository root (`app:app`) | `/api` and `/api/*` |

The browser calls `/api/...` on the same host as the UI. No
`NEXT_PUBLIC_API_BASE` is required for that deploy. Set it only when the UI
should call an API on a different host. `make demo` still sets it to the
local uvicorn origin.

```bash
npx vercel login
npx vercel link    # from the repository root
npx vercel dev     # both services
npx vercel --prod
```

The API process is still in-memory, and a ranking run can outlive one
function instance. MHCflurry's downloaded model weights are not part of the
deploy; until they are present the API reports that presentation models are
missing. The static page in `share/` is not one of these services.
