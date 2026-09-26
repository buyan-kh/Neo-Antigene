#!/usr/bin/env bash
# Bring up the API and the web UI together, and shut both down on Ctrl-C.
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$PWD"

API_PORT="${API_PORT:-8000}"
WEB_PORT="${WEB_PORT:-3000}"

if [[ ! -d apps/web/node_modules ]]; then
  echo "installing web dependencies ..."
  (cd apps/web && npm install)
fi

cleanup() {
  # Kill the whole process group so uvicorn's reloader and next's workers go too.
  trap - EXIT INT TERM
  [[ -n "${API_PID:-}" ]] && kill -- "-${API_PID}" 2>/dev/null || true
  [[ -n "${WEB_PID:-}" ]] && kill -- "-${WEB_PID}" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

echo "API  -> http://127.0.0.1:${API_PORT}"
set -m
uv run --extra api --extra presentation uvicorn neoantigene_api.main:app \
  --app-dir apps/api --port "${API_PORT}" &
API_PID=$!

echo "Web  -> http://127.0.0.1:${WEB_PORT}"
(
  cd "${ROOT}/apps/web"
  NEXT_PUBLIC_API_BASE="http://127.0.0.1:${API_PORT}" npm run dev -- --port "${WEB_PORT}"
) &
WEB_PID=$!
set +m

# Plain `wait`, not `wait -n`: macOS still ships bash 3.2, where -n is a
# syntax error and the script would exit immediately, killing both servers.
wait "${API_PID}" "${WEB_PID}"
