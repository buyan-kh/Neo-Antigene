import type { NextConfig } from "next";

// Decided when the client bundle is built. `VERCEL=1` covers `vercel dev`
// and a Vercel deployment: the API is the other service on this same host.
// `make demo` sets NEXT_PUBLIC_API_BASE to the local uvicorn origin, and a
// plain `npm run dev` keeps the laptop default.
const configuredApiBase = process.env.NEXT_PUBLIC_API_BASE;
const apiBase =
  configuredApiBase !== undefined
    ? configuredApiBase
    : process.env.VERCEL === "1"
      ? ""
      : "http://127.0.0.1:8000";

const nextConfig: NextConfig = {
  // Without this, visiting the dev server as 127.0.0.1 rather than localhost
  // makes Next refuse to serve its own client assets, so React never hydrates:
  // the page renders but every button is inert. Both hostnames are listed
  // because either is a reasonable thing to type during a screen share.
  allowedDevOrigins: ["127.0.0.1", "localhost"],
  env: {
    NEXT_PUBLIC_API_BASE: apiBase,
  },
};

export default nextConfig;
