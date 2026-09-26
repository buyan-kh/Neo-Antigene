import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  // Without this, visiting the dev server as 127.0.0.1 rather than localhost
  // makes Next refuse to serve its own client assets, so React never hydrates:
  // the page renders but every button is inert. Both hostnames are listed
  // because either is a reasonable thing to type during a screen share.
  allowedDevOrigins: ["127.0.0.1", "localhost"],
};

export default nextConfig;
