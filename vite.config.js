import { defineConfig } from "vite";

// https://vitejs.dev/config/
export default defineConfig({
  // So the build can be served from an arbitrary path
  base: "./",
  build: {
    chunkSizeWarningLimit: 2000,
    rollupOptions: {
      output: {
        manualChunks: (id) => {
          // Vendor libraries - split large dependencies
          if (id.includes("node_modules")) {
            if (id.includes("scenerystack")) {
              return "vendor-scenery";
            }
            if (id.includes("katex")) {
              return "vendor-katex";
            }
            if (id.includes("three")) {
              return "vendor-three";
            }
            // Other vendor dependencies
            return "vendor-other";
          }

          // Core quantum solvers (shared across screens)
          if (
            id.includes("/common/model/") &&
            (id.includes("Solver") ||
              id.includes("QuantumBoundState") ||
              id.includes("Schrodinger1D"))
          ) {
            return "quantum-core";
          }

          // Analytical solutions
          if (id.includes("/analytical-solutions/")) {
            return "analytical";
          }

          // Common view components
          if (id.includes("/common/view/")) {
            return "common-view";
          }

          // Screen-specific chunks
          if (id.includes("/intro/")) {
            return "screen-intro";
          }
          if (id.includes("/one-well/")) {
            return "screen-one-well";
          }
          if (id.includes("/two-wells/")) {
            return "screen-two-wells";
          }
          if (id.includes("/many-wells/")) {
            return "screen-many-wells";
          }
        },
      },
    },
  },
});
