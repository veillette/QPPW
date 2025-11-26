import { defineConfig } from "vite";
import { visualizer } from "rollup-plugin-visualizer";

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [
    // Bundle analysis - generates stats.html in the project root
    // Set to false in normal builds, enable with ANALYZE=true npm run build
    process.env.ANALYZE === "true" &&
      visualizer({
        filename: "./dist/stats.html",
        open: false,
        gzipSize: true,
        brotliSize: true,
        template: "treemap", // treemap, sunburst, network
      }),
  ].filter(Boolean),
  // So the build can be served from an arbitrary path
  base: "./",
  build: {
    chunkSizeWarningLimit: 2000,
    // Enable terser for better minification and tree shaking
    minify: "terser",
    terserOptions: {
      compress: {
        // Remove unused code
        dead_code: true,
        drop_console: false, // Keep console for debugging
        drop_debugger: true,
        passes: 2, // Multiple passes for better optimization
        unused: true, // Remove unused variables and functions
      },
      mangle: {
        safari10: true,
      },
      format: {
        comments: false, // Remove comments
      },
    },
    rollupOptions: {
      // Enable tree shaking - Vite will respect the sideEffects field in package.json
      treeshake: {
        preset: "recommended",
        // Vite will automatically handle side effects based on package.json
      },
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
