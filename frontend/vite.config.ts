import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// See https://vitejs.dev/config/
export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    host: "0.0.0.0",
    allowedHosts: true,
    // allowedHosts: [
    //   "localhost", // Allow all subdomains
    // ],
    proxy: {
      '/api': {
        target: process.env.VITE_API_TARGET || 'http://localhost:8001',
        changeOrigin: true,
        secure: false,
        rewrite: (path) => path,
      }
    }
  },
  preview: {
    port: 4173,
    host: "0.0.0.0",
  },
});
