# Omicsomics Frontend

This directory contains the web UI for the Omicsomics unified omics platform. It is built with React, TypeScript, and Vite.

## Features (planned)

- Dashboard for global metrics and alerts
- Project and sample management views
- Pipeline run monitoring and logs
- Data catalog with metadata search and QC summaries
- Configurable visualizations (UMAP, volcano plots, network diagrams)

## Getting Started

```bash
# Install dependencies
npm install

# Start development server
npm run dev -- --host

# Run lint checks
npm run lint
```

# Build production bundle

npm run build

The frontend talks to the FastAPI backend (see `../backend`). Adjust API base URLs and authentication settings in `src/lib` (to be added).
Environment variables (see `infrastructure/docker-compose.yml` for defaults):

- `VITE_API_BASE_URL` — defaults to `http://localhost:8001`
- Additional secrets (auth, analytics) can be added via `.env.local`

The frontend talks to the FastAPI backend (see `../backend`). REST endpoints are currently served under `/api/v1`.
