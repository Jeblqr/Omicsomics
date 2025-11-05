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
npm run dev

# Run lint checks
npm run lint
```

The frontend talks to the FastAPI backend (see `../backend`). Adjust API base URLs and authentication settings in `src/lib` (to be added).
