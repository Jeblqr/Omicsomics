# Omicsomics Frontend

This directory contains the web UI for the Omicsomics unified omics platform. It is built with React, TypeScript, and Vite.

## Features

### ✅ Implemented

- **Authentication**: Login/register with JWT token-based auth
- **Projects Manager**: Full CRUD operations for projects
  - Create new projects with name and description
  - Edit existing projects
  - Delete projects
  - Select active project for data/runs operations
- **Runs Manager**: Pipeline run tracking
  - Create new runs linked to projects
  - View run status (pending/running/completed/failed)
  - Track start/finish times
  - Filter runs by selected project
- **Data Catalog**: Encrypted file storage and management
  - Upload files with automatic per-user encryption (HKDF + AES-GCM)
  - View file metadata (size, checksum, upload time)
  - Download decrypted files
  - Files are isolated by project
- **Dashboard**: Global metrics and system overview
- **Settings**: User preferences and configuration

### 🔒 Security

- All data files are encrypted server-side with per-user keys
- Master key derivation using HKDF-SHA256
- AES-GCM encryption with unique nonces
- Encrypted storage path: `sandbox/{project_id}/{user_id}/{uuid}.enc`

### 🚀 Architecture

- **State Management**: React Context API
  - `AuthContext`: User authentication state
  - `ProjectsContext`: Project management with localStorage persistence
- **API Client**: Centralized axios configuration with JWT interceptors
- **Routing**: React Router v6 with protected routes
- **Styling**: Inline styles with consistent design system

## Getting Started

```bash
# Install dependencies
npm install

# Start development server
npm run dev -- --host

# Run lint checks
npm run lint

# Build production bundle
npm run build
```

## Environment Variables

Configure via `.env` or `.env.local` (see `infrastructure/docker-compose.yml` for defaults):

- `VITE_API_BASE_URL` — defaults to `http://localhost:8001`
- Additional secrets (auth, analytics) can be added as needed

## API Integration

The frontend talks to the FastAPI backend (see `../backend`). REST endpoints:

- `/auth/register` - User registration
- `/auth/login` - User login
- `/projects/` - Projects CRUD
- `/runs/` - Pipeline runs CRUD
- `/data/upload` - Encrypted file upload
- `/data/` - List files
- `/data/{id}/download-decrypted` - Download decrypted file

All API calls include JWT authentication via `Authorization: Bearer <token>` header.
