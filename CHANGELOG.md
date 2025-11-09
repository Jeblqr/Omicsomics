# Changelog - Omicsomics Platform

## [Unreleased]

### Added - 2024-01 (Major Feature Release: Projects, Runs, and Data with Encryption)

#### Backend
- **Encryption System**: Per-user data encryption using HKDF + AES-GCM
  - `backend/app/core/crypto.py`: Key derivation and encryption utilities
  - MASTER_KEY-based key derivation with HKDF-SHA256
  - AES-GCM encryption with 128-bit nonces
  - Functions: `derive_user_key`, `encrypt_for_user`, `decrypt_for_user`

- **Database Models**:
  - `Run` model for tracking pipeline executions
  - `DataFile` model for encrypted file metadata
  - Updated `User` model to SQLAlchemy 2.x Mapped notation
  - Migration: `0002_add_runs_and_datafiles.py`

- **Storage Service** (`backend/app/services/storage_service.py`):
  - S3/MinIO integration with aioboto3
  - Encrypted file upload with automatic key derivation
  - Decrypted file download with streaming support
  - Storage path convention: `sandbox/{project_id}/{user_id}/{uuid}.enc`

- **API Endpoints**:
  - `/data/upload` - POST multipart file upload with encryption
  - `/data/` - GET list files by project
  - `/data/{id}/download` - GET presigned URL for encrypted file
  - `/data/{id}/download-decrypted` - GET streaming decrypted file
  - `/runs/` - POST create run, GET list runs by project
  - `/runs/{id}` - GET run details

- **Tests**:
  - `test_crypto.py`: 5 tests for encryption/decryption (✅ all passing)
  - `test_runs.py`: 7 tests for run CRUD operations (✅ all passing)
  - `test_data.py`: 5 tests for data operations (requires MinIO integration)

#### Frontend
- **Projects Manager** (`frontend/src/pages/projects/ProjectsPage.tsx`):
  - Full CRUD interface (Create, Read, Update, Delete)
  - Inline editing with form toggle
  - Delete confirmation dialog
  - Project selection for data/runs operations

- **Runs Manager** (`frontend/src/pages/runs/RunsPage.tsx`):
  - Create new runs with name and description
  - List runs filtered by current project
  - Color-coded status badges (pending/running/completed/failed)
  - Timestamp display for created_at, started_at, finished_at

- **Data Manager** (`frontend/src/pages/data/DataPage.tsx`):
  - File upload with automatic encryption
  - File listing with metadata (size, checksum, timestamps)
  - Download decrypted files
  - Format file sizes (Bytes/KB/MB/GB)

- **State Management**:
  - `ProjectsContext`: Global project state with Context API
  - localStorage persistence for current project
  - Auto-select newly created projects
  - CRUD operations: create, update, delete, fetch

- **Components**:
  - `ProjectSwitcher`: Dropdown for project selection with visual feedback
  - Updated `Sidebar`: Links to all new pages
  - `Layout`: Consistent page wrapper

- **Developer Experience**:
  - TypeScript strict type checking
  - Proper error handling with user-friendly messages
  - Loading states during async operations
  - Disabled states when no project selected

#### Infrastructure
- **Scripts**:
  - `scripts/start_frontend.sh`: Frontend development server startup
  - Updated `scripts/start_backend.sh`: Backend startup with checks

- **Documentation**:
  - `IMPLEMENTATION_COMPLETE.md`: Complete feature documentation
  - `QUICKSTART.md`: Step-by-step setup guide
  - Updated `frontend/README.md`: Frontend features and architecture

### Changed
- User model refactored from `Column` to `Mapped[type]` notation for SQLAlchemy 2.x compatibility
- DataFile uses `metadata_` attribute to avoid SQLAlchemy declarative base conflict
- Frontend API client centralized with axios configuration
- ProjectsPage redesigned from simple list to full CRUD interface

### Fixed
- TypeScript type errors in ProjectsPage (data→projects, removed unused imports)
- SQLAlchemy metadata attribute conflict in DataFile model
- Test fixture naming consistency (async_db, async_client)
- React prop-types and no-explicit-any linting issues

### Security
- All user data files encrypted at rest with per-user keys
- Master key derivation with HKDF ensures key isolation
- Project ownership verification on all data/run operations
- JWT authentication required for all protected endpoints
- SHA256 checksums for data integrity verification

## [0.1.0] - Initial Release

### Added
- FastAPI backend with async SQLAlchemy
- React + TypeScript frontend with Vite
- User authentication with JWT
- PostgreSQL database
- Docker Compose infrastructure
- Basic project structure and scaffolding
- Health check endpoints
- CORS middleware
- Alembic migrations setup

---

## Notes

### Breaking Changes
None yet - this is pre-1.0 development

### Deprecations
None

### Migration Guide
If upgrading from earlier version:
1. Run `alembic upgrade head` to apply new migrations
2. Set `MASTER_KEY` environment variable (64 hex characters)
3. Configure MinIO/S3 storage settings
4. Restart backend and frontend servers

### Known Issues
- Data upload tests require MinIO to be running
- Some TypeScript linting warnings remain (non-blocking)
- No pagination for large file/run lists yet
- No file preview functionality yet

### Contributors
- Development: AI-assisted implementation
- Testing: Manual and automated testing
- Documentation: Comprehensive guides and API docs
