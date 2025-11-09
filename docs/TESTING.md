# Testing Guide

This document describes the testing strategy and CI/CD pipeline for Omicsomics.

## Test Structure

```
backend/tests/
├── conftest.py              # Shared fixtures
├── test_auth.py             # Authentication tests
├── test_projects.py         # Project management tests
├── test_samples.py          # Sample management tests
├── test_files.py            # File operations tests
└── integration/
    └── test_end_to_end.py   # End-to-end integration tests
```

## Running Tests Locally

### Backend Tests

```bash
cd backend

# Run all tests
pytest

# Run with coverage
pytest --cov=app --cov-report=html --cov-report=term-missing

# Run specific test file
pytest tests/test_auth.py -v

# Run specific test
pytest tests/test_auth.py::test_register_user -v

# Run with more output
pytest -vv -s
```

### Frontend Tests

```bash
cd frontend

# Lint check
npm run lint

# Type check
npx tsc --noEmit

# Build test
npm run build
```

## CI/CD Pipeline

The project uses GitHub Actions for continuous integration and deployment.

### Pipeline Overview

```
┌─────────────────────────────────────────────────────────┐
│                    Pull Request / Push                   │
└─────────────────────┬───────────────────────────────────┘
                      │
        ┌─────────────┴─────────────┐
        │                           │
┌───────▼────────┐         ┌────────▼───────┐
│  Backend Tests │         │ Frontend Tests │
├────────────────┤         ├────────────────┤
│ • Linting      │         │ • ESLint       │
│ • Type Check   │         │ • TypeScript   │
│ • Unit Tests   │         │ • Build        │
│ • Coverage     │         └────────┬───────┘
└───────┬────────┘                  │
        │                           │
        └─────────────┬─────────────┘
                      │
              ┌───────▼────────┐
              │ Integration    │
              │ Tests          │
              └───────┬────────┘
                      │
              ┌───────▼────────┐
              │ Docker Build   │
              │ Test           │
              └───────┬────────┘
                      │
                ✅ All Checks
                   Passed
```

### Workflow Jobs

#### 1. Backend Linting (`backend-lint`)
- **Runs**: Ruff linter and formatter
- **Purpose**: Code quality and style consistency
- **Commands**:
  ```bash
  ruff check app/ tests/
  ruff format --check app/ tests/
  ```

#### 2. Backend Tests (`backend-tests`)
- **Services**: PostgreSQL, Redis, MinIO
- **Runs**: Full test suite with coverage
- **Environment**:
  - PostgreSQL 16 (test database)
  - Redis 7 (task queue)
  - MinIO (object storage)
- **Coverage**: Uploads to Codecov

#### 3. Frontend Linting (`frontend-lint`)
- **Runs**: ESLint + TypeScript type checking
- **Commands**:
  ```bash
  npm run lint
  npx tsc --noEmit
  ```

#### 4. Frontend Build (`frontend-build`)
- **Runs**: Production build
- **Artifacts**: Uploads `dist/` folder
- **Purpose**: Ensures frontend builds successfully

#### 5. Docker Build (`docker-build`)
- **Trigger**: Pull requests only
- **Builds**: Backend and frontend Docker images
- **Purpose**: Validates Dockerfiles
- **Cache**: Uses GitHub Actions cache

#### 6. Integration Tests (`integration-tests`)
- **Depends on**: backend-tests, frontend-build
- **Services**: Full stack (db, redis, minio)
- **Runs**: End-to-end integration tests
- **Purpose**: Validates full system behavior

### Environment Variables

CI pipeline uses the following environment variables:

```yaml
# Database
DATABASE_URL: postgresql+asyncpg://omicsomics:testpassword123@localhost:5432/omicsomics_test

# Celery
CELERY_BROKER_URL: redis://localhost:6379/0
CELERY_RESULT_BACKEND: redis://localhost:6379/0

# Object Storage
OBJECT_STORAGE_ENDPOINT: http://localhost:9000
OBJECT_STORAGE_ACCESS_KEY: minioadmin
OBJECT_STORAGE_SECRET_KEY: minioadmin123
OBJECT_STORAGE_BUCKET: omicsomics-test

# Security (testing only)
SECRET_KEY: test_secret_key_for_ci_testing_only_12345678901234567890
ENCRYPTION_KEY: test_encryption_key_for_ci_32bytes_long_string_here
```

**⚠️ Production Note**: These are test credentials only. Production uses secure secrets from GitHub Secrets.

## Test Coverage

### Current Coverage

View coverage reports:
- **Badge**: See README.md
- **Codecov**: https://codecov.io/gh/Jeblqr/Omicsomics
- **Local**: Run `pytest --cov=app --cov-report=html` and open `htmlcov/index.html`

### Coverage Goals

| Component | Current | Target |
|-----------|---------|--------|
| Core API  | 85%+    | 90%+   |
| Services  | 70%+    | 85%+   |
| Models    | 90%+    | 95%+   |
| Overall   | 75%+    | 85%+   |

## Writing Tests

### Example: API Test

```python
import pytest
from httpx import AsyncClient

@pytest.mark.asyncio
async def test_create_project(async_client: AsyncClient, auth_headers: dict):
    """Test project creation."""
    response = await async_client.post(
        "/api/v1/projects/",
        headers=auth_headers,
        json={
            "name": "My Project",
            "description": "Test project"
        }
    )
    assert response.status_code == 201
    data = response.json()
    assert data["name"] == "My Project"
```

### Example: Service Test

```python
import pytest
from app.services.projects import ProjectService

@pytest.mark.asyncio
async def test_get_project_by_id(db_session, test_project):
    """Test getting project by ID."""
    service = ProjectService(db_session)
    project = await service.get_project(test_project.id)
    assert project is not None
    assert project.name == "Test Project"
```

### Example: Integration Test

```python
import pytest
from httpx import AsyncClient

@pytest.mark.asyncio
async def test_upload_process_workflow(
    async_client: AsyncClient,
    auth_headers: dict,
    test_project
):
    """Test full upload and processing workflow."""
    # Upload file
    files = {"file": ("test.csv", b"gene,count\nGene1,100", "text/csv")}
    response = await async_client.post(
        "/api/v1/data/upload",
        headers=auth_headers,
        files=files,
        data={"project_id": test_project.id}
    )
    assert response.status_code == 201
    
    # Check processed data
    data_id = response.json()["id"]
    response = await async_client.get(
        f"/api/v1/data/{data_id}/processed",
        headers=auth_headers
    )
    assert response.status_code == 200
```

## Pre-commit Hooks

Pre-commit hooks run automatically before each commit to ensure code quality.

### Setup

```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install

# Run manually on all files
pre-commit run --all-files
```

### Hooks Included

1. **trailing-whitespace**: Removes trailing whitespace
2. **end-of-file-fixer**: Ensures files end with newline
3. **check-yaml/json/toml**: Validates config files
4. **ruff**: Python linting and formatting
5. **mypy**: Type checking
6. **eslint**: JavaScript/TypeScript linting

## Troubleshooting

### Tests Fail Locally but Pass in CI

1. **Check environment variables**: Ensure `.env` is configured correctly
2. **Database state**: Run `alembic downgrade base && alembic upgrade head`
3. **Clean test database**: Remove `tests/test_database.sqlite3`
4. **Dependencies**: Run `pip install -e ".[dev]"` to update

### CI Pipeline Fails

1. **Check logs**: Click on failed job in GitHub Actions
2. **Service health**: Ensure PostgreSQL/Redis/MinIO containers are healthy
3. **Database migrations**: Check if Alembic migrations are up to date
4. **Secrets**: Verify GitHub Secrets are configured (for production)

### Coverage Drops

1. **Add tests**: Write tests for new features
2. **Remove dead code**: Delete unused code
3. **Integration tests**: Add end-to-end tests for full workflows

## Best Practices

### Do's ✅

- ✅ Write tests for new features before merging
- ✅ Maintain >80% test coverage for critical paths
- ✅ Use fixtures for common test data
- ✅ Mock external services (S3, external APIs)
- ✅ Test error cases and edge conditions
- ✅ Run tests locally before pushing

### Don'ts ❌

- ❌ Commit without running tests
- ❌ Skip flaky tests (fix them instead)
- ❌ Use production credentials in tests
- ❌ Hard-code test data (use fixtures)
- ❌ Ignore linter warnings
- ❌ Push code that breaks CI

## Resources

- **pytest documentation**: https://docs.pytest.org/
- **pytest-asyncio**: https://pytest-asyncio.readthedocs.io/
- **GitHub Actions**: https://docs.github.com/actions
- **Ruff**: https://docs.astral.sh/ruff/
- **Codecov**: https://docs.codecov.com/

## Contact

For CI/CD issues or questions:
- Open an issue: https://github.com/Jeblqr/Omicsomics/issues
- Check workflow runs: https://github.com/Jeblqr/Omicsomics/actions
