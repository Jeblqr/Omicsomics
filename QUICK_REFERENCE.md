# Quick Reference Guide

## 🚀 Getting Started

### Start Development Environment

```bash
# Start all services
./manage.sh start

# View logs
./manage.sh logs

# Stop all services
./manage.sh stop

# Restart services
./manage.sh restart
```

### Access Points

- **Frontend**: http://localhost:5173
- **Backend API**: http://localhost:8000
- **API Docs**: http://localhost:8000/docs
- **MinIO Console**: http://localhost:9001

---

## 📁 Project Structure

```
Omicsomics/
├── backend/                 # Backend API
│   ├── app/
│   │   ├── models/         # Database models (12+ models)
│   │   ├── services/       # Business logic (6+ services)
│   │   ├── api/            # API endpoints (60+ endpoints)
│   │   └── core/           # Configuration & utilities
│   ├── alembic/            # Database migrations
│   └── tests/              # Backend tests
│
├── frontend/                # Frontend React app
│   ├── src/
│   │   ├── pages/          # Main pages (5 new pages)
│   │   ├── components/     # Reusable components
│   │   └── config.ts       # Configuration
│   └── public/             # Static assets
│
├── docs/                    # Documentation (3,100+ lines)
│   ├── DATA_EXPORT.md
│   ├── DATA_EDITING.md
│   ├── CUSTOM_SCRIPT_TOOLS.md
│   └── IMPLEMENTATION_SUMMARY.md
│
├── scripts/                 # Utility scripts & tests
│   ├── test_data_export.py
│   ├── test_custom_scripts.py
│   └── test_integration.py
│
└── infrastructure/          # Docker configuration
    └── docker-compose.yml
```

---

## 🔧 Common Commands

### Backend

```bash
cd backend

# Install dependencies
pip install -r requirements.txt

# Run migrations
alembic upgrade head

# Create migration
alembic revision --autogenerate -m "Description"

# Start server
uvicorn app.main:app --reload

# Start Celery worker
celery -A app.celery_worker worker --loglevel=info

# Run tests
pytest tests/
```

### Frontend

```bash
cd frontend

# Install dependencies
npm install

# Start dev server
npm run dev

# Build for production
npm run build

# Preview production build
npm run preview

# Type check
npm run type-check
```

### Database

```bash
# Connect to PostgreSQL
docker exec -it omicsomics-postgres psql -U postgres -d omicsomics

# Backup database
docker exec omicsomics-postgres pg_dump -U postgres omicsomics > backup.sql

# Restore database
docker exec -i omicsomics-postgres psql -U postgres omicsomics < backup.sql
```

### MinIO

```bash
# Access MinIO CLI
docker exec -it omicsomics-minio mc

# List buckets
docker exec omicsomics-minio mc ls local/

# Copy file to bucket
docker exec omicsomics-minio mc cp file.txt local/omicsomics/
```

---

## 🎯 API Quick Reference

### Datasets

```bash
# Create dataset
POST /api/datasets
{
  "name": "My Dataset",
  "description": "Dataset description",
  "data_type": "genomics",
  "tags": ["test"]
}

# List datasets
GET /api/datasets?data_type=genomics&search=test

# Get dataset
GET /api/datasets/{id}

# Update dataset
PUT /api/datasets/{id}

# Delete dataset
DELETE /api/datasets/{id}
```

### Data Export

```bash
# Create export job
POST /api/data-export/jobs
{
  "name": "My Export",
  "file_ids": [1, 2, 3],
  "export_format": "csv",
  "include_metadata": true
}

# Get job status
GET /api/data-export/jobs/{id}

# List jobs
GET /api/data-export/jobs

# Cancel job
POST /api/data-export/jobs/{id}/cancel

# Delete job
DELETE /api/data-export/jobs/{id}
```

### Data Editing

```bash
# Create edit session
POST /api/data-editing/sessions
{
  "name": "Clean Data",
  "file_id": 123
}

# Add operation
POST /api/data-editing/sessions/{id}/operations
{
  "operation_type": "filter_rows",
  "parameters": {
    "column": "age",
    "operator": "gt",
    "value": 18
  }
}

# Preview operations
POST /api/data-editing/sessions/{id}/preview

# Apply operations
POST /api/data-editing/sessions/{id}/apply

# Revert operations
POST /api/data-editing/sessions/{id}/revert
```

### Custom Scripts

```bash
# Create script
POST /api/custom-scripts/scripts
{
  "name": "My Script",
  "description": "Script description",
  "language": "python",
  "script_content": "...",
  "parameters_schema": {...}
}

# Execute script
POST /api/custom-scripts/scripts/{id}/execute
{
  "parameters": {"param": "value"},
  "input_file_ids": [1, 2]
}

# Get execution status
GET /api/custom-scripts/executions/{id}

# List executions
GET /api/custom-scripts/executions

# Get templates
GET /api/custom-scripts/templates
```

### Visualization Workspace

```bash
# Create workspace
POST /api/viz-workspace/workspaces
{
  "name": "My Dashboard",
  "layout": {"type": "grid", "rows": 2, "cols": 2}
}

# Add panel
POST /api/viz-workspace/workspaces/{id}/panels
{
  "panel_type": "scatter",
  "title": "Scatter Plot",
  "config": {"file_id": 123, "x_column": "x", "y_column": "y"},
  "position": {"row": 0, "col": 0}
}

# Get workspace
GET /api/viz-workspace/workspaces/{id}
```

---

## 🧪 Testing

### Run All Tests

```bash
# Backend unit tests
cd backend
pytest tests/

# Feature tests
cd scripts
python test_data_export.py --token YOUR_TOKEN
python test_custom_scripts.py --token YOUR_TOKEN
python test_integration.py --token YOUR_TOKEN

# Run all tests
./scripts/run_integration_tests.sh
```

### Test with Authentication

```bash
# Get auth token (replace with your credentials)
TOKEN=$(curl -X POST "http://localhost:8000/api/auth/login" \
  -H "Content-Type: application/json" \
  -d '{"username":"user","password":"pass"}' | jq -r .access_token)

# Run tests with token
python test_data_export.py --token $TOKEN
```

---

## 📊 Monitoring

### Check Service Health

```bash
# API health
curl http://localhost:8000/health

# Database connectivity
docker exec omicsomics-postgres pg_isready

# MinIO status
curl http://localhost:9000/minio/health/live

# Redis status
docker exec omicsomics-redis redis-cli ping
```

### View Logs

```bash
# All services
docker-compose logs -f

# Specific service
docker-compose logs -f backend
docker-compose logs -f frontend
docker-compose logs -f postgres
docker-compose logs -f minio
docker-compose logs -f redis

# Celery worker logs
docker-compose logs -f celery-worker
```

### Monitor Resources

```bash
# Docker stats
docker stats

# Database size
docker exec omicsomics-postgres psql -U postgres -d omicsomics \
  -c "SELECT pg_size_pretty(pg_database_size('omicsomics'));"

# MinIO storage usage
docker exec omicsomics-minio du -sh /data
```

---

## 🐛 Troubleshooting

### Backend Issues

```bash
# Check backend logs
docker-compose logs backend

# Restart backend
docker-compose restart backend

# Check database connection
docker exec omicsomics-postgres psql -U postgres -d omicsomics -c "SELECT 1;"

# Run migrations
docker-compose exec backend alembic upgrade head
```

### Frontend Issues

```bash
# Clear cache and rebuild
cd frontend
rm -rf node_modules package-lock.json
npm install
npm run dev

# Check for TypeScript errors
npm run type-check

# Check build
npm run build
```

### Database Issues

```bash
# Reset database (WARNING: deletes all data)
docker-compose down -v
docker-compose up -d postgres
docker-compose exec backend alembic upgrade head

# Check database connections
docker exec omicsomics-postgres psql -U postgres -c "SELECT * FROM pg_stat_activity;"
```

### Celery Issues

```bash
# Check Celery logs
docker-compose logs celery-worker

# Restart Celery
docker-compose restart celery-worker

# Check Redis
docker exec omicsomics-redis redis-cli ping

# Clear Redis queue
docker exec omicsomics-redis redis-cli FLUSHALL
```

---

## 🔍 Useful Queries

### Database Queries

```sql
-- Count datasets
SELECT COUNT(*) FROM datasets;

-- Recent export jobs
SELECT * FROM export_jobs ORDER BY created_at DESC LIMIT 10;

-- Script execution statistics
SELECT
  cs.name,
  COUNT(se.id) as total_executions,
  SUM(CASE WHEN se.status = 'completed' THEN 1 ELSE 0 END) as successful
FROM custom_scripts cs
LEFT JOIN script_executions se ON cs.id = se.script_id
GROUP BY cs.id, cs.name;

-- Active edit sessions
SELECT * FROM edit_sessions WHERE status IN ('draft', 'previewing', 'applying');
```

### API Queries

```bash
# Get all datasets
curl http://localhost:8000/api/datasets

# Get export job statistics
curl http://localhost:8000/api/data-export/statistics

# Get script templates
curl http://localhost:8000/api/custom-scripts/templates | jq

# Get workspace list
curl http://localhost:8000/api/viz-workspace/workspaces | jq
```

---

## 📚 Documentation Links

| Topic          | Link                                                             |
| -------------- | ---------------------------------------------------------------- |
| Data Export    | [docs/DATA_EXPORT.md](docs/DATA_EXPORT.md)                       |
| Data Editing   | [docs/DATA_EDITING.md](docs/DATA_EDITING.md)                     |
| Custom Scripts | [docs/CUSTOM_SCRIPT_TOOLS.md](docs/CUSTOM_SCRIPT_TOOLS.md)       |
| Implementation | [docs/IMPLEMENTATION_SUMMARY.md](docs/IMPLEMENTATION_SUMMARY.md) |
| Architecture   | [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)                     |
| Quick Start    | [QUICKSTART.md](QUICKSTART.md)                                   |

---

## 💡 Tips & Tricks

### Performance

- Use pagination for large datasets
- Enable database indexing for frequent queries
- Use Celery for long-running tasks
- Cache static data with Redis

### Development

- Use `--reload` flag for hot reloading
- Enable debug mode in development
- Use TypeScript strict mode
- Write tests before features

### Debugging

- Check logs first
- Use browser DevTools for frontend issues
- Use FastAPI's interactive docs (`/docs`)
- Enable verbose logging in development

### Security

- Never commit secrets to git
- Use environment variables for config
- Keep dependencies updated
- Enable HTTPS in production
- Validate all user inputs

---

## 🚀 Deployment Checklist

### Pre-Deployment

- [ ] All tests passing
- [ ] Documentation updated
- [ ] Environment variables configured
- [ ] Database migrations ready
- [ ] Backup strategy in place

### Production Setup

- [ ] SSL/TLS certificates installed
- [ ] CORS configured correctly
- [ ] Rate limiting enabled
- [ ] Monitoring set up
- [ ] Logging configured
- [ ] Backup jobs scheduled

### Post-Deployment

- [ ] Smoke tests passed
- [ ] Performance metrics baseline
- [ ] Error monitoring active
- [ ] User documentation available
- [ ] Support channels ready

---

_Last Updated: 2025-11-10_
