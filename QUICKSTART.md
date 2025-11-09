# 🚀 Quick Start Guide - Omicsomics Platform

## Prerequisites

- Python 3.10+
- Node.js 18+
- PostgreSQL 14+
- MinIO or S3-compatible storage

## 1. Backend Setup

### Environment Configuration

```bash
cd backend
cp .env.example .env

# Edit .env with your configuration:
# MASTER_KEY=<generate-64-hex-chars>  # python -c "import secrets; print(secrets.token_hex(32))"
# DATABASE_URL=postgresql+asyncpg://user:pass@localhost:5432/omicsomics
# S3_ENDPOINT_URL=http://localhost:9000
# S3_ACCESS_KEY=minioadmin
# S3_SECRET_KEY=minioadmin
# S3_BUCKET_NAME=omicsomics-data
```

### Install Dependencies

```bash
pip install -e .
# or
poetry install
```

### Initialize Database

```bash
# Run migrations
alembic upgrade head
```

### Start Backend Server

```bash
uvicorn app.main:app --reload --host 0.0.0.0 --port 8001
# or use convenience script
../scripts/start_backend.sh
```

Backend will be available at: http://localhost:8001
API docs at: http://localhost:8001/docs

## 2. MinIO Setup

### Using Docker

```bash
docker run -d \
  -p 9000:9000 \
  -p 9001:9001 \
  --name minio \
  -e "MINIO_ROOT_USER=minioadmin" \
  -e "MINIO_ROOT_PASSWORD=minioadmin" \
  minio/minio server /data --console-address ":9001"
```

### Create Bucket

```bash
# Using mc (MinIO client)
mc alias set local http://localhost:9000 minioadmin minioadmin
mc mb local/omicsomics-data

# Or via scripts
python scripts/init_minio.py
```

MinIO Console: http://localhost:9001

## 3. Frontend Setup

### Install Dependencies

```bash
cd frontend
npm install
```

### Environment Configuration

```bash
# Create .env file (optional, defaults work for local dev)
echo "VITE_API_BASE_URL=http://localhost:8001" > .env
```

### Start Development Server

```bash
npm run dev -- --host
# or use convenience script
../scripts/start_frontend.sh
```

Frontend will be available at: http://localhost:5173

## 4. First Steps

### Create Account

1. Navigate to http://localhost:5173
2. Click "Register" if not logged in
3. Create account with username and password

### Create First Project

1. Click "Projects" in sidebar
2. Click "Create New Project"
3. Enter project name and description
4. Click "Create Project"
5. Click "Select" to activate the project

### Upload Data Files

1. With a project selected, click "Data" in sidebar
2. Click "Upload File"
3. Select a file from your computer
4. Click "Upload"
5. File will be encrypted and uploaded to MinIO

### Create Pipeline Run

1. With a project selected, click "Runs" in sidebar
2. Click "New Run"
3. Enter run name and description
4. Click "Create Run"
5. Run will appear in the list with "pending" status

## 5. Verify Encryption

### Check Encrypted Storage

```bash
# Connect to MinIO
mc alias set local http://localhost:9000 minioadmin minioadmin

# List encrypted files
mc ls local/omicsomics-data/sandbox/

# Files are stored as:
# sandbox/{project_id}/{user_id}/{uuid}.enc
```

### Download Decrypted File

1. In Data page, click "Download" button
2. Browser will download decrypted file with original name
3. Verify file content matches original upload

## 6. Testing

### Backend Tests

```bash
cd backend
pytest tests/ -v
# ✅ test_crypto.py: 5 passed
# ✅ test_runs.py: 7 passed
# ⚠️  test_data.py: requires MinIO running
```

### Frontend Linting

```bash
cd frontend
npm run lint
```

### Frontend Build

```bash
cd frontend
npm run build
# Should complete without errors
```

## 7. Docker Compose (All-in-One)

### Start All Services

```bash
cd infrastructure
docker-compose up -d
```

Services:

- PostgreSQL: localhost:5432
- MinIO: localhost:9000 (API), localhost:9001 (Console)
- Backend: localhost:8001
- Frontend: localhost:3000

### Stop All Services

```bash
docker-compose down
```

## 8. Common Tasks

### Generate MASTER_KEY

```bash
python -c "import secrets; print(secrets.token_hex(32))"
# Copy output to .env as MASTER_KEY
```

### Reset Database

```bash
# Downgrade all migrations
alembic downgrade base

# Re-apply migrations
alembic upgrade head
```

### Check Backend Logs

```bash
# If using script
tail -f logs/backend.log

# If using docker-compose
docker-compose logs -f backend
```

### View MinIO Data

```bash
# Using mc client
mc ls --recursive local/omicsomics-data/
```

## 9. API Usage Examples

### Register User

```bash
curl -X POST http://localhost:8001/auth/register \
  -H "Content-Type: application/json" \
  -d '{"username":"alice","password":"secret123"}'
```

### Login

```bash
TOKEN=$(curl -X POST http://localhost:8001/auth/login \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=alice&password=secret123" | jq -r '.access_token')
```

### Create Project

```bash
curl -X POST http://localhost:8001/projects/ \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"name":"My Project","description":"Test project"}'
```

### Upload File

```bash
curl -X POST http://localhost:8001/data/upload \
  -H "Authorization: Bearer $TOKEN" \
  -F "file=@data.csv" \
  -F "project_id=1"
```

### Create Run

```bash
curl -X POST http://localhost:8001/runs/ \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"name":"Analysis Run 1","description":"Test run","project_id":1}'
```

## 10. Troubleshooting

### Backend won't start

- Check DATABASE_URL is correct
- Ensure PostgreSQL is running: `psql -U postgres -l`
- Check MASTER_KEY is set (64 hex characters)

### Frontend can't connect to backend

- Check VITE_API_BASE_URL in frontend/.env
- Verify backend is running: `curl http://localhost:8001/health`
- Check browser console for CORS errors

### File upload fails

- Ensure MinIO is running: `curl http://localhost:9000/minio/health/live`
- Verify bucket exists: `mc ls local/omicsomics-data/`
- Check S3 credentials in backend/.env

### Migration fails

- Check database connection
- Verify alembic.ini has correct database URL
- Try: `alembic stamp head` then `alembic upgrade head`

### Decryption fails

- Verify MASTER_KEY hasn't changed
- Check user_id matches between encryption and decryption
- Review backend logs for detailed error messages

## 📚 Additional Resources

- Backend API docs: http://localhost:8001/docs
- Frontend README: frontend/README.md
- Architecture docs: docs/architecture/overview.md
- Complete implementation details: IMPLEMENTATION_COMPLETE.md

## 🎉 Success Indicators

✅ Backend starts without errors  
✅ Frontend builds successfully  
✅ Can create account and login  
✅ Can create projects  
✅ Can upload files (encrypted storage)  
✅ Can download files (decrypted)  
✅ Can create runs  
✅ All pages load without errors

Enjoy using Omicsomics! 🧬🔬
