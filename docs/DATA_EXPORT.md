# Data Export Documentation

## Overview

The **Data Export** feature provides comprehensive batch export capabilities with format conversion, metadata inclusion, and progress tracking. Users can export multiple files at once, convert formats on-the-fly, and download packaged results.

## Features

### Core Capabilities

1. **Batch Export**

   - Export multiple files in a single operation
   - Select files from Data Browser or datasets
   - Queue-based processing for large exports

2. **Format Conversion**

   - **CSV**: Comma-separated values
   - **TSV**: Tab-separated values
   - **JSON**: JavaScript Object Notation
   - **Excel**: Microsoft Excel (.xlsx)
   - **Parquet**: Apache Parquet columnar format
   - **HDF5**: Hierarchical Data Format
   - **ZIP**: Archive multiple files

3. **Metadata & Lineage**

   - Include export metadata JSON
   - Track data lineage and provenance
   - File-level metadata preservation
   - Processing history

4. **Progress Tracking**

   - Real-time progress updates
   - Status indicators (pending, processing, completed, failed, cancelled)
   - Detailed error reporting
   - File-level success/failure tracking

5. **Auto-Management**
   - TTL-based auto-cleanup
   - Configurable expiration (1 hour to 1 week)
   - Automatic storage management
   - Download URL generation

## Architecture

### Database Models

#### ExportJob Model

```python
class ExportJob(Base):
    id: int
    job_key: str              # Unique identifier
    name: str
    description: str
    user_id: int
    project_id: int

    # Configuration
    export_format: ExportFormat
    include_metadata: bool
    include_lineage: bool
    compress: bool

    # File selection
    file_ids: List[int]       # Files to export
    dataset_ids: List[int]    # Optional: entire datasets
    export_options: dict      # Format-specific options

    # Output
    output_path: str          # Storage path
    output_size: int          # Size in bytes
    download_url: str         # Pre-signed URL

    # Status tracking
    status: ExportStatus      # pending/processing/completed/failed/cancelled
    progress: int             # 0-100%
    error_message: str
    total_files: int
    processed_files: int
    failed_files: int

    # Timestamps
    created_at: datetime
    started_at: datetime
    completed_at: datetime
    expires_at: datetime
```

#### ExportStatus Enum

```python
class ExportStatus(str, enum.Enum):
    PENDING = "pending"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
```

#### ExportFormat Enum

```python
class ExportFormat(str, enum.Enum):
    CSV = "csv"
    TSV = "tsv"
    JSON = "json"
    EXCEL = "excel"
    PARQUET = "parquet"
    HDF5 = "hdf5"
    ZIP = "zip"
```

### Service Layer

**DataExportService** (`backend/app/services/data_export.py`)

Main operations:

- `create_export_job()` - Create and queue export job
- `get_export_job()` - Get job by ID
- `get_export_job_by_key()` - Get job by unique key
- `list_export_jobs()` - List with filters
- `update_job_status()` - Update status and progress
- `update_job_progress()` - Update progress counters
- `cancel_export_job()` - Cancel pending/processing job
- `delete_export_job()` - Delete job and files
- `cleanup_expired_jobs()` - Auto-maintenance
- `get_export_statistics()` - Export statistics

### Async Processing

**Celery Task**: `process_export_job`

Processing pipeline:

1. Update status to PROCESSING
2. Create temporary directory
3. Download files from storage
4. Apply format conversions
5. Generate metadata manifest
6. Compress into ZIP (if requested)
7. Upload to storage
8. Generate download URL
9. Update status to COMPLETED

## API Reference

### Base URL

`/api/data-export`

### Endpoints

#### Create Export Job

```http
POST /jobs
Content-Type: application/json

{
  "name": "RNA-seq Results Export",
  "description": "Export processed RNA-seq data",
  "file_ids": [123, 124, 125],
  "export_format": "csv",
  "include_metadata": true,
  "include_lineage": false,
  "compress": true,
  "ttl_hours": 48
}

Response 200:
{
  "id": 1,
  "job_key": "export_1_20251110120000",
  "name": "RNA-seq Results Export",
  "status": "pending",
  "progress": 0,
  "total_files": 3,
  "created_at": "2025-11-10T12:00:00Z",
  ...
}
```

#### Get Export Job

```http
GET /jobs/{job_id}

Response 200:
{
  "id": 1,
  "job_key": "export_1_20251110120000",
  "name": "RNA-seq Results Export",
  "status": "completed",
  "progress": 100,
  "total_files": 3,
  "processed_files": 3,
  "failed_files": 0,
  "output_size": 1048576,
  "download_url": "https://storage.example.com/exports/...",
  "expires_at": "2025-11-12T12:00:00Z",
  ...
}
```

#### List Export Jobs

```http
GET /jobs?status=completed&limit=20

Response 200:
[
  {
    "id": 1,
    "job_key": "export_1_20251110120000",
    "name": "RNA-seq Results Export",
    "status": "completed",
    "progress": 100,
    ...
  },
  ...
]
```

#### Cancel Export Job

```http
POST /jobs/{job_id}/cancel

Response 200:
{
  "message": "Export job cancelled",
  "job_id": 1
}
```

#### Delete Export Job

```http
DELETE /jobs/{job_id}

Response 204: No Content
```

#### Get Statistics

```http
GET /statistics

Response 200:
{
  "total_jobs": 50,
  "by_status": {
    "pending": 2,
    "processing": 3,
    "completed": 40,
    "failed": 3,
    "cancelled": 2
  },
  "by_format": {
    "csv": 25,
    "json": 15,
    "excel": 10
  },
  "total_size_bytes": 104857600,
  "total_files_exported": 150
}
```

#### Cleanup Expired Jobs

```http
POST /cleanup

Response 200:
{
  "message": "Cleaned up 5 expired export jobs",
  "count": 5
}
```

## Usage Scenarios

### Scenario 1: Export Analysis Results

**User Story**: Researcher wants to export processed RNA-seq results for publication.

**Steps**:

1. **Select Files in Data Browser**

   - Browse to results directory
   - Select differential expression files
   - Note file IDs

2. **Create Export Job**

```javascript
const exportJob = await axios.post("/api/data-export/jobs", {
  name: "RNA-seq Publication Data",
  description: "Differential expression results for paper",
  file_ids: [123, 124, 125],
  export_format: "excel",
  include_metadata: true,
  include_lineage: true,
  compress: true,
  ttl_hours: 72, // 3 days
});
```

3. **Monitor Progress**

   - UI automatically refreshes every 5 seconds
   - Progress bar shows completion percentage
   - Status updates: pending → processing → completed

4. **Download Results**
   - Click "Download" button
   - Opens pre-signed URL in new tab
   - ZIP file contains:
     - `file_123_degs.xlsx`
     - `file_124_counts.xlsx`
     - `file_125_metadata.xlsx`
     - `export_metadata.json`

### Scenario 2: Batch Export Project Data

**User Story**: Export all results from a project for archival.

**Steps**:

1. **Query Project Files**

```javascript
const projectFiles = await axios.get("/api/files", {
  params: { project_id: 5 },
});
const fileIds = projectFiles.data.map((f) => f.id);
```

2. **Create Large Export**

```javascript
const exportJob = await axios.post("/api/data-export/jobs", {
  name: "Project 5 - Complete Archive",
  description: "All project files for archival",
  file_ids: fileIds, // May be 100+ files
  export_format: "zip",
  include_metadata: true,
  include_lineage: true,
  compress: true,
  ttl_hours: 168, // 1 week
});
```

3. **Background Processing**

   - Celery worker processes asynchronously
   - User can navigate away
   - Email notification on completion (optional)

4. **Download When Ready**
   - Return to Data Export page
   - See completed export
   - Download large ZIP archive

### Scenario 3: Format Conversion Export

**User Story**: Convert CSV files to Parquet for efficient storage.

**Steps**:

1. **Select CSV Files**

   - Identify CSV files for conversion
   - Note file IDs

2. **Export with Conversion**

```javascript
const exportJob = await axios.post("/api/data-export/jobs", {
  name: "CSV to Parquet Conversion",
  description: "Convert large CSV datasets to Parquet",
  file_ids: [200, 201, 202],
  export_format: "parquet",
  export_options: {
    compression: "snappy",
    row_group_size: 100000,
  },
  compress: false, // Don't ZIP parquet files
  ttl_hours: 24,
});
```

3. **Process and Download**
   - System converts CSV → Parquet
   - Preserves data types and metadata
   - Downloads optimized Parquet files

### Scenario 4: Cancel Long-Running Export

**User Story**: Cancel an export that's taking too long or was created by mistake.

**Steps**:

1. **Identify Problem Export**

   - See export stuck at 10% for 30 minutes
   - Or realize wrong files were selected

2. **Cancel Export**

```javascript
await axios.post(`/api/data-export/jobs/${jobId}/cancel`);
```

3. **Verify Cancellation**
   - Status changes to "cancelled"
   - Processing stops immediately
   - Can delete and recreate with correct settings

## Frontend Integration

### Component: DataExportPage

**Location**: `frontend/src/pages/DataExportPage.tsx`

**Key Features**:

- Export job list with cards
- Create export dialog
- Real-time progress tracking
- Download management
- Auto-refresh for active jobs

**State Management**:

```typescript
interface State {
  jobs: ExportJob[]; // All export jobs
  loading: boolean; // Loading state
  error: string | null; // Error message
  showCreateDialog: boolean; // Dialog visibility
  formData: ExportJobForm; // Create form data
}
```

**Key Operations**:

```typescript
// Create export
const handleCreateJob = async () => {
  await axios.post("/api/data-export/jobs", formData);
  loadJobs();
};

// Download export
const handleDownload = (job: ExportJob) => {
  window.open(job.download_url, "_blank");
};

// Cancel export
const handleCancelJob = async (jobId: number) => {
  await axios.post(`/api/data-export/jobs/${jobId}/cancel`);
  loadJobs();
};

// Delete export
const handleDeleteJob = async (jobId: number) => {
  await axios.delete(`/api/data-export/jobs/${jobId}`);
  loadJobs();
};
```

**Auto-Refresh**:

```typescript
useEffect(() => {
  const interval = setInterval(() => {
    const hasActiveJobs = jobs.some(
      (job) => job.status === "pending" || job.status === "processing"
    );
    if (hasActiveJobs) {
      loadJobs();
    }
  }, 5000); // Refresh every 5 seconds

  return () => clearInterval(interval);
}, [jobs]);
```

## Configuration

### Export Format Options

#### CSV Options

```json
{
  "export_format": "csv",
  "export_options": {
    "delimiter": ",",
    "quote_char": "\"",
    "encoding": "utf-8",
    "include_header": true
  }
}
```

#### Excel Options

```json
{
  "export_format": "excel",
  "export_options": {
    "sheet_name": "Data",
    "include_index": false,
    "freeze_panes": [1, 0]
  }
}
```

#### Parquet Options

```json
{
  "export_format": "parquet",
  "export_options": {
    "compression": "snappy",
    "row_group_size": 100000,
    "use_dictionary": true
  }
}
```

### TTL Configuration

```python
# Minimum: 1 hour
# Maximum: 168 hours (1 week)
# Default: 48 hours (2 days)

ttl_hours = 48
expires_at = datetime.utcnow() + timedelta(hours=ttl_hours)
```

### Metadata Manifest Structure

```json
{
  "export_info": {
    "job_key": "export_1_20251110120000",
    "name": "RNA-seq Results Export",
    "created_at": "2025-11-10T12:00:00Z",
    "format": "csv"
  },
  "files": [
    {
      "id": 123,
      "filename": "degs.csv",
      "size": 102400,
      "uploaded_at": "2025-11-09T10:00:00Z",
      "metadata": {
        "experiment": "RNA-seq",
        "sample_count": 12
      }
    },
    ...
  ],
  "lineage": {
    "source_datasets": [5, 6],
    "pipeline_runs": [10, 11],
    "processing_steps": [...]
  }
}
```

## Best Practices

### 1. File Selection

**DO**:

- Group related files in single export
- Use meaningful export names
- Include descriptions for context
- Verify file IDs before exporting

**DON'T**:

- Export too many files at once (>1000)
- Mix unrelated data types
- Use generic names like "Export 1"
- Export without checking file sizes

### 2. Format Selection

**DO**:

- Use CSV for maximum compatibility
- Use Parquet for large datasets
- Use Excel for end-user reports
- Use JSON for structured data

**DON'T**:

- Export large files to Excel (>100MB)
- Use HDF5 for small files
- Mix formats unnecessarily

### 3. TTL Management

**DO**:

- Set appropriate TTL (24-72 hours typical)
- Download exports promptly
- Delete after downloading
- Use longer TTL (1 week) for archival

**DON'T**:

- Use very short TTL (<4 hours)
- Let exports accumulate
- Rely on exports as backup storage

### 4. Progress Monitoring

**DO**:

- Monitor large exports
- Check error messages
- Cancel if stuck
- Retry failed exports

**DON'T**:

- Create duplicate exports
- Ignore failed files count
- Let failed exports accumulate

## Troubleshooting

### Issue: Export Stuck at 0%

**Symptoms**:

- Status is "processing"
- Progress stays at 0%
- No error message

**Solutions**:

```bash
# Check Celery worker status
celery -A app.celery_app inspect active

# Check Celery logs
docker-compose logs celery

# Restart Celery worker
docker-compose restart celery

# Cancel and recreate export
curl -X POST http://localhost:8000/api/data-export/jobs/{job_id}/cancel
```

### Issue: Download URL Expired

**Symptoms**:

- "URL expired" error when downloading
- 403 Forbidden from storage

**Solutions**:

```javascript
// Regenerate download URL (API enhancement needed)
const job = await axios.get(`/api/data-export/jobs/${jobId}`);
// Download URL is regenerated on each GET

// Or download before expiry
// Check expires_at field
if (new Date(job.expires_at) < new Date()) {
  alert("Export has expired. Please create a new export.");
}
```

### Issue: Export Failed with Error

**Symptoms**:

- Status is "failed"
- Error message displayed
- Some files may have processed

**Solutions**:

```javascript
// Check error message
const job = await axios.get(`/api/data-export/jobs/${jobId}`);
console.log("Error:", job.error_message);
console.log("Failed files:", job.failed_files);

// Common issues:
// 1. File not found - Check file_ids are valid
// 2. Format conversion error - Check file format compatibility
// 3. Storage error - Check storage configuration

// Retry with correct settings
await axios.post("/api/data-export/jobs", {
  ...correctedData,
});
```

### Issue: Large Export Takes Too Long

**Symptoms**:

- Export processing for hours
- Progress moving very slowly

**Solutions**:

```javascript
// Solution 1: Split into smaller batches
const fileIdsBatch1 = fileIds.slice(0, 100);
const fileIdsBatch2 = fileIds.slice(100, 200);

await axios.post('/api/data-export/jobs', {
  name: 'Export Batch 1',
  file_ids: fileIdsBatch1,
  ...
});

// Solution 2: Disable compression for faster processing
await axios.post('/api/data-export/jobs', {
  ...data,
  compress: false
});

// Solution 3: Use simpler format
await axios.post('/api/data-export/jobs', {
  ...data,
  export_format: 'zip'  // No conversion
});
```

## Maintenance

### Automated Cleanup

Schedule periodic cleanup:

```python
# Using APScheduler or cron
from apscheduler.schedulers.background import BackgroundScheduler

scheduler = BackgroundScheduler()

def cleanup_exports():
    from app.database import SessionLocal
    from app.services.data_export import DataExportService

    db = SessionLocal()
    try:
        service = DataExportService(db)
        count = service.cleanup_expired_jobs()
        print(f"Cleaned up {count} expired exports")
    finally:
        db.close()

# Run daily at 2 AM
scheduler.add_job(cleanup_exports, 'cron', hour=2)
scheduler.start()
```

### Storage Monitoring

Monitor export storage usage:

```python
def get_export_storage_usage():
    from app.models.data_export import ExportJob

    jobs = db.query(ExportJob).filter(
        ExportJob.status == ExportStatus.COMPLETED
    ).all()

    total_size = sum(job.output_size or 0 for job in jobs)
    return {
        'total_exports': len(jobs),
        'total_size_gb': total_size / (1024**3),
        'avg_size_mb': (total_size / len(jobs)) / (1024**2) if jobs else 0
    }
```

## Future Enhancements

1. **Advanced Features**

   - Dataset-level export (export all files in dataset)
   - Selective column export (choose columns from CSV)
   - Custom transformation scripts
   - Split large files during export

2. **Format Support**

   - NetCDF for scientific data
   - Feather for fast dataframes
   - Avro for streaming data
   - ORC for Hadoop ecosystems

3. **Notifications**

   - Email on completion
   - Slack/Teams integration
   - Browser push notifications
   - Webhook callbacks

4. **Sharing**

   - Generate shareable links
   - Set access permissions
   - Track download analytics
   - Temporary public URLs

5. **Performance**
   - Parallel file processing
   - Streaming large files
   - Incremental compression
   - CDN integration for downloads

## References

- [Pandas Export Formats](https://pandas.pydata.org/docs/user_guide/io.html)
- [Apache Parquet](https://parquet.apache.org/)
- [HDF5 Format](https://www.hdfgroup.org/)
- [Celery Best Practices](https://docs.celeryproject.org/en/stable/userguide/tasks.html)

---

**Version**: 1.0  
**Last Updated**: 2025-11-10  
**Author**: Omicsomics Development Team
