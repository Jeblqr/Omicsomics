# Archive Processing Feature

## Overview

The archive processing feature allows users to upload compressed archives containing multiple omics data files, preview the contents, and selectively extract and process individual files without having to manually decompress the archive.

## Supported Formats

- **ZIP** (`.zip`)
- **TAR** (`.tar`)
- **TAR.GZ** (`.tar.gz`, `.tgz`)
- **TAR.BZ2** (`.tar.bz2`, `.tbz`)

## Features

### 1. Secure Archive Extraction

- **Path Traversal Protection**: Prevents malicious archives from accessing files outside the extraction directory
- **Sanitization**: Removes dangerous path components (`..,`, absolute paths, drive letters)
- **Safe Extraction**: Validates all file paths before extraction

### 2. Archive Preview

List all files in an archive without extracting:

```bash
GET /api/v1/data/{datafile_id}/preview
```

**Response:**

```json
{
  "is_archive": true,
  "archive_filename": "omics_data.zip",
  "total_files": 15,
  "total_size": 52428800,
  "files": [
    {
      "name": "proteins.csv",
      "path": "proteomics/proteins.csv",
      "size": 1024000,
      "mime_type": "text/csv",
      "is_dir": false
    },
    ...
  ]
}
```

### 3. Selective File Processing

Extract and process a specific file from an archive:

```bash
POST /api/v1/data/{datafile_id}/process-from-archive
Content-Type: multipart/form-data

file_path=proteomics/proteins.csv
sample_id=sample_001
```

**Response:**

```json
{
  "datafile_id": 123,
  "processed_file_id": 124,
  "extracted_filename": "proteins.csv",
  "archive_source": {
    "archive_id": 100,
    "archive_filename": "omics_data.zip",
    "file_path": "proteomics/proteins.csv"
  },
  "processing_info": {
    "success": true,
    "format_detected": "csv",
    "data_type": "proteomics",
    "features_count": 1500
  }
}
```

## Usage Examples

### Web UI

1. **Upload Archive**

   - Navigate to project data page
   - Upload a ZIP or TAR.GZ file
   - Disable "Process file" option (to preview first)

2. **Preview Contents**

   - Open the uploaded archive file
   - Click "Preview Archive" button
   - Browse files in the archive

3. **Process Files**
   - Select files from the preview list
   - Click "Process" for individual files
   - View processed unified format data

### API Usage

```python
import requests

# Upload archive
files = {'file': open('omics_data.zip', 'rb')}
data = {'project_id': 1, 'process_file': 'false'}
headers = {'Authorization': f'Bearer {token}'}

response = requests.post(
    'http://api.example.com/api/v1/data/upload',
    files=files,
    data=data,
    headers=headers
)
archive_id = response.json()['id']

# Preview contents
preview = requests.get(
    f'http://api.example.com/api/v1/data/{archive_id}/preview',
    headers=headers
).json()

print(f"Archive contains {preview['total_files']} files")
for file in preview['files']:
    print(f"  - {file['path']} ({file['size']} bytes)")

# Process specific file
process_data = {'file_path': 'proteomics/proteins.csv'}
result = requests.post(
    f'http://api.example.com/api/v1/data/{archive_id}/process-from-archive',
    data=process_data,
    headers=headers
).json()

print(f"Processed file ID: {result['datafile_id']}")
print(f"Format: {result['processing_info']['format_detected']}")
```

## Security Considerations

### Path Traversal Prevention

The archive utilities implement multiple layers of security:

1. **Name Sanitization**: Removes dangerous path components

   - Parent directory references (`../`)
   - Absolute paths (`/etc/passwd`, `C:\Windows\...`)
   - Drive letters (Windows)

2. **Path Validation**: Verifies extraction targets

   - Resolves paths to absolute form
   - Checks paths remain within extraction directory
   - Rejects paths that escape the base directory

3. **File Filtering**: Skips unsafe files during listing and extraction

### Example Attack Prevention

```python
# Malicious archive attempt:
# File path: "../../../etc/passwd"

# Sanitization result:
# UnsafeArchiveError: Archive contains unsafe path: ../../../etc/passwd

# The file is skipped and never extracted
```

## Architecture

### Backend Components

```
backend/
├── app/
│   ├── utils/
│   │   └── archive.py          # Core archive utilities
│   ├── api/
│   │   └── routers/
│   │       └── data.py          # Archive API endpoints
│   └── services/
│       └── file_processor.py    # File processing logic
└── tests/
    ├── test_archive.py          # Unit tests
    └── integration/
        └── test_archive_processing.py  # Integration tests
```

### Frontend Components

```
frontend/
└── src/
    └── components/
        └── ArchivePreview.tsx   # Archive preview UI
```

## Performance Considerations

- **Streaming**: Archives are processed in memory using temporary files
- **Cleanup**: Temporary files are automatically deleted after use
- **Size Limits**: Follow standard file upload limits (default: 100MB sync, unlimited async)
- **Large Archives**: For archives > 100MB, consider async processing

## Future Enhancements

1. **Batch Processing**: Process multiple selected files simultaneously
2. **Archive Streaming**: Stream large archives without full download
3. **Nested Archives**: Support archives within archives
4. **Preview Content**: Show file previews without full extraction
5. **Smart Filtering**: Filter files by type or pattern
6. **Archive Creation**: Create archives from processed files

## Testing

### Unit Tests

```bash
pytest tests/test_archive.py -v
```

Tests:

- Path safety validation
- Archive format detection
- File listing
- Extraction safety
- Error handling

### Integration Tests

```bash
pytest tests/integration/test_archive_processing.py -v
```

Tests:

- Complete workflow (upload → preview → extract → process)
- Multiple archive formats
- Multiple data types
- Error scenarios

## Error Handling

| Error             | HTTP Code | Description                             |
| ----------------- | --------- | --------------------------------------- |
| Not an archive    | 200       | File is not a supported archive format  |
| Invalid archive   | 400       | Archive is corrupted or invalid         |
| Unsafe path       | 400       | Archive contains path traversal attempt |
| File not found    | 404       | Requested file not in archive           |
| Processing failed | 500       | Failed to process extracted file        |

## Dependencies

- Python standard library:
  - `zipfile` - ZIP archive handling
  - `tarfile` - TAR/TAR.GZ/TAR.BZ2 handling
  - `pathlib` - Path manipulation
  - `mimetypes` - MIME type detection
  - `tempfile` - Temporary file handling

No external dependencies required for core functionality.

## License

This feature is part of the Omicsomics platform and follows the same license.
