# Omics Platform Test Scripts

This directory contains comprehensive test scripts for the Omics analysis platform.

## Test Scripts

### 1. Quick Test (`quick_test.py`)

Fast API health check and basic functionality test.

```bash
python scripts/quick_test.py
```

**Tests:**

- Health endpoint
- User authentication (register/login)
- Project CRUD operations
- Pipeline template listing

**Duration:** ~5 seconds

### 2. Comprehensive Test (`comprehensive_omics_test.py`)

Full end-to-end testing of all omics modules with real data workflows.

```bash
# Set environment variables (optional)
export API_URL=http://localhost:8001/api/v1
export TEST_EMAIL=test@example.com
export TEST_PASSWORD=testpass123

# Run comprehensive tests
python scripts/comprehensive_omics_test.py
```

**Tests:**

- ✅ Genomics workflow (variant calling, FASTQ/VCF processing)
- ✅ Transcriptomics workflow (RNA-seq, differential expression)
- ✅ Proteomics workflow (protein identification, quantification)
- ✅ Metabolomics workflow (metabolite profiling)
- ✅ Epigenomics workflow (methylation analysis)
- ✅ Single-cell workflow (scRNA-seq clustering)
- ✅ Multi-omics integration (custom pipeline creation)

**Features:**

- Automatically generates test data
- Creates projects for each omics type
- Uploads data files
- Creates and executes analysis runs
- Generates detailed logs and JSON results

**Duration:** ~2-5 minutes

**Output:**

- Log file: `comprehensive_test.log`
- Results JSON: `test_results/test_results.json`
- Test data: `test_data/` directory

## Test Data

Test data is automatically generated for each omics type:

| Omics Type      | Generated Files                                | Description                |
| --------------- | ---------------------------------------------- | -------------------------- |
| Genomics        | `test_genomics.fastq`, `test_variants.vcf`     | 1000 reads, 100 variants   |
| Transcriptomics | `test_rna_counts.csv`, `test_rna_metadata.csv` | 1000 genes, 4 samples      |
| Proteomics      | `test_protein_abundance.csv`                   | 500 proteins, 4 samples    |
| Metabolomics    | `test_metabolite_abundance.csv`                | 300 metabolites, 4 samples |
| Epigenomics     | `test_methylation.csv`                         | 1000 CpG sites, 4 samples  |
| Single-cell     | `test_singlecell_matrix.csv`                   | 500 genes, 100 cells       |

## Environment Variables

```bash
# API endpoint (default: http://localhost:8001/api/v1)
export API_URL=http://localhost:8001/api/v1

# Test user credentials
export TEST_EMAIL=test@example.com
export TEST_PASSWORD=testpass123
```

## Docker Testing

To test against the Docker deployment:

```bash
# Quick test
docker exec infrastructure-backend-1 python /app/../scripts/quick_test.py

# Or from host with port forwarding
API_URL=http://localhost:8001/api/v1 python scripts/comprehensive_omics_test.py
```

## Continuous Integration

Add to CI/CD pipeline:

```yaml
# .github/workflows/test.yml
- name: Run Quick Tests
  run: python scripts/quick_test.py

- name: Run Comprehensive Tests
  run: python scripts/comprehensive_omics_test.py
  env:
    API_URL: http://localhost:8001/api/v1
```

## Interpreting Results

### Success

```
TEST SUMMARY
============================================================
genomics: ✅ PASSED
transcriptomics: ✅ PASSED
proteomics: ✅ PASSED
metabolomics: ✅ PASSED
epigenomics: ✅ PASSED
singlecell: ✅ PASSED
multiomics: ✅ PASSED

Total: 7/7 tests passed
Success rate: 100.0%
```

### Failure Investigation

1. Check `comprehensive_test.log` for detailed error messages
2. Verify backend is running: `docker ps`
3. Check backend logs: `docker logs infrastructure-backend-1`
4. Verify database migrations: `docker exec infrastructure-backend-1 alembic current`
5. Test individual endpoints with curl

## Troubleshooting

### Connection Refused

```bash
# Check backend is running
docker ps | grep backend

# Check backend logs
docker logs infrastructure-backend-1 --tail 50
```

### Authentication Errors

```bash
# Test login manually
curl -X POST http://localhost:8001/api/v1/auth/login \
  -d "username=test@example.com&password=testpass123"
```

### File Upload Errors

```bash
# Check MinIO is running
docker ps | grep minio

# Check storage configuration
docker logs infrastructure-minio-1
```

## Development

To add new tests:

1. Add test data generator in `DataDownloader` class
2. Implement workflow test function
3. Add to workflows list in `main()`
4. Update this README

Example:

```python
def download_new_omics_data(self) -> List[Path]:
    files = []
    # Generate test data
    return files

def test_new_omics_workflow(client, project_id, data_files) -> bool:
    # Test workflow
    return True
```

## Performance Benchmarks

Expected execution times on typical hardware:

| Test Suite    | Duration | API Calls | File Uploads |
| ------------- | -------- | --------- | ------------ |
| Quick Test    | ~5s      | 8         | 0            |
| Comprehensive | ~2-5min  | ~50       | 12           |

## Requirements

```bash
pip install requests
```

All other dependencies are Python standard library.
