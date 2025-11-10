# Pipeline Dataset Integration Documentation

## Overview

Pipeline Dataset Integration seamlessly connects the **Dataset Manager** with the **Pipeline Builder**, enabling:

- **Use datasets as pipeline inputs**: Select organized datasets instead of individual files
- **Auto-create output datasets**: Automatically package pipeline results into versioned datasets
- **Complete lineage tracking**: Record full provenance from input datasets through pipeline operations to outputs
- **Bidirectional linking**: Track which pipelines used which datasets and which datasets were created by which runs

This integration closes the loop between data processing and data organization, making workflows more reproducible and data more discoverable.

## Architecture

### Components

1. **PipelineDatasetService** (`backend/app/services/pipeline_dataset_service.py`)

   - Core business logic for dataset-pipeline integration
   - 15+ methods for dataset operations, validation, and lineage

2. **Pipeline Dataset API** (`backend/app/api/pipeline_datasets.py`)

   - REST endpoints for frontend integration
   - 7 main endpoints for all integration scenarios

3. **DatasetInputSelector** (`frontend/src/components/DatasetInputSelector.tsx`)

   - UI component for selecting datasets as pipeline inputs
   - Shows dataset cards with files, validates compatibility

4. **CreateDatasetFromRunDialog** (`frontend/src/components/CreateDatasetFromRunDialog.tsx`)
   - Dialog for creating datasets from completed runs
   - Auto-fill mode and manual configuration

### Data Flow

```
Input Phase:
  Dataset → DatasetInputSelector → Pipeline Node Input → Run Execution

Output Phase:
  Run Completion → Auto-create Dataset → Dataset Manager → Lineage Records
```

## API Reference

### Base URL

All endpoints are prefixed with `/api/pipeline-datasets`

### Endpoints

#### 1. Get Dataset Input Files

Get files from a dataset suitable for pipeline input.

```http
POST /input-files
```

**Request Body:**

```json
{
  "dataset_id": 123,
  "file_role": "primary" // optional: filter by role
}
```

**Response:**

```json
{
  "dataset_id": 123,
  "files": [
    {
      "id": 456,
      "path": "/data/processed/sample1.csv",
      "name": "sample1.csv",
      "type": "csv",
      "role": "primary",
      "size": 1048576
    }
  ]
}
```

#### 2. Create Dataset from Run

Create a dataset from pipeline run outputs.

```http
POST /create-from-run
```

**Request Body:**

```json
{
  "run_id": 789,
  "dataset_name": "RNA-seq Analysis Results",
  "dataset_description": "Processed RNA-seq data with DEG analysis",
  "data_type": "transcriptomics",
  "file_roles": {
    "101": "primary",
    "102": "metadata"
  },
  "tags": ["pipeline-output", "deg-analysis", "qc-passed"]
}
```

**Response:**

```json
{
  "dataset": {
    "id": 200,
    "name": "RNA-seq Analysis Results",
    "description": "Processed RNA-seq data with DEG analysis",
    "data_type": "transcriptomics",
    "status": "active",
    "file_count": 5,
    "created_at": "2025-11-10T10:30:00Z"
  },
  "message": "Dataset created successfully with 5 files"
}
```

#### 3. Link Datasets to Run

Link input datasets to a pipeline run for tracking.

```http
POST /link-to-run
```

**Request Body:**

```json
{
  "run_id": 789,
  "dataset_ids": [123, 124],
  "input_mapping": {
    "node_1_input_data": 123,
    "node_3_reference": 124
  }
}
```

**Response:**

```json
{
  "message": "Linked 2 datasets to run 789",
  "run_id": 789,
  "dataset_ids": [123, 124]
}
```

#### 4. Get Run Output Datasets

Get all datasets created from a pipeline run.

```http
GET /run/{run_id}/output-datasets
```

**Response:**

```json
[
  {
    "id": 200,
    "name": "RNA-seq Analysis Results",
    "data_type": "transcriptomics",
    "status": "active",
    "file_count": 5,
    "created_at": "2025-11-10T10:30:00Z"
  }
]
```

#### 5. Get Dataset Usage in Runs

Get all pipeline runs that used a dataset as input.

```http
GET /dataset/{dataset_id}/usage
```

**Response:**

```json
{
  "dataset_id": 123,
  "runs": [
    {
      "run_id": 789,
      "run_name": "RNA-seq Pipeline v2",
      "status": "completed",
      "started_at": "2025-11-10T10:00:00Z",
      "finished_at": "2025-11-10T10:30:00Z"
    }
  ]
}
```

#### 6. Validate Dataset for Pipeline

Validate dataset compatibility with pipeline requirements.

```http
POST /validate
```

**Request Body:**

```json
{
  "dataset_id": 123,
  "required_file_types": ["fastq", "fasta"] // optional
}
```

**Response:**

```json
{
  "is_valid": true,
  "errors": []
}
```

Or with errors:

```json
{
  "is_valid": false,
  "errors": [
    "Missing required file types: fastq",
    "File not found: reference.fasta"
  ]
}
```

#### 7. Auto-Create Dataset from Run

Automatically create dataset with inferred settings.

```http
POST /auto-create
```

**Request Body:**

```json
{
  "run_id": 789,
  "auto_tags": true
}
```

**Response:**

```json
{
  "dataset": {
    "id": 201,
    "name": "RNA-seq Pipeline v2_output",
    "data_type": "transcriptomics",
    "status": "active",
    "file_count": 5,
    "created_at": "2025-11-10T10:35:00Z"
  },
  "message": "Dataset auto-created successfully with 5 files"
}
```

## Usage Scenarios

### Scenario 1: Use Existing Dataset in New Pipeline

**User Story**: Researcher wants to run quality control on a previously uploaded RNA-seq dataset.

**Steps**:

1. **Open Pipeline Builder** and create a new QC pipeline

2. **Add QC tool node** (e.g., FastQC)

3. **Click on input port** → Select "Use Dataset"

4. **DatasetInputSelector opens**:

   - Filter by data type: "transcriptomics"
   - Select the RNA-seq dataset
   - Review and select specific files
   - Validate compatibility

5. **Run pipeline** with dataset files as inputs

6. **System automatically**:
   - Links dataset to run in metadata
   - Records lineage for traceability

**API Calls**:

```javascript
// Get dataset files
const response = await axios.post("/api/pipeline-datasets/input-files", {
  dataset_id: datasetId,
});

// Link to run
await axios.post("/api/pipeline-datasets/link-to-run", {
  run_id: runId,
  dataset_ids: [datasetId],
});
```

### Scenario 2: Auto-Create Dataset After Pipeline Completion

**User Story**: After running a variant calling pipeline, automatically save outputs as a versioned dataset.

**Steps**:

1. **Pipeline completes successfully**

2. **RunDetailsPage shows** "Create Dataset" button

3. **User clicks button** → CreateDatasetFromRunDialog opens

4. **Two options**:

   **Option A - Auto Create**:

   - Click "Auto Create" button
   - System infers data type from pipeline
   - Auto-adds tags: "pipeline-output", "pipeline-template"
   - Creates dataset instantly

   **Option B - Manual Configure**:

   - Edit dataset name and description
   - Select data type
   - Add custom tags
   - Click "Create Dataset"

5. **Result**:
   - Dataset created with all output files
   - MD5/SHA256 hashes computed
   - Lineage records created linking run → dataset
   - Run metadata updated with dataset reference

**API Calls**:

```javascript
// Auto-create
const response = await axios.post("/api/pipeline-datasets/auto-create", {
  run_id: runId,
  auto_tags: true,
});

// Or manual create
const response = await axios.post("/api/pipeline-datasets/create-from-run", {
  run_id: runId,
  dataset_name: "Variant Calling Results",
  data_type: "genomics",
  tags: ["vcf", "filtered", "hg38"],
});
```

### Scenario 3: Track Complete Data Lineage

**User Story**: Trace the full history of a final analysis result back through multiple pipeline stages.

**Steps**:

1. **Open Dataset Manager** → Select final result dataset

2. **Go to "Lineage" tab** in detail view

3. **View lineage tree**:

   ```
   Input Dataset A (Raw FASTQ)
       ↓ [QC Pipeline]
   Dataset B (QC-passed FASTQ)
       ↓ [Alignment Pipeline]
   Dataset C (BAM files)
       ↓ [Variant Calling Pipeline]
   Dataset D (VCF files) ← Current dataset
   ```

4. **Each lineage record shows**:

   - Operation type: "pipeline"
   - Pipeline name and version
   - Tool chain used
   - Execution timestamps
   - Parameters used

5. **Click on any intermediate dataset** to:
   - View its files
   - See which runs used it as input
   - Navigate to those run details

**API Calls**:

```javascript
// Get dataset lineage
const lineage = await axios.get(`/api/datasets/${datasetId}/lineage`);

// Get full file lineage tree
const fileLineage = await axios.get(`/api/datasets/files/${fileId}/lineage`);

// Get runs that used this dataset
const usage = await axios.get(
  `/api/pipeline-datasets/dataset/${datasetId}/usage`
);
```

### Scenario 4: Validate Dataset Before Pipeline Execution

**User Story**: Before starting an expensive pipeline run, verify that the input dataset meets all requirements.

**Steps**:

1. **User selects dataset** in Pipeline Builder

2. **System automatically validates**:

   - Dataset status is "active"
   - All files exist on disk
   - Required file types are present
   - File integrity (hashes match)

3. **Validation results displayed**:

   ```
   ✓ Dataset is active
   ✓ All 10 files found
   ✓ Required types present: fastq, fasta
   ✗ Warning: 2 files missing hash verification
   ```

4. **If validation fails**:

   - Pipeline execution blocked
   - User prompted to fix issues
   - Suggestions provided

5. **If validation passes**:
   - User can proceed with confidence
   - Run has better chance of success

**API Calls**:

```javascript
// Validate dataset
const validation = await axios.post("/api/pipeline-datasets/validate", {
  dataset_id: datasetId,
  required_file_types: ["fastq", "fasta"],
});

if (!validation.data.is_valid) {
  alert("Dataset validation failed: " + validation.data.errors.join(", "));
}
```

## Service Layer Methods

### PipelineDatasetService Class

Located in `backend/app/services/pipeline_dataset_service.py`

#### Input Operations

**`get_dataset_files_for_input(dataset_id, file_role=None)`**

- Get files from dataset suitable for pipeline input
- Optionally filter by file role
- Returns list of file info dicts

**`prepare_dataset_files_for_pipeline(dataset_id, node_id, input_name)`**

- Prepare dataset files for specific pipeline node input
- Returns list of file paths ready for execution

**`validate_dataset_for_pipeline_input(dataset_id, required_file_types=None)`**

- Validate dataset meets pipeline requirements
- Returns (is_valid, error_list) tuple

#### Output Operations

**`create_dataset_from_run(run, dataset_name, ...)`**

- Create dataset from pipeline run outputs
- Add all output files with hash computation
- Apply tags and metadata
- Returns created Dataset object

**`auto_create_output_dataset(run, auto_tags=True)`**

- Automatically create dataset with inferred settings
- Smart data type detection
- Auto-tagging based on pipeline type
- Returns Dataset or None

#### Lineage Operations

**`record_pipeline_lineage(dataset_id, run, input_file_ids, output_file_ids)`**

- Record complete lineage for pipeline operations
- Link input files → operation → output files
- Store tool names, versions, and parameters
- Returns list of lineage records

**`link_input_datasets_to_run(run_id, dataset_ids, input_mapping=None)`**

- Link input datasets to run for tracking
- Store dataset references in run metadata
- Optional input mapping for multi-dataset runs
- Returns updated Run object

#### Query Operations

**`get_run_output_datasets(run_id)`**

- Get all datasets created from a run
- Returns list of Dataset objects

**`get_dataset_usage_in_runs(dataset_id)`**

- Get all runs that used dataset as input
- Returns list of run info dicts

**`_infer_data_type_from_run(run)`**

- Internal method to infer data type from run
- Checks pipeline template, parameters, file types
- Returns data type string

## Frontend Integration

### DatasetInputSelector Component

**Props**:

```typescript
interface DatasetInputSelectorProps {
  open: boolean; // Dialog open state
  onClose: () => void; // Close handler
  onSelect: (
    // Selection handler
    datasetId: number,
    files: DatasetFile[]
  ) => void;
  nodeId: string; // Pipeline node ID
  inputName: string; // Input parameter name
  requiredFileTypes?: string[]; // Required file types
  projectId?: number; // Filter by project
}
```

**Features**:

- Dataset browsing with filters (data type, search)
- File list with role and type chips
- Multi-file selection with checkboxes
- Real-time validation
- Visual feedback for compatibility

**Usage Example**:

```tsx
<DatasetInputSelector
  open={showDatasetSelector}
  onClose={() => setShowDatasetSelector(false)}
  onSelect={(datasetId, files) => {
    // Update node input with dataset files
    updateNodeInput(nodeId, {
      dataset_id: datasetId,
      file_paths: files.map((f) => f.path),
    });
  }}
  nodeId="alignment_node_1"
  inputName="input_fastq"
  requiredFileTypes={["fastq", "fastq.gz"]}
  projectId={currentProject.id}
/>
```

### CreateDatasetFromRunDialog Component

**Props**:

```typescript
interface CreateDatasetFromRunDialogProps {
  open: boolean; // Dialog open state
  onClose: () => void; // Close handler
  run: Run | null; // Run object
  onSuccess?: (
    // Success callback
    datasetId: number
  ) => void;
}
```

**Features**:

- Auto-fill from run information
- Data type inference
- Auto-tagging option
- Manual configuration mode
- Real-time validation

**Usage Example**:

```tsx
<CreateDatasetFromRunDialog
  open={showCreateDataset}
  onClose={() => setShowCreateDataset(false)}
  run={completedRun}
  onSuccess={(datasetId) => {
    // Navigate to dataset or show success message
    navigate(`/datasets/${datasetId}`);
  }}
/>
```

## Integration Points

### 1. Pipeline Builder Integration

**File**: `frontend/src/pages/PipelineBuilderPage.tsx`

**Changes**:

- Add "Use Dataset" button in node input configuration
- Show DatasetInputSelector when clicked
- Store dataset reference in node data
- Display dataset name in node UI

**Code Example**:

```tsx
// In node configuration panel
{
  inputType === "file" && (
    <Box display="flex" gap={1}>
      <Button onClick={() => setShowFileSelector(true)}>Select Files</Button>
      <Button
        onClick={() => setShowDatasetSelector(true)}
        startIcon={<DatasetOutlined />}
      >
        Use Dataset
      </Button>
    </Box>
  );
}
```

### 2. Run Details Integration

**File**: `frontend/src/pages/RunDetailsPage.tsx`

**Changes**:

- Add "Create Dataset" button for completed runs
- Show CreateDatasetFromRunDialog
- Display linked input datasets
- Show created output datasets

**Code Example**:

```tsx
// In run actions
{
  run.status === "completed" && run.output_files?.length > 0 && (
    <Button
      variant="outlined"
      onClick={() => setShowCreateDataset(true)}
      startIcon={<SaveOutlined />}
    >
      Save as Dataset
    </Button>
  );
}

<CreateDatasetFromRunDialog
  open={showCreateDataset}
  onClose={() => setShowCreateDataset(false)}
  run={run}
  onSuccess={(datasetId) => {
    showSnackbar("Dataset created successfully!");
    loadOutputDatasets();
  }}
/>;
```

### 3. Dataset Manager Integration

**File**: `frontend/src/pages/DatasetManagerPage.tsx`

**Changes**:

- Add "Usage" tab showing runs that used this dataset
- Add "Created By" field showing source run
- Link to run details from dataset view

**Code Example**:

```tsx
// In dataset detail dialog
<Tab label="Usage" />;

{
  activeTab === "usage" && (
    <Box>
      <Typography variant="subtitle2" gutterBottom>
        Used in Pipeline Runs
      </Typography>
      {usageRuns.map((run) => (
        <Card key={run.run_id}>
          <CardContent>
            <Typography variant="h6">{run.run_name}</Typography>
            <Chip label={run.status} />
            <Button onClick={() => navigate(`/runs/${run.run_id}`)}>
              View Run
            </Button>
          </CardContent>
        </Card>
      ))}
    </Box>
  );
}
```

### 4. Data Browser Integration

**File**: `frontend/src/pages/DataBrowserPage.tsx`

**Changes**:

- Add "Create Dataset" action for selected files
- Quick-create dataset from file selection
- Show dataset membership for files

**Code Example**:

```tsx
// In file actions
<Button
  onClick={() => createDatasetFromFiles(selectedFiles)}
  startIcon={<DatasetOutlined />}
>
  Create Dataset from Selection
</Button>
```

## Best Practices

### 1. Dataset Organization

**DO**:

- Create datasets for all significant pipeline outputs
- Use consistent naming: `{project}_{stage}_{date}`
- Add descriptive tags for easy filtering
- Include run metadata in dataset description

**DON'T**:

- Create datasets for temporary intermediate files
- Use generic names like "output1", "data"
- Skip data type classification
- Forget to add tags

### 2. Lineage Tracking

**DO**:

- Always link input datasets to runs
- Record lineage for all dataset-creating operations
- Include tool versions and parameters
- Maintain bidirectional links (dataset↔run)

**DON'T**:

- Skip lineage for "simple" operations
- Forget to update run metadata
- Lose track of intermediate datasets
- Break lineage chains

### 3. Validation

**DO**:

- Validate datasets before pipeline execution
- Check file existence and integrity
- Verify required file types present
- Provide clear error messages

**DON'T**:

- Skip validation to save time
- Ignore validation warnings
- Allow pipelines to run with invalid inputs
- Provide vague error messages

### 4. Auto-Creation

**DO**:

- Use auto-create for standard workflows
- Review and configure for important results
- Add meaningful descriptions and tags
- Set correct data types

**DON'T**:

- Over-rely on auto-create for all scenarios
- Accept default names without review
- Skip manual tagging for special cases
- Create datasets from failed runs

## Troubleshooting

### Issue: Dataset Files Not Appearing in Pipeline Builder

**Symptoms**:

- DatasetInputSelector shows dataset but no files
- "No files found" message

**Possible Causes**:

1. Dataset has no files added
2. Files filtered by role
3. File paths are broken

**Solutions**:

```bash
# Check dataset files
curl -X GET "/api/datasets/{dataset_id}/files"

# Validate dataset
curl -X POST "/api/pipeline-datasets/validate" \
  -H "Content-Type: application/json" \
  -d '{"dataset_id": 123}'

# Re-add files if needed
curl -X POST "/api/datasets/{dataset_id}/files" \
  -H "Content-Type: application/json" \
  -d '{"file_path": "/path/to/file.csv", "compute_hash": true}'
```

### Issue: Auto-Create Dataset Fails

**Symptoms**:

- "Failed to create dataset" error
- Empty dataset created

**Possible Causes**:

1. Run is not completed
2. Run has no output files
3. Output files don't exist on disk
4. Permission issues

**Solutions**:

```python
# Check run status
run = db.query(Run).filter(Run.id == run_id).first()
print(f"Status: {run.status}")
print(f"Output files: {run.output_files}")

# Verify files exist
for file_id in run.output_files:
    df = db.query(DataFile).filter(DataFile.id == file_id).first()
    print(f"{df.filename}: exists={Path(df.file_path).exists()}")

# Manual create with error handling
try:
    dataset = service.create_dataset_from_run(run)
except Exception as e:
    print(f"Error: {e}")
```

### Issue: Lineage Not Recorded

**Symptoms**:

- Lineage tab empty in Dataset Manager
- Missing provenance information

**Possible Causes**:

1. Lineage recording skipped
2. Input/output file IDs not provided
3. Dataset created outside pipeline system

**Solutions**:

```python
# Manually record lineage
from app.services.pipeline_dataset_service import PipelineDatasetService

service = PipelineDatasetService(db)
service.record_pipeline_lineage(
    dataset_id=dataset_id,
    run=run,
    input_file_ids=[...],  # DatasetFileEntry IDs
    output_file_ids=[...]  # DatasetFileEntry IDs
)
```

### Issue: Validation Fails for Valid Dataset

**Symptoms**:

- "Dataset not valid" despite files being present
- False positive errors

**Possible Causes**:

1. File paths changed after dataset creation
2. File type mismatch (e.g., "fastq.gz" vs "fastq")
3. Symbolic link issues

**Solutions**:

```python
# Check actual file paths
dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
for file in dataset.files:
    path = Path(file.file_path)
    print(f"{file.file_name}: exists={path.exists()}, resolved={path.resolve()}")

# Update file paths if needed
file_entry.file_path = str(new_path)
db.commit()

# Re-compute hashes if needed
service = DatasetManager(db)
service.add_file(dataset_id, new_path, compute_hash=True)
```

## Future Enhancements

1. **Smart Dataset Recommendations**

   - Recommend compatible datasets based on pipeline type
   - Suggest similar datasets used in past runs
   - Auto-detect format compatibility

2. **Dataset Versioning in Pipelines**

   - Support dataset version pinning
   - Auto-upgrade to newer dataset versions
   - Compare results across dataset versions

3. **Batch Pipeline Execution**

   - Run same pipeline on multiple datasets
   - Parallel execution with dataset queue
   - Aggregate results into single output dataset

4. **Dataset Provenance Visualization**

   - Interactive lineage graph
   - Timeline view of dataset history
   - Visual diff between dataset versions

5. **Dataset Quality Metrics**

   - Track dataset usage frequency
   - Monitor pipeline success rates by dataset
   - Flag problematic datasets

6. **Cross-Project Dataset Sharing**

   - Share datasets between projects
   - Dataset permissions and access control
   - Dataset marketplace/library

7. **Smart Caching**

   - Cache dataset file listings
   - Prefetch commonly used datasets
   - Optimize validation queries

8. **Integration with External Sources**
   - Import datasets from URLs
   - Sync with cloud storage (S3, GCS)
   - Connect to data repositories (GEO, SRA)

## References

- [Dataset Manager Documentation](./DATASET_MANAGER.md)
- [Pipeline Builder Documentation](./PIPELINE_BUILDER.md)
- [Data Browser Documentation](./DATA_FLOW_AND_REQUIREMENTS.md)
- [Format Conversion System](./FORMAT_CONVERSION_SYSTEM.md)

---

**Version**: 1.0  
**Last Updated**: 2025-11-10  
**Author**: Omicsomics Development Team
