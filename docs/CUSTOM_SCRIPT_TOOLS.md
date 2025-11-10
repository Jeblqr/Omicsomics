# Custom Script Tools

## Overview

The Custom Script Tools feature allows users to upload, manage, and execute custom Python, R, and Bash scripts within the Omicsomics platform. This enables researchers to integrate their own analysis tools and workflows seamlessly into the platform.

## Features

### Core Functionality

- **Script Management**: Upload, view, edit, and delete custom scripts
- **Multi-Language Support**: Python, R, and Bash scripts
- **Parameter Configuration**: JSON Schema-based parameter validation
- **Sandboxed Execution**: Secure script execution with resource limits
- **Result Tracking**: Complete execution history with outputs and logs
- **Visibility Control**: Private, project, or public script sharing
- **Script Templates**: Pre-built templates to get started quickly
- **Script Verification**: Mark trusted scripts as verified

### Security Features

- **Sandboxed Environment**: Scripts run in isolated environments
- **Resource Limits**: Configurable timeout and memory limits
- **Input Validation**: Parameter validation against JSON Schema
- **Output Isolation**: Separate output directories for each execution

### Collaboration Features

- **Visibility Levels**:
  - **Private**: Only visible to owner
  - **Project**: Visible to project members
  - **Public**: Visible to all users
- **Script Verification**: Platform admins can verify trusted scripts
- **Execution Statistics**: Track usage and success rates

## Architecture

### Database Models

#### CustomScript

Stores script metadata and configuration:

```python
- id: Primary key
- script_key: Unique identifier
- name: Script name
- description: Script description
- user_id: Owner user ID
- language: Script language (python, r, bash)
- script_content: The actual script code
- entry_point: Entry point function/method (optional)
- parameters_schema: JSON Schema for parameters
- requirements: List of dependencies
- timeout: Execution timeout (seconds)
- max_memory: Maximum memory (MB)
- visibility: Visibility level
- category: Script category
- tags: Script tags
- is_verified: Verification status
- total_executions: Total execution count
- successful_executions: Successful execution count
- failed_executions: Failed execution count
- avg_duration: Average execution duration
- last_executed_at: Last execution timestamp
- created_at: Creation timestamp
- updated_at: Update timestamp
```

#### ScriptExecution

Tracks script execution results:

```python
- id: Primary key
- execution_key: Unique identifier
- script_id: Foreign key to CustomScript
- user_id: Executing user ID
- parameters: Execution parameters (JSON)
- input_file_ids: List of input file IDs
- output_file_ids: List of output file IDs
- status: Execution status
- description: Execution description
- result_data: Result data (JSON)
- output_text: Standard output
- error_text: Standard error
- exit_code: Exit code
- duration: Execution duration (seconds)
- memory_usage: Memory usage (MB)
- started_at: Start timestamp
- completed_at: Completion timestamp
- created_at: Creation timestamp
```

### Service Layer

**CustomScriptService** provides:

- **Script Management**:

  - `create_script()`: Upload new script
  - `get_script()`: Get script with visibility check
  - `list_scripts()`: List scripts with filters
  - `update_script()`: Update script (owner only)
  - `delete_script()`: Delete script (owner only)

- **Execution Management**:

  - `execute_script()`: Execute script with parameters
  - `get_execution()`: Get execution details
  - `list_executions()`: List executions with filters
  - `cancel_execution()`: Cancel running execution

- **Validation**:

  - `_validate_json_schema()`: Validate JSON Schema
  - `_validate_parameters()`: Validate parameters against schema

- **Execution Engine**:
  - `_run_script_execution()`: Run script in sandboxed environment

### API Endpoints

All endpoints are prefixed with `/api/custom-scripts`.

#### Script Management

**POST /scripts**
Create a new script.

Request:

```json
{
  "name": "My Analysis Script",
  "description": "Performs statistical analysis",
  "language": "python",
  "script_content": "import pandas as pd\n...",
  "parameters_schema": {
    "type": "object",
    "properties": {
      "input_file": { "type": "string" },
      "threshold": { "type": "number", "default": 0.05 }
    },
    "required": ["input_file"]
  },
  "requirements": ["pandas", "numpy"],
  "timeout": 300,
  "max_memory": 1024,
  "visibility": "private",
  "category": "statistics",
  "tags": ["analysis", "statistics"]
}
```

Response: `201 Created` with script object

**GET /scripts/{script_id}**
Get script details (full content).

Response:

```json
{
  "id": 1,
  "script_key": "script_my_analysis_1699123456789",
  "name": "My Analysis Script",
  "description": "Performs statistical analysis",
  "language": "python",
  "script_content": "import pandas as pd\n...",
  "parameters_schema": {...},
  "requirements": ["pandas", "numpy"],
  "timeout": 300,
  "visibility": "private",
  "total_executions": 5,
  "successful_executions": 4,
  "failed_executions": 1,
  "created_at": "2024-11-10T10:00:00Z"
}
```

**GET /scripts**
List scripts with filters.

Query Parameters:

- `language`: Filter by language (python, r, bash)
- `category`: Filter by category
- `visibility`: Filter by visibility
- `tags`: Comma-separated tags
- `search`: Search in name and description
- `verified_only`: Show only verified scripts
- `skip`: Pagination offset
- `limit`: Page size (max 100)

Response:

```json
{
  "items": [...],
  "total": 10,
  "skip": 0,
  "limit": 50
}
```

**PUT /scripts/{script_id}**
Update script (owner only).

Request: Same as POST, all fields optional

Response: Updated script object

**DELETE /scripts/{script_id}**
Delete script (owner only).

Response: `204 No Content`

#### Script Execution

**POST /scripts/{script_id}/execute**
Execute a script.

Request:

```json
{
  "parameters": {
    "input_file": "data.csv",
    "threshold": 0.05
  },
  "input_file_ids": [1, 2, 3],
  "description": "Analysis run for experiment X"
}
```

Response: `202 Accepted` with execution object

**GET /executions/{execution_id}**
Get execution details.

Response:

```json
{
  "id": 1,
  "execution_key": "exec_my_analysis_1699123456789",
  "script_id": 1,
  "status": "completed",
  "parameters": {...},
  "output_text": "Analysis complete\nResults saved",
  "error_text": null,
  "exit_code": 0,
  "duration": 12.5,
  "output_file_ids": [10, 11],
  "started_at": "2024-11-10T10:00:00Z",
  "completed_at": "2024-11-10T10:00:12Z"
}
```

**GET /executions**
List executions.

Query Parameters:

- `script_id`: Filter by script
- `status`: Filter by status
- `skip`: Pagination offset
- `limit`: Page size

Response: Paginated execution list

**POST /executions/{execution_id}/cancel**
Cancel a running execution.

Response: `204 No Content`

#### Utilities

**POST /scripts/{script_id}/validate-parameters**
Validate parameters against schema.

Request:

```json
{
  "parameters": {
    "input_file": "data.csv",
    "threshold": 0.05
  }
}
```

Response:

```json
{
  "valid": true,
  "errors": null
}
```

**GET /languages**
Get supported languages.

Response: `["python", "r", "bash"]`

**GET /categories**
Get all categories in use.

Response: `["statistics", "preprocessing", "visualization"]`

**GET /templates**
Get script templates.

Response: Array of template objects with example scripts

### Frontend

**CustomScriptToolsPage** features:

- **Three Tabs**:

  1. **Scripts**: Browse and manage scripts
  2. **Executions**: View execution history
  3. **Templates**: Quick start with templates

- **Script Management**:

  - Create script dialog with form
  - Script cards with key information
  - View script details with full content
  - Execute scripts with parameter input
  - Delete scripts with confirmation

- **Execution Monitoring**:

  - Execution table with status indicators
  - Auto-refresh for active executions
  - View detailed execution logs
  - Cancel running executions

- **Filters**:
  - Search by name/description
  - Filter by language
  - Filter by category

## Usage Examples

### Example 1: Creating a Python Analysis Script

```python
# Step 1: Create script
POST /api/custom-scripts/scripts
{
  "name": "Gene Expression Analysis",
  "description": "Analyzes gene expression data using DESeq2",
  "language": "python",
  "script_content": """
import pandas as pd
import json
import sys
import os

def analyze_expression(counts_file, metadata_file, threshold):
    # Load data
    counts = pd.read_csv(counts_file, index_col=0)
    metadata = pd.read_csv(metadata_file)

    # Perform analysis (simplified)
    results = counts[counts.mean(axis=1) > threshold]

    # Save results
    os.makedirs('output', exist_ok=True)
    results.to_csv('output/filtered_genes.csv')

    return {
        'total_genes': len(counts),
        'filtered_genes': len(results),
        'output_file': 'filtered_genes.csv'
    }

if __name__ == '__main__':
    with open(sys.argv[2], 'r') as f:
        params = json.load(f)

    result = analyze_expression(
        params['counts_file'],
        params['metadata_file'],
        params['threshold']
    )
    print(json.dumps(result))
""",
  "parameters_schema": {
    "type": "object",
    "properties": {
      "counts_file": {
        "type": "string",
        "description": "Path to counts CSV file"
      },
      "metadata_file": {
        "type": "string",
        "description": "Path to metadata CSV file"
      },
      "threshold": {
        "type": "number",
        "default": 10,
        "description": "Expression threshold"
      }
    },
    "required": ["counts_file", "metadata_file"]
  },
  "requirements": ["pandas"],
  "timeout": 600,
  "visibility": "private",
  "category": "genomics",
  "tags": ["gene-expression", "analysis"]
}
```

### Example 2: Executing a Script

```python
# Step 1: Execute script
POST /api/custom-scripts/scripts/1/execute
{
  "parameters": {
    "counts_file": "counts.csv",
    "metadata_file": "metadata.csv",
    "threshold": 15
  },
  "input_file_ids": [100, 101],
  "description": "Analysis for experiment ABC"
}

# Response: Execution created with status "pending"

# Step 2: Check execution status
GET /api/custom-scripts/executions/1

# Response: Execution details with status, output, and results
```

### Example 3: Using R Script

```r
# Create R script for statistical analysis
POST /api/custom-scripts/scripts
{
  "name": "ANOVA Analysis",
  "description": "Performs ANOVA on experimental data",
  "language": "r",
  "script_content": """
library(jsonlite)

# Parse parameters
args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(args[2])

# Load data
data <- read.csv(params$input_file)

# Perform ANOVA
model <- aov(value ~ group, data=data)
summary_data <- summary(model)

# Save results
dir.create('output', showWarnings = FALSE)
write.csv(summary_data, 'output/anova_results.csv')

# Return results
result <- list(
  p_value = summary_data[[1]]$'Pr(>F)'[1],
  output_file = 'anova_results.csv'
)

cat(toJSON(result))
""",
  "parameters_schema": {
    "type": "object",
    "properties": {
      "input_file": {"type": "string"}
    },
    "required": ["input_file"]
  },
  "requirements": ["jsonlite"],
  "timeout": 300,
  "language": "r"
}
```

### Example 4: Bash Script for File Processing

```bash
# Create bash script for preprocessing
POST /api/custom-scripts/scripts
{
  "name": "FASTQ Quality Filter",
  "description": "Filters FASTQ files by quality score",
  "language": "bash",
  "script_content": """
#!/bin/bash

PARAMS_FILE=$2
INPUT_FILE=$(jq -r '.input_file' $PARAMS_FILE)
MIN_QUALITY=$(jq -r '.min_quality // 20' $PARAMS_FILE)

mkdir -p output

# Filter by quality (simplified example)
seqtk seq -Q $MIN_QUALITY $INPUT_FILE > output/filtered.fastq

echo "{\"output_file\": \"filtered.fastq\"}"
""",
  "parameters_schema": {
    "type": "object",
    "properties": {
      "input_file": {"type": "string"},
      "min_quality": {"type": "integer", "default": 20}
    },
    "required": ["input_file"]
  },
  "requirements": ["seqtk", "jq"],
  "timeout": 1800
}
```

## Script Development Guidelines

### Parameter Handling

Scripts receive parameters via a JSON file passed as command-line argument:

**Python**:

```python
import json
import sys

with open(sys.argv[2], 'r') as f:
    params = json.load(f)
```

**R**:

```r
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(args[2])
```

**Bash**:

```bash
PARAMS_FILE=$2
VALUE=$(jq -r '.param_name' $PARAMS_FILE)
```

### Output Files

Save output files to the `output/` directory:

```python
import os
os.makedirs('output', exist_ok=True)
result.to_csv('output/results.csv')
```

The system will automatically:

- Upload output files to MinIO
- Create File records in database
- Attach file IDs to execution record

### Return Values

Print JSON to stdout for structured results:

**Python**:

```python
import json
result = {'key': 'value'}
print(json.dumps(result))
```

**R**:

```r
library(jsonlite)
result <- list(key = 'value')
cat(toJSON(result))
```

**Bash**:

```bash
echo '{"key": "value"}'
```

### Input Files

Access input files by their IDs. The system downloads them to the working directory:

```python
# Input file IDs: [1, 2, 3]
# Files available in working directory
df1 = pd.read_csv('file1.csv')
df2 = pd.read_csv('file2.csv')
```

### Error Handling

Use appropriate exit codes:

```python
import sys

try:
    # Your code
    pass
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    sys.exit(1)
```

## Configuration

### Resource Limits

Configure resource limits per script:

- **timeout**: Maximum execution time (1-3600 seconds)
- **max_memory**: Maximum memory usage (100-8192 MB)

### Visibility Levels

- **private**: Only visible to owner
- **project**: Visible to project members (future)
- **public**: Visible to all users

### Script Verification

Platform administrators can mark scripts as "verified" to indicate they've been reviewed for:

- Security (no malicious code)
- Functionality (works as described)
- Best practices (follows guidelines)

## Best Practices

### Security

1. **Input Validation**: Always validate input parameters
2. **Resource Limits**: Set appropriate timeout and memory limits
3. **Error Handling**: Handle errors gracefully with clear messages
4. **No Hardcoded Credentials**: Never hardcode API keys or passwords

### Performance

1. **Efficient Algorithms**: Use efficient algorithms for large datasets
2. **Memory Management**: Clean up large objects when done
3. **Parallel Processing**: Use multiprocessing for CPU-intensive tasks
4. **Incremental Output**: Write results incrementally for long runs

### Maintainability

1. **Documentation**: Include clear description and parameter docs
2. **Comments**: Comment complex logic in script
3. **Versioning**: Use tags to indicate versions (e.g., "v1.0")
4. **Requirements**: List all dependencies explicitly

### Testing

1. **Local Testing**: Test scripts locally before uploading
2. **Parameter Validation**: Test with various parameter combinations
3. **Edge Cases**: Test with empty inputs, large files, etc.
4. **Error Scenarios**: Test error handling

## Troubleshooting

### Common Issues

**Script fails with timeout error**:

- Increase timeout value
- Optimize script performance
- Process data in chunks

**Script fails with memory error**:

- Increase max_memory
- Process data in smaller batches
- Use memory-efficient data structures

**Parameters not validated correctly**:

- Check JSON Schema syntax
- Use online JSON Schema validators
- Test parameters with validate endpoint

**Output files not appearing**:

- Ensure files are saved to `output/` directory
- Check file permissions
- Verify file creation in script logs

**Script fails immediately**:

- Check script syntax
- Verify all requirements are available
- Check error logs in execution details

### Debugging

1. **Check Execution Logs**:

   - View output_text for stdout
   - View error_text for stderr
   - Check exit_code

2. **Validate Parameters**:

   ```bash
   POST /api/custom-scripts/scripts/{id}/validate-parameters
   ```

3. **Test Locally**:

   ```bash
   # Create params file
   echo '{"input_file": "test.csv"}' > params.json

   # Run script
   python script.py --params params.json
   ```

4. **Check Requirements**:
   ```bash
   pip list  # Python
   installed.packages()  # R
   which tool  # Bash
   ```

## Migration Notes

### From Standalone Scripts

If migrating existing scripts:

1. **Add Parameter Parsing**: Add JSON parameter file parsing
2. **Add Output Directory**: Save results to `output/`
3. **Add Return JSON**: Print results as JSON
4. **List Requirements**: Document all dependencies
5. **Define Schema**: Create parameters_schema

### From Other Platforms

- **Galaxy**: Convert XML tools to scripts
- **Nextflow**: Extract process scripts
- **Snakemake**: Extract rule scripts

## Future Enhancements

- **Docker Containers**: Run scripts in Docker containers
- **GPU Support**: Enable GPU-accelerated scripts
- **Scheduled Executions**: Schedule periodic script runs
- **Script Marketplace**: Share scripts in marketplace
- **Version Control**: Git integration for script versions
- **Collaborative Editing**: Real-time collaborative script editing
- **Code Review**: Peer review system for public scripts
- **Performance Profiling**: Built-in profiling tools

## API Reference Summary

| Endpoint                            | Method | Description             |
| ----------------------------------- | ------ | ----------------------- |
| `/scripts`                          | POST   | Create script           |
| `/scripts/{id}`                     | GET    | Get script details      |
| `/scripts`                          | GET    | List scripts            |
| `/scripts/{id}`                     | PUT    | Update script           |
| `/scripts/{id}`                     | DELETE | Delete script           |
| `/scripts/{id}/execute`             | POST   | Execute script          |
| `/scripts/{id}/validate-parameters` | POST   | Validate parameters     |
| `/executions/{id}`                  | GET    | Get execution           |
| `/executions`                       | GET    | List executions         |
| `/executions/{id}/cancel`           | POST   | Cancel execution        |
| `/languages`                        | GET    | Get supported languages |
| `/categories`                       | GET    | Get categories          |
| `/templates`                        | GET    | Get templates           |
