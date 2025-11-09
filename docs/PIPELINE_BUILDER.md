# Visual Pipeline Builder - User Guide

## Overview

Omicsomics now features a powerful visual pipeline builder that allows you to create, customize, and merge analysis pipelines through an intuitive drag-and-drop interface.

## Key Features

### 1. Visual Pipeline Editor

- **Drag-and-Drop Interface**: Build pipelines visually using React Flow
- **6 Node Types**:
  - 🟢 **Input**: Data input nodes
  - 🔵 **Process**: Processing and transformation steps
  - 🟠 **Filter**: Data filtering operations
  - 🟣 **Transform**: Data transformation steps
  - 🔴 **Analysis**: Statistical and computational analysis
  - ⚫ **Output**: Results and export nodes

### 2. Pipeline Management

- **Create Custom Pipelines**: Design your own analysis workflows
- **Edit Existing Pipelines**: Modify saved pipelines
- **Delete Pipelines**: Remove unwanted pipelines
- **Public/Private Sharing**: Control pipeline visibility
- **Categories**: Organize by omics type (genomics, proteomics, metabolomics, etc.)

### 3. Pipeline Merging

- **Multi-Select**: Select 2+ pipelines to merge
- **Sequential Connection**: Automatically connects pipelines in sequence
- **Smart Node Remapping**: Avoids ID conflicts (p{i}\_n{offset} format)
- **Bridge Edges**: Creates connections between merged segments
- **Parameter Namespacing**: Isolates parameters (pipeline*{i}*{key})
- **Metadata Tracking**: Records merge source and strategy

## How to Use

### Creating a Custom Pipeline

1. Navigate to **Custom Pipelines** in the sidebar
2. Click **+ Create New Pipeline**
3. Enter pipeline details:
   - **Name**: Descriptive pipeline name
   - **Description**: Optional detailed description
   - **Category**: Choose from predefined categories
   - **Make Public**: Toggle to share with other users
4. Build your pipeline:
   - Select node type from dropdown
   - Click **+ Add Node** to add to canvas
   - Drag nodes to position them
   - Click and drag between nodes to create connections
5. Click **Save Pipeline** when finished
6. Click **← Back to List** to return

### Editing a Pipeline

1. Go to **Custom Pipelines**
2. Find your pipeline in the list
3. Click **Edit** button
4. Modify nodes, connections, or metadata
5. Click **Save Pipeline** to apply changes

### Merging Pipelines

1. Navigate to **Custom Pipelines**
2. Check the checkboxes for 2 or more pipelines
3. Click **Merge Selected (N)** button
4. Review the merged pipeline in the editor
5. Optionally rename and adjust
6. Click **Save Pipeline** to save the merged result

The merge algorithm will:

- Combine all nodes from selected pipelines
- Renumber node IDs to prevent conflicts
- Connect the last node of pipeline N to the first node of pipeline N+1
- Preserve all original edges within each pipeline segment
- Add metadata tracking the source pipelines

### Using Pipelines in Runs

1. Go to **Runs** page
2. Click **Create New Run**
3. Select a custom pipeline from the pipeline selector
4. Configure run parameters
5. Submit the run

The pipeline configuration will be saved with the run for reproducibility.

## API Endpoints

### Custom Pipelines API

```bash
# Get authentication token
TOKEN=$(curl -s -X POST http://localhost:8001/api/v1/auth/login/access-token \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=your@email.com&password=yourpassword" | jq -r .access_token)

# List all accessible pipelines (own + public)
curl -H "Authorization: Bearer $TOKEN" \
  http://localhost:8001/api/v1/custom-pipelines/

# Get specific pipeline
curl -H "Authorization: Bearer $TOKEN" \
  http://localhost:8001/api/v1/custom-pipelines/1

# Create new pipeline
curl -X POST -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  http://localhost:8001/api/v1/custom-pipelines/ \
  -d '{
    "name": "My Pipeline",
    "description": "Pipeline description",
    "category": "genomics",
    "is_public": false,
    "definition": {
      "nodes": [...],
      "edges": [...],
      "parameters": {}
    }
  }'

# Update pipeline
curl -X PUT -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  http://localhost:8001/api/v1/custom-pipelines/1 \
  -d '{...}'

# Delete pipeline
curl -X DELETE -H "Authorization: Bearer $TOKEN" \
  http://localhost:8001/api/v1/custom-pipelines/1

# Merge pipelines
curl -X POST -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  http://localhost:8001/api/v1/custom-pipelines/merge \
  -d '{
    "pipeline_ids": [1, 2, 3]
  }'
```

## Pipeline Definition Schema

### Node Structure

```json
{
  "id": "node_1",
  "type": "process",
  "label": "Quality Check",
  "data": {},
  "position": { "x": 100, "y": 100 }
}
```

### Edge Structure

```json
{
  "id": "edge_1",
  "source": "node_1",
  "target": "node_2",
  "sourceHandle": null,
  "targetHandle": null
}
```

### Complete Definition

```json
{
  "nodes": [
    {
      "id": "node_1",
      "type": "input",
      "label": "Raw Data",
      "data": {},
      "position": { "x": 100, "y": 100 }
    },
    {
      "id": "node_2",
      "type": "process",
      "label": "QC",
      "data": {},
      "position": { "x": 300, "y": 100 }
    },
    {
      "id": "node_3",
      "type": "output",
      "label": "Clean Data",
      "data": {},
      "position": { "x": 500, "y": 100 }
    }
  ],
  "edges": [
    { "id": "e1", "source": "node_1", "target": "node_2" },
    { "id": "e2", "source": "node_2", "target": "node_3" }
  ],
  "parameters": {}
}
```

## Merge Algorithm Details

When merging pipelines, the system performs:

1. **Node ID Remapping**:

   - Format: `p{pipeline_index}_n{node_offset}`
   - Example: `p0_n0`, `p0_n1`, `p1_n2`, `p1_n3`

2. **Edge ID Remapping**:

   - Format: `p{pipeline_index}_e{edge_offset}`
   - Example: `p0_e0`, `p0_e1`, `p1_e2`

3. **Bridge Creation**:

   - Connects last node of pipeline N to first node of pipeline N+1
   - Format: `bridge_p{N}_to_p{N+1}_{offset}`
   - Type marked as "bridge"

4. **Parameter Namespacing**:

   - Format: `pipeline_{index}_{original_key}`
   - Prevents parameter conflicts

5. **Metadata Addition**:
   - `merged_from`: Array of source pipeline IDs
   - `merge_strategy`: "sequential"

## Database Schema

### Custom Pipelines Table

```sql
CREATE TABLE custom_pipelines (
    id SERIAL PRIMARY KEY,
    name VARCHAR NOT NULL,
    description TEXT,
    definition JSONB NOT NULL,
    category VARCHAR,
    is_public BOOLEAN DEFAULT FALSE,
    owner_id INTEGER REFERENCES users(id),
    template_id VARCHAR,
    created_at TIMESTAMP DEFAULT NOW(),
    updated_at TIMESTAMP DEFAULT NOW()
);
```

### Runs Table Enhancement

```sql
ALTER TABLE runs ADD COLUMN pipeline_config JSONB;
```

## Best Practices

1. **Naming**: Use descriptive names for pipelines and nodes
2. **Description**: Document what each pipeline does
3. **Categories**: Use appropriate categories for organization
4. **Public Sharing**: Only make pipelines public if they're complete and tested
5. **Merging**: Test individual pipelines before merging
6. **Node Count**: Keep pipelines focused (5-10 nodes is ideal)
7. **Documentation**: Add notes in node labels about what each step does

## Troubleshooting

### Pipeline Won't Save

- Ensure all required fields are filled (name, definition)
- Check that nodes have unique IDs
- Verify edges reference existing node IDs

### Merge Fails

- Ensure you've selected at least 2 pipelines
- Check that all pipelines have valid definitions
- Verify you have access to all selected pipelines

### Edges Not Connecting

- Make sure source and target nodes exist
- Try re-dragging the connection
- Check that nodes are on the canvas (not deleted)

## Future Enhancements

Planned features:

- [ ] Real-time collaboration on pipeline editing
- [ ] Pipeline templates library
- [ ] Version control for pipelines
- [ ] Pipeline validation and testing
- [ ] Auto-layout for better visualization
- [ ] Export/import pipeline definitions
- [ ] Pipeline execution preview
- [ ] Custom node types with parameters
- [ ] Conditional branching in pipelines
- [ ] Parallel execution paths

## Support

For issues or questions:

- Check the API documentation at `/docs` endpoint
- Review error messages in browser console
- Check backend logs: `docker logs infrastructure-backend-1`
- Verify database migration applied: `docker exec infrastructure-backend-1 alembic current`
