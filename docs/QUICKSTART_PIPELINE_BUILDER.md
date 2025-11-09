# Quick Start: Using the Visual Pipeline Builder

## Starting the Application

1. **Start all services:**
   ```bash
   cd /home/jeblqr/data1/projects/Omicsomics
   docker-compose -f infrastructure/docker-compose.yml up -d
   ```

2. **Verify backend is running:**
   ```bash
   curl http://localhost:8001/api/v1/pipelines/
   ```

3. **Start frontend:**
   ```bash
   cd frontend
   npm run dev
   ```

4. **Access the application:**
   - Open browser to `http://localhost:5173`
   - Login with your credentials

## Creating Your First Pipeline

### Step 1: Navigate to Custom Pipelines
- Click **Custom Pipelines** in the left sidebar

### Step 2: Create New Pipeline
1. Click **+ Create New Pipeline** button
2. Fill in the form:
   - **Name**: "My First Pipeline"
   - **Description**: "A simple quality control workflow"
   - **Category**: Select "genomics"
   - **Make Public**: Leave unchecked

### Step 3: Build the Pipeline
1. Select node type: **Input**
2. Click **+ Add Node**
3. Drag the node to desired position
4. Repeat for more nodes:
   - Add a **Process** node (Quality Check)
   - Add a **Filter** node (Remove Low Quality)
   - Add an **Output** node (Clean Data)

### Step 4: Connect Nodes
1. Click and hold on the edge of the first node
2. Drag to the next node
3. Release to create connection
4. Repeat to connect all nodes in sequence

### Step 5: Save Pipeline
1. Click **Save Pipeline** button
2. You'll see a success message
3. Click **← Back to List** to return

## Merging Pipelines

### Step 1: Create Second Pipeline
1. Follow the same steps to create another pipeline
2. Example: "Analysis Pipeline" with Analysis and Output nodes

### Step 2: Merge Pipelines
1. Go to **Custom Pipelines**
2. Check the boxes next to the 2 pipelines you want to merge
3. Click **Merge Selected (2)** button
4. The editor opens with the merged pipeline:
   - All nodes from both pipelines
   - Automatic connection between pipelines
   - Unique node IDs (p0_n0, p1_n3, etc.)

### Step 3: Save Merged Pipeline
1. Optionally rename to "Complete Workflow"
2. Click **Save Pipeline**
3. Your merged pipeline is now ready to use

## Using Pipelines in Runs

1. Navigate to **Runs** page
2. Click **Create New Run**
3. Select your custom pipeline
4. Configure parameters
5. Submit the run

The complete pipeline configuration is stored with the run for reproducibility.

## Quick API Examples

### Get All Pipelines
```bash
# Login and get token
TOKEN=$(curl -s -X POST http://localhost:8001/api/v1/auth/login/access-token \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=your@email.com&password=yourpassword" | jq -r .access_token)

# List pipelines
curl -H "Authorization: Bearer $TOKEN" \
  http://localhost:8001/api/v1/custom-pipelines/ | jq .
```

### Create Pipeline via API
```bash
curl -X POST -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  http://localhost:8001/api/v1/custom-pipelines/ \
  -d '{
    "name": "API Pipeline",
    "description": "Created via API",
    "category": "custom",
    "is_public": false,
    "definition": {
      "nodes": [
        {
          "id": "node_1",
          "type": "input",
          "label": "Data Input",
          "data": {},
          "position": {"x": 100, "y": 100}
        },
        {
          "id": "node_2",
          "type": "output",
          "label": "Results",
          "data": {},
          "position": {"x": 300, "y": 100}
        }
      ],
      "edges": [
        {
          "id": "e1",
          "source": "node_1",
          "target": "node_2"
        }
      ],
      "parameters": {}
    }
  }' | jq .
```

### Merge Pipelines via API
```bash
curl -X POST -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  http://localhost:8001/api/v1/custom-pipelines/merge \
  -d '{
    "pipeline_ids": [1, 2]
  }' | jq .
```

## Troubleshooting

### Can't See Custom Pipelines Menu
- **Solution**: Refresh the page, ensure you're logged in

### Nodes Won't Connect
- **Solution**: 
  - Make sure both nodes are on the canvas
  - Try clicking precisely on the node edge
  - Check that you're not in read-only mode

### Save Button Not Working
- **Solution**: 
  - Ensure pipeline name is filled in
  - Check browser console for errors
  - Verify backend is running: `docker ps | grep backend`

### Merge Button Disabled
- **Solution**: 
  - You need to select at least 2 pipelines
  - Click the checkboxes next to pipeline cards

### API Returns 404
- **Solution**:
  - Verify backend is running: `docker logs infrastructure-backend-1`
  - Check migration applied: `docker exec infrastructure-backend-1 alembic current`
  - Restart backend: `docker restart infrastructure-backend-1`

### Frontend Not Loading
- **Solution**:
  ```bash
  cd frontend
  npm install  # Reinstall dependencies including reactflow
  npm run dev  # Restart dev server
  ```

## Common Tasks

### Delete All Test Pipelines
```bash
# Get all pipeline IDs
curl -H "Authorization: Bearer $TOKEN" \
  http://localhost:8001/api/v1/custom-pipelines/ | jq '.[].id'

# Delete each one
curl -X DELETE -H "Authorization: Bearer $TOKEN" \
  http://localhost:8001/api/v1/custom-pipelines/1

curl -X DELETE -H "Authorization: Bearer $TOKEN" \
  http://localhost:8001/api/v1/custom-pipelines/2
```

### Export Pipeline Definition
```bash
# Get pipeline and save to file
curl -H "Authorization: Bearer $TOKEN" \
  http://localhost:8001/api/v1/custom-pipelines/1 | \
  jq '.definition' > my_pipeline.json
```

### Import Pipeline Definition
```bash
# Load from file and create new pipeline
curl -X POST -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  http://localhost:8001/api/v1/custom-pipelines/ \
  -d "{
    \"name\": \"Imported Pipeline\",
    \"description\": \"Loaded from file\",
    \"category\": \"custom\",
    \"is_public\": false,
    \"definition\": $(cat my_pipeline.json)
  }"
```

### Check Database Migration Status
```bash
# See current migration version
docker exec infrastructure-backend-1 alembic current

# Should show:
# 0003_add_custom_pipelines (head)

# If not at head, upgrade:
docker exec infrastructure-backend-1 alembic upgrade head
```

### View Pipeline in Database
```bash
# Connect to PostgreSQL
docker exec -it infrastructure-postgres-1 psql -U postgres -d omicsomics

# Query pipelines
SELECT id, name, category, is_public, owner_id 
FROM custom_pipelines;

# View full definition (formatted)
SELECT id, name, jsonb_pretty(definition) 
FROM custom_pipelines 
WHERE id = 1;

# Exit psql
\q
```

## Next Steps

1. **Explore Templates**: Check out the built-in pipeline templates in the Pipelines page
2. **Merge with Templates**: Try merging your custom pipeline with a template
3. **Run Pipelines**: Execute your custom pipeline from the Runs page
4. **Share Pipelines**: Toggle "Make Public" to share with other users
5. **Advanced Features**: Check the full documentation in `docs/PIPELINE_BUILDER.md`

## Need Help?

- **User Guide**: `docs/PIPELINE_BUILDER.md`
- **Implementation Details**: `docs/IMPLEMENTATION_PIPELINE_BUILDER.md`
- **API Docs**: Visit `http://localhost:8001/docs` when backend is running
- **Logs**: `docker logs infrastructure-backend-1`
