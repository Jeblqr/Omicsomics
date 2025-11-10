# Visualization Workspace Documentation

## Overview

The **Visualization Workspace** is a flexible dashboard builder that allows users to create, organize, and share multi-panel visualizations. It provides a drag-and-drop interface for arranging visualization panels, supporting multiple chart types, and enabling data-driven insights.

## Features

### Core Capabilities

1. **Dashboard Management**

   - Create unlimited dashboards per project
   - Save and load dashboard configurations
   - Duplicate dashboards for quick replication
   - Public/private sharing options
   - Template system for reusable layouts

2. **Panel Types**

   - **Line Chart** 📈: Time series and trend visualization
   - **Bar Chart** 📊: Categorical comparisons
   - **Scatter Plot** ⚫: Correlation analysis
   - **Heatmap** 🟥: Matrix visualization
   - **Table** 📋: Tabular data display
   - **Metric** 🔢: Single value display
   - **Text** 📝: Annotations and descriptions

3. **Layout System**

   - Responsive grid layout (12-column)
   - Drag-and-drop panel positioning
   - Resize panels dynamically
   - Automatic layout adjustment

4. **Data Integration**

   - Link panels to datasets
   - Connect to data files
   - Real-time data updates
   - Auto-refresh capabilities

5. **Collaboration**
   - Share dashboards publicly
   - Export dashboard configurations
   - Create templates from dashboards
   - Import shared dashboards

## Architecture

### Database Models

#### Dashboard Model

```python
class Dashboard(Base):
    id: int
    name: str
    description: str
    layout: dict              # Grid configuration
    metadata: dict            # Custom properties
    is_public: bool
    is_template: bool
    project_id: int
    owner_id: int
    panels: List[Panel]       # One-to-many relationship
    created_at: datetime
    updated_at: datetime
```

#### Panel Model

```python
class Panel(Base):
    id: int
    dashboard_id: int
    panel_key: str            # Unique within dashboard
    title: str
    description: str
    viz_type: str             # line, bar, scatter, etc.
    data_source: dict         # File, dataset, or URL
    viz_config: dict          # Chart options
    position: dict            # {x, y, w, h}
    auto_refresh: bool
    refresh_interval: int     # Seconds
    created_at: datetime
    updated_at: datetime
```

### Service Layer

**VisualizationWorkspaceService** (`backend/app/services/visualization_workspace.py`)

Main operations:

- `create_dashboard()` - Create new dashboard
- `update_dashboard()` - Modify dashboard properties
- `duplicate_dashboard()` - Clone with all panels
- `add_panel()` - Add visualization panel
- `update_panel()` - Modify panel configuration
- `update_panel_positions()` - Batch update layout
- `create_template_from_dashboard()` - Convert to template
- `create_dashboard_from_template()` - Instantiate template

## API Reference

### Base URL

`/api/viz-workspace`

### Dashboard Endpoints

#### Create Dashboard

```http
POST /dashboards
Content-Type: application/json

{
  "project_id": 1,
  "name": "RNA-seq Analysis Dashboard",
  "description": "Comprehensive view of RNA-seq results",
  "layout": {
    "cols": 12,
    "rowHeight": 100,
    "breakpoints": {"lg": 1200, "md": 996, "sm": 768}
  },
  "metadata": {},
  "is_public": false
}

Response:
{
  "id": 1,
  "name": "RNA-seq Analysis Dashboard",
  "panel_count": 0,
  "created_at": "2025-11-10T10:00:00Z",
  ...
}
```

#### Get Dashboard

```http
GET /dashboards/{dashboard_id}?include_panels=true

Response:
{
  "id": 1,
  "name": "RNA-seq Analysis Dashboard",
  "description": "...",
  "layout": {...},
  "panels": [
    {
      "id": 1,
      "panel_key": "panel_1",
      "title": "Expression Heatmap",
      "viz_type": "heatmap",
      ...
    }
  ],
  ...
}
```

#### List Dashboards

```http
GET /dashboards?project_id=1&search=rna

Response:
[
  {
    "id": 1,
    "name": "RNA-seq Analysis Dashboard",
    "panel_count": 5,
    ...
  }
]
```

#### Update Dashboard

```http
PUT /dashboards/{dashboard_id}
Content-Type: application/json

{
  "name": "Updated Dashboard Name",
  "layout": {...},
  "is_public": true
}
```

#### Delete Dashboard

```http
DELETE /dashboards/{dashboard_id}

Response: 204 No Content
```

#### Duplicate Dashboard

```http
POST /dashboards/{dashboard_id}/duplicate
Content-Type: application/json

{
  "new_name": "Copy of RNA-seq Dashboard"
}

Response:
{
  "id": 2,
  "name": "Copy of RNA-seq Dashboard",
  "panels": [...],
  ...
}
```

#### Get Dashboard Statistics

```http
GET /dashboards/{dashboard_id}/stats

Response:
{
  "dashboard_id": 1,
  "panel_count": 5,
  "viz_types": {
    "heatmap": 2,
    "line": 2,
    "metric": 1
  },
  "auto_refresh_panels": 2,
  "created_at": "2025-11-10T10:00:00Z",
  "updated_at": "2025-11-10T11:30:00Z"
}
```

### Panel Endpoints

#### Add Panel

```http
POST /dashboards/{dashboard_id}/panels
Content-Type: application/json

{
  "panel_key": "expression_heatmap",
  "title": "Gene Expression Heatmap",
  "description": "Top 50 differentially expressed genes",
  "viz_type": "heatmap",
  "data_source": {
    "type": "dataset",
    "dataset_id": 123,
    "file_path": "/data/expression_matrix.csv"
  },
  "viz_config": {
    "colormap": "RdBu",
    "show_legend": true,
    "normalize": "row"
  },
  "position": {
    "x": 0,
    "y": 0,
    "w": 6,
    "h": 4
  },
  "auto_refresh": false,
  "refresh_interval": null
}

Response:
{
  "id": 1,
  "dashboard_id": 1,
  "panel_key": "expression_heatmap",
  "title": "Gene Expression Heatmap",
  ...
}
```

#### Get Panel

```http
GET /panels/{panel_id}
```

#### List Dashboard Panels

```http
GET /dashboards/{dashboard_id}/panels

Response:
[
  {
    "id": 1,
    "panel_key": "expression_heatmap",
    "title": "Gene Expression Heatmap",
    ...
  },
  ...
]
```

#### Update Panel

```http
PUT /panels/{panel_id}
Content-Type: application/json

{
  "title": "Updated Panel Title",
  "viz_config": {...},
  "position": {"x": 6, "y": 0, "w": 6, "h": 4}
}
```

#### Delete Panel

```http
DELETE /panels/{panel_id}

Response: 204 No Content
```

#### Batch Update Panel Positions

```http
PUT /dashboards/{dashboard_id}/panels/positions
Content-Type: application/json

{
  "positions": {
    "expression_heatmap": {"x": 0, "y": 0, "w": 6, "h": 4},
    "volcano_plot": {"x": 6, "y": 0, "w": 6, "h": 4},
    "gene_counts": {"x": 0, "y": 4, "w": 12, "h": 2}
  }
}
```

### Template Endpoints

#### List Templates

```http
GET /templates

Response:
[
  {
    "id": 10,
    "name": "RNA-seq Analysis Template",
    "is_template": true,
    "is_public": true,
    "panel_count": 5,
    "panels": [...]
  }
]
```

#### Create Template from Dashboard

```http
POST /dashboards/{dashboard_id}/make-template
Content-Type: application/json

{
  "template_name": "My Analysis Template"
}
```

#### Create Dashboard from Template

```http
POST /templates/create-dashboard
Content-Type: application/json

{
  "template_id": 10,
  "name": "New Analysis Dashboard",
  "project_id": 1
}

Response:
{
  "id": 15,
  "name": "New Analysis Dashboard",
  "panels": [...],
  ...
}
```

## Usage Scenarios

### Scenario 1: Create RNA-seq Analysis Dashboard

**User Story**: Researcher wants to visualize RNA-seq results in a single dashboard.

**Steps**:

1. **Create Dashboard**

```javascript
const dashboard = await axios.post("/api/viz-workspace/dashboards", {
  project_id: 1,
  name: "RNA-seq QC Dashboard",
  description: "Quality control metrics for RNA-seq experiment",
});
```

2. **Add Expression Heatmap**

```javascript
await axios.post(`/api/viz-workspace/dashboards/${dashboard.data.id}/panels`, {
  panel_key: "expression_heatmap",
  title: "Top 50 DEGs",
  viz_type: "heatmap",
  data_source: {
    type: "dataset",
    dataset_id: 123,
  },
  position: { x: 0, y: 0, w: 6, h: 4 },
});
```

3. **Add Volcano Plot**

```javascript
await axios.post(`/api/viz-workspace/dashboards/${dashboard.data.id}/panels`, {
  panel_key: "volcano_plot",
  title: "Differential Expression",
  viz_type: "scatter",
  position: { x: 6, y: 0, w: 6, h: 4 },
});
```

4. **Add Sample Metrics**

```javascript
await axios.post(`/api/viz-workspace/dashboards/${dashboard.data.id}/panels`, {
  panel_key: "sample_counts",
  title: "Sample Statistics",
  viz_type: "table",
  position: { x: 0, y: 4, w: 12, h: 3 },
});
```

5. **Save Dashboard**
   Dashboard auto-saves on each operation.

### Scenario 2: Create and Use Template

**User Story**: Team wants to standardize QC dashboards across projects.

**Steps**:

1. **Create Master Dashboard**

   - Create dashboard with standard QC panels
   - Configure all visualizations
   - Test with sample data

2. **Convert to Template**

```javascript
const template = await axios.post(
  `/api/viz-workspace/dashboards/${dashboardId}/make-template`,
  { template_name: "Standard RNA-seq QC" }
);
```

3. **Use Template for New Project**

```javascript
const newDashboard = await axios.post(
  "/api/viz-workspace/templates/create-dashboard",
  {
    template_id: template.data.id,
    name: "Project XYZ - RNA-seq QC",
    project_id: 5,
  }
);
```

4. **Customize for Project**
   - Update panel data sources
   - Adjust configurations
   - Add project-specific panels

### Scenario 3: Real-time Monitoring Dashboard

**User Story**: Monitor pipeline execution with auto-refreshing metrics.

**Steps**:

1. **Create Monitoring Dashboard**

```javascript
const dashboard = await axios.post("/api/viz-workspace/dashboards", {
  name: "Pipeline Monitor",
  description: "Real-time pipeline execution metrics",
});
```

2. **Add Auto-refresh Metric Panels**

```javascript
// Success rate metric
await axios.post(`/api/viz-workspace/dashboards/${dashboard.data.id}/panels`, {
  panel_key: "success_rate",
  title: "Success Rate",
  viz_type: "metric",
  data_source: {
    type: "api",
    url: "/api/runs/stats/success-rate",
  },
  auto_refresh: true,
  refresh_interval: 30, // 30 seconds
  position: { x: 0, y: 0, w: 3, h: 2 },
});

// Active runs count
await axios.post(`/api/viz-workspace/dashboards/${dashboard.data.id}/panels`, {
  panel_key: "active_runs",
  title: "Active Runs",
  viz_type: "metric",
  auto_refresh: true,
  refresh_interval: 10,
  position: { x: 3, y: 0, w: 3, h: 2 },
});
```

3. **Add Run Timeline**

```javascript
await axios.post(`/api/viz-workspace/dashboards/${dashboard.data.id}/panels`, {
  panel_key: "run_timeline",
  title: "Recent Runs",
  viz_type: "line",
  auto_refresh: true,
  refresh_interval: 60,
  position: { x: 0, y: 2, w: 12, h: 4 },
});
```

### Scenario 4: Share Dashboard Publicly

**User Story**: Researcher wants to share results dashboard with collaborators.

**Steps**:

1. **Create Dashboard** (standard creation)

2. **Configure for Sharing**

```javascript
await axios.put(`/api/viz-workspace/dashboards/${dashboardId}`, {
  is_public: true,
  metadata: {
    shared_with: ["team_a", "team_b"],
    publish_date: new Date().toISOString(),
  },
});
```

3. **Generate Shareable Link**

   - System generates unique URL
   - Share link with collaborators
   - Track access in metadata

4. **Collaborators View Dashboard**
   - Read-only access for viewers
   - Owner can update anytime
   - Changes reflect immediately

## Frontend Integration

### Component Structure

**VisualizationWorkspacePage** (`frontend/src/pages/VisualizationWorkspacePage.tsx`)

Key features:

- AppBar with dashboard controls
- Grid layout for panels
- Create/Edit dialogs
- Load dashboard picker
- Settings menu

### Usage Example

```tsx
import VisualizationWorkspacePage from "./pages/VisualizationWorkspacePage";

// In router configuration
<Route path="/workspace" element={<VisualizationWorkspacePage />} />;
```

### State Management

```typescript
interface WorkspaceState {
  dashboards: Dashboard[]; // Available dashboards
  currentDashboard: Dashboard | null; // Active dashboard
  panels: Panel[]; // Current dashboard panels
  loading: boolean;
  error: string | null;
}
```

### Key Operations

```typescript
// Create dashboard
const createDashboard = async (data: DashboardCreate) => {
  const response = await axios.post("/api/viz-workspace/dashboards", data);
  setCurrentDashboard(response.data);
};

// Add panel
const addPanel = async (dashboardId: number, panel: PanelCreate) => {
  const response = await axios.post(
    `/api/viz-workspace/dashboards/${dashboardId}/panels`,
    panel
  );
  setPanels([...panels, response.data]);
};

// Update layout
const updateLayout = async (positions: Record<string, Position>) => {
  await axios.put(
    `/api/viz-workspace/dashboards/${dashboardId}/panels/positions`,
    { positions }
  );
};
```

## Configuration

### Panel Configuration Examples

#### Heatmap Panel

```json
{
  "viz_type": "heatmap",
  "viz_config": {
    "colormap": "RdBu",
    "normalize": "row",
    "cluster_rows": true,
    "cluster_cols": false,
    "show_legend": true,
    "show_row_names": true,
    "show_col_names": true,
    "font_size": 10
  }
}
```

#### Line Chart Panel

```json
{
  "viz_type": "line",
  "viz_config": {
    "x_column": "time",
    "y_columns": ["value1", "value2"],
    "line_style": "solid",
    "marker_style": "circle",
    "legend_position": "top-right",
    "grid": true,
    "title": "Time Series Plot"
  }
}
```

#### Metric Panel

```json
{
  "viz_type": "metric",
  "viz_config": {
    "value_format": ".2f",
    "prefix": "$",
    "suffix": "",
    "color": "#4CAF50",
    "comparison": {
      "enabled": true,
      "previous_value": 85.2,
      "show_trend": true
    }
  }
}
```

### Layout Configuration

```json
{
  "cols": 12,
  "rowHeight": 100,
  "breakpoints": {
    "lg": 1200,
    "md": 996,
    "sm": 768,
    "xs": 480
  },
  "margin": [10, 10],
  "containerPadding": [10, 10],
  "responsive": true
}
```

## Best Practices

### 1. Dashboard Organization

**DO**:

- Use clear, descriptive dashboard names
- Group related visualizations
- Maintain consistent panel sizes
- Add descriptions to panels
- Use templates for standard layouts

**DON'T**:

- Create too many small panels (cluttered)
- Mix unrelated visualizations
- Use generic names like "Dashboard 1"
- Overload with too many panels (>12)

### 2. Panel Configuration

**DO**:

- Configure appropriate refresh intervals
- Use panel-specific data sources
- Add clear titles and descriptions
- Choose visualization types suited to data
- Test panel configurations before saving

**DON'T**:

- Set refresh intervals too short (<10s)
- Reuse same data source unnecessarily
- Leave default panel titles
- Force data into unsuitable viz types

### 3. Performance

**DO**:

- Limit auto-refresh panels
- Use appropriate data sampling
- Cache frequently accessed data
- Optimize queries for panel data
- Test with production-size data

**DON'T**:

- Refresh all panels simultaneously
- Load entire datasets into panels
- Ignore loading states
- Create circular data dependencies

### 4. Collaboration

**DO**:

- Create templates for common patterns
- Document dashboard purpose
- Use public sharing judiciously
- Maintain version history
- Test shared dashboards as viewer

**DON'T**:

- Share sensitive data publicly
- Create too many templates
- Modify templates without communication
- Delete dashboards without backup

## Troubleshooting

### Issue: Panel Not Updating

**Symptoms**:

- Panel shows stale data
- Auto-refresh not working

**Solutions**:

```javascript
// Check panel configuration
const panel = await axios.get(`/api/viz-workspace/panels/${panelId}`);
console.log("Auto-refresh:", panel.data.auto_refresh);
console.log("Interval:", panel.data.refresh_interval);

// Force panel update
await axios.put(`/api/viz-workspace/panels/${panelId}`, {
  auto_refresh: true,
  refresh_interval: 30,
});

// Check data source
console.log("Data source:", panel.data.data_source);
```

### Issue: Layout Not Saving

**Symptoms**:

- Panel positions reset after reload
- Layout changes not persisted

**Solutions**:

```javascript
// Verify positions update
const positions = {
  panel_1: { x: 0, y: 0, w: 6, h: 4 },
  panel_2: { x: 6, y: 0, w: 6, h: 4 },
};

const response = await axios.put(
  `/api/viz-workspace/dashboards/${dashboardId}/panels/positions`,
  { positions }
);

// Check dashboard update
await axios.put(`/api/viz-workspace/dashboards/${dashboardId}`, {
  layout: newLayoutConfig,
});
```

### Issue: Template Not Applying

**Symptoms**:

- Dashboard created from template is empty
- Panels not copied

**Solutions**:

```javascript
// Verify template exists
const template = await axios.get(`/api/viz-workspace/dashboards/${templateId}`);
console.log("Is template:", template.data.is_template);
console.log("Panel count:", template.data.panel_count);

// Create with include_panels
const newDashboard = await axios.post(
  "/api/viz-workspace/templates/create-dashboard",
  {
    template_id: templateId,
    name: "New Dashboard",
    project_id: projectId,
  }
);

// Verify panels copied
const panels = await axios.get(
  `/api/viz-workspace/dashboards/${newDashboard.data.id}/panels`
);
console.log("Panels copied:", panels.data.length);
```

## Future Enhancements

1. **Advanced Visualizations**

   - 3D plots
   - Network diagrams
   - Sankey diagrams
   - Geographic maps

2. **Real-time Collaboration**

   - Multi-user editing
   - Live cursor sharing
   - Comment system
   - Change history

3. **Data Connections**

   - Direct database queries
   - API endpoints
   - Streaming data sources
   - Cloud storage integration

4. **Export Options**

   - PDF export
   - PNG/SVG export
   - Interactive HTML
   - PowerPoint slides

5. **Advanced Features**

   - Conditional formatting
   - Drill-down capabilities
   - Cross-panel filtering
   - Custom JavaScript panels

6. **Performance**
   - Virtual scrolling for large dashboards
   - Lazy loading panels
   - Caching strategies
   - Progressive rendering

## References

- [Grid Layout System](https://github.com/react-grid-layout/react-grid-layout)
- [Chart.js Documentation](https://www.chartjs.org/)
- [Material-UI Grid](https://mui.com/components/grid/)
- [Dashboard Design Patterns](https://www.nngroup.com/articles/dashboard-design/)

---

**Version**: 1.0  
**Last Updated**: 2025-11-10  
**Author**: Omicsomics Development Team
