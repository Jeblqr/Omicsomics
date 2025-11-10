"""
Visualization Workspace API

Provides endpoints for dashboard and panel management:
- CRUD operations for dashboards
- Panel management within dashboards
- Template operations
- Dashboard export/import
"""

from typing import Dict, List, Optional
from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from app.database import get_db
from app.services.visualization_workspace import VisualizationWorkspaceService

router = APIRouter(prefix="/api/viz-workspace")


# ==================== Pydantic Models ====================


class DashboardCreate(BaseModel):
    """Request to create a dashboard."""

    project_id: int
    name: str
    description: str = ""
    layout: Optional[Dict] = None
    metadata: Optional[Dict] = None
    is_public: bool = False


class DashboardUpdate(BaseModel):
    """Request to update a dashboard."""

    name: Optional[str] = None
    description: Optional[str] = None
    layout: Optional[Dict] = None
    metadata: Optional[Dict] = None
    is_public: Optional[bool] = None


class DashboardResponse(BaseModel):
    """Dashboard information."""

    id: int
    project_id: int
    name: str
    description: str
    layout: Dict
    metadata: Dict
    is_public: bool
    is_template: bool
    owner_id: int
    panel_count: int
    created_at: str
    updated_at: str


class PanelCreate(BaseModel):
    """Request to create a panel."""

    panel_key: str = Field(..., description="Unique key within dashboard")
    title: str
    description: str = ""
    viz_type: str = Field(
        ...,
        description="Visualization type: line, bar, scatter, heatmap, table, metric",
    )
    data_source: Dict = Field(..., description="Data source configuration")
    viz_config: Dict = Field(
        default_factory=dict, description="Visualization configuration"
    )
    position: Dict = Field(..., description="Panel position: {x, y, w, h}")
    auto_refresh: bool = False
    refresh_interval: Optional[int] = None


class PanelUpdate(BaseModel):
    """Request to update a panel."""

    title: Optional[str] = None
    description: Optional[str] = None
    viz_type: Optional[str] = None
    data_source: Optional[Dict] = None
    viz_config: Optional[Dict] = None
    position: Optional[Dict] = None
    auto_refresh: Optional[bool] = None
    refresh_interval: Optional[int] = None


class PanelResponse(BaseModel):
    """Panel information."""

    id: int
    dashboard_id: int
    panel_key: str
    title: str
    description: str
    viz_type: str
    data_source: Dict
    viz_config: Dict
    position: Dict
    auto_refresh: bool
    refresh_interval: Optional[int]
    created_at: str
    updated_at: str


class PositionUpdate(BaseModel):
    """Request to update panel positions."""

    positions: Dict[str, Dict] = Field(..., description="Map of panel_key to position")


class DuplicateRequest(BaseModel):
    """Request to duplicate a dashboard."""

    new_name: str


class TemplateCreateRequest(BaseModel):
    """Request to create template from dashboard."""

    template_name: str


class DashboardFromTemplateRequest(BaseModel):
    """Request to create dashboard from template."""

    template_id: int
    name: str
    project_id: int


# ==================== Helper Functions ====================


def format_dashboard_response(dashboard, include_panels: bool = False) -> Dict:
    """Format dashboard for response."""
    response = {
        "id": dashboard.id,
        "project_id": dashboard.project_id,
        "name": dashboard.name,
        "description": dashboard.description,
        "layout": dashboard.layout,
        "metadata": dashboard.metadata,
        "is_public": dashboard.is_public,
        "is_template": dashboard.is_template,
        "owner_id": dashboard.owner_id,
        "panel_count": len(dashboard.panels),
        "created_at": dashboard.created_at.isoformat(),
        "updated_at": dashboard.updated_at.isoformat(),
    }

    if include_panels:
        response["panels"] = [format_panel_response(p) for p in dashboard.panels]

    return response


def format_panel_response(panel) -> Dict:
    """Format panel for response."""
    return {
        "id": panel.id,
        "dashboard_id": panel.dashboard_id,
        "panel_key": panel.panel_key,
        "title": panel.title,
        "description": panel.description,
        "viz_type": panel.viz_type,
        "data_source": panel.data_source,
        "viz_config": panel.viz_config,
        "position": panel.position,
        "auto_refresh": panel.auto_refresh,
        "refresh_interval": panel.refresh_interval,
        "created_at": panel.created_at.isoformat(),
        "updated_at": panel.updated_at.isoformat(),
    }


# ==================== Dashboard Endpoints ====================


@router.post(
    "/dashboards", response_model=DashboardResponse, status_code=status.HTTP_201_CREATED
)
def create_dashboard(
    request: DashboardCreate,
    db: Session = Depends(get_db),
    owner_id: int = 1,  # TODO: Get from auth
):
    """Create a new dashboard."""
    service = VisualizationWorkspaceService(db)

    try:
        dashboard = service.create_dashboard(
            project_id=request.project_id,
            name=request.name,
            description=request.description,
            owner_id=owner_id,
            layout=request.layout,
            metadata=request.metadata,
            is_public=request.is_public,
        )
        return format_dashboard_response(dashboard)
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e)
        )


@router.get("/dashboards/{dashboard_id}", response_model=DashboardResponse)
def get_dashboard(
    dashboard_id: int,
    include_panels: bool = True,
    db: Session = Depends(get_db),
):
    """Get dashboard by ID."""
    service = VisualizationWorkspaceService(db)
    dashboard = service.get_dashboard(dashboard_id)

    if not dashboard:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Dashboard {dashboard_id} not found",
        )

    return format_dashboard_response(dashboard, include_panels=include_panels)


@router.get("/dashboards", response_model=List[DashboardResponse])
def list_dashboards(
    project_id: Optional[int] = None,
    owner_id: Optional[int] = None,
    is_public: Optional[bool] = None,
    is_template: Optional[bool] = None,
    search: Optional[str] = None,
    db: Session = Depends(get_db),
):
    """List dashboards with optional filters."""
    service = VisualizationWorkspaceService(db)
    dashboards = service.list_dashboards(
        project_id=project_id,
        owner_id=owner_id,
        is_public=is_public,
        is_template=is_template,
        search=search,
    )
    return [format_dashboard_response(d) for d in dashboards]


@router.put("/dashboards/{dashboard_id}", response_model=DashboardResponse)
def update_dashboard(
    dashboard_id: int,
    request: DashboardUpdate,
    db: Session = Depends(get_db),
):
    """Update dashboard properties."""
    service = VisualizationWorkspaceService(db)
    dashboard = service.update_dashboard(
        dashboard_id=dashboard_id,
        name=request.name,
        description=request.description,
        layout=request.layout,
        metadata=request.metadata,
        is_public=request.is_public,
    )

    if not dashboard:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Dashboard {dashboard_id} not found",
        )

    return format_dashboard_response(dashboard)


@router.delete("/dashboards/{dashboard_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_dashboard(
    dashboard_id: int,
    db: Session = Depends(get_db),
):
    """Delete a dashboard."""
    service = VisualizationWorkspaceService(db)
    success = service.delete_dashboard(dashboard_id)

    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Dashboard {dashboard_id} not found",
        )


@router.post("/dashboards/{dashboard_id}/duplicate", response_model=DashboardResponse)
def duplicate_dashboard(
    dashboard_id: int,
    request: DuplicateRequest,
    db: Session = Depends(get_db),
    owner_id: int = 1,  # TODO: Get from auth
):
    """Duplicate a dashboard."""
    service = VisualizationWorkspaceService(db)
    new_dashboard = service.duplicate_dashboard(
        dashboard_id, request.new_name, owner_id
    )

    if not new_dashboard:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Dashboard {dashboard_id} not found",
        )

    return format_dashboard_response(new_dashboard, include_panels=True)


@router.get("/dashboards/{dashboard_id}/stats")
def get_dashboard_statistics(
    dashboard_id: int,
    db: Session = Depends(get_db),
):
    """Get dashboard statistics."""
    service = VisualizationWorkspaceService(db)
    stats = service.get_dashboard_statistics(dashboard_id)

    if not stats:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Dashboard {dashboard_id} not found",
        )

    return stats


# ==================== Panel Endpoints ====================


@router.post(
    "/dashboards/{dashboard_id}/panels",
    response_model=PanelResponse,
    status_code=status.HTTP_201_CREATED,
)
def create_panel(
    dashboard_id: int,
    request: PanelCreate,
    db: Session = Depends(get_db),
):
    """Add a panel to a dashboard."""
    service = VisualizationWorkspaceService(db)
    panel = service.add_panel(
        dashboard_id=dashboard_id,
        panel_key=request.panel_key,
        title=request.title,
        description=request.description,
        viz_type=request.viz_type,
        data_source=request.data_source,
        viz_config=request.viz_config,
        position=request.position,
        auto_refresh=request.auto_refresh,
        refresh_interval=request.refresh_interval,
    )

    if not panel:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Failed to create panel (dashboard not found or duplicate panel_key)",
        )

    return format_panel_response(panel)


@router.get("/panels/{panel_id}", response_model=PanelResponse)
def get_panel(
    panel_id: int,
    db: Session = Depends(get_db),
):
    """Get panel by ID."""
    service = VisualizationWorkspaceService(db)
    panel = service.get_panel(panel_id)

    if not panel:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Panel {panel_id} not found"
        )

    return format_panel_response(panel)


@router.get("/dashboards/{dashboard_id}/panels", response_model=List[PanelResponse])
def list_panels(
    dashboard_id: int,
    db: Session = Depends(get_db),
):
    """List all panels in a dashboard."""
    service = VisualizationWorkspaceService(db)
    panels = service.list_panels(dashboard_id)
    return [format_panel_response(p) for p in panels]


@router.put("/panels/{panel_id}", response_model=PanelResponse)
def update_panel(
    panel_id: int,
    request: PanelUpdate,
    db: Session = Depends(get_db),
):
    """Update panel properties."""
    service = VisualizationWorkspaceService(db)
    panel = service.update_panel(
        panel_id=panel_id,
        title=request.title,
        description=request.description,
        viz_type=request.viz_type,
        data_source=request.data_source,
        viz_config=request.viz_config,
        position=request.position,
        auto_refresh=request.auto_refresh,
        refresh_interval=request.refresh_interval,
    )

    if not panel:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Panel {panel_id} not found"
        )

    return format_panel_response(panel)


@router.delete("/panels/{panel_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_panel(
    panel_id: int,
    db: Session = Depends(get_db),
):
    """Delete a panel."""
    service = VisualizationWorkspaceService(db)
    success = service.delete_panel(panel_id)

    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Panel {panel_id} not found"
        )


@router.put("/dashboards/{dashboard_id}/panels/positions")
def update_panel_positions(
    dashboard_id: int,
    request: PositionUpdate,
    db: Session = Depends(get_db),
):
    """Batch update panel positions."""
    service = VisualizationWorkspaceService(db)
    success = service.update_panel_positions(dashboard_id, request.positions)

    return {"message": f"Updated positions for {len(request.positions)} panels"}


# ==================== Template Endpoints ====================


@router.get("/templates", response_model=List[DashboardResponse])
def list_templates(
    db: Session = Depends(get_db),
):
    """List all available dashboard templates."""
    service = VisualizationWorkspaceService(db)
    templates = service.list_templates()
    return [format_dashboard_response(t, include_panels=True) for t in templates]


@router.post(
    "/dashboards/{dashboard_id}/make-template", response_model=DashboardResponse
)
def create_template_from_dashboard(
    dashboard_id: int,
    request: TemplateCreateRequest,
    db: Session = Depends(get_db),
):
    """Convert a dashboard to a template."""
    service = VisualizationWorkspaceService(db)
    template = service.create_template_from_dashboard(
        dashboard_id, request.template_name
    )

    if not template:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Dashboard {dashboard_id} not found",
        )

    return format_dashboard_response(template, include_panels=True)


@router.post("/templates/create-dashboard", response_model=DashboardResponse)
def create_dashboard_from_template(
    request: DashboardFromTemplateRequest,
    db: Session = Depends(get_db),
    owner_id: int = 1,  # TODO: Get from auth
):
    """Create a new dashboard from a template."""
    service = VisualizationWorkspaceService(db)
    dashboard = service.create_dashboard_from_template(
        template_id=request.template_id,
        name=request.name,
        project_id=request.project_id,
        owner_id=owner_id,
    )

    if not dashboard:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Template {request.template_id} not found",
        )

    return format_dashboard_response(dashboard, include_panels=True)


# ==================== Info Endpoint ====================


@router.get("/")
def visualization_workspace_info():
    """Get information about visualization workspace."""
    return {
        "service": "Visualization Workspace",
        "version": "1.0",
        "description": "Dashboard builder for multi-panel visualizations",
        "features": [
            "Create and manage dashboards",
            "Drag-and-drop panel layout",
            "Multiple visualization types",
            "Dashboard templates",
            "Auto-refresh panels",
            "Public/private sharing",
        ],
        "viz_types": [
            "line",
            "bar",
            "scatter",
            "heatmap",
            "table",
            "metric",
            "text",
        ],
        "endpoints": {
            "POST /dashboards": "Create dashboard",
            "GET /dashboards/{id}": "Get dashboard",
            "GET /dashboards": "List dashboards",
            "PUT /dashboards/{id}": "Update dashboard",
            "DELETE /dashboards/{id}": "Delete dashboard",
            "POST /dashboards/{id}/duplicate": "Duplicate dashboard",
            "POST /dashboards/{id}/panels": "Add panel",
            "GET /panels/{id}": "Get panel",
            "PUT /panels/{id}": "Update panel",
            "DELETE /panels/{id}": "Delete panel",
            "GET /templates": "List templates",
            "POST /templates/create-dashboard": "Create from template",
        },
    }
