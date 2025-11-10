"""
Visualization Workspace Service

Handles dashboard and panel management:
- Create, read, update, delete dashboards
- Manage panels within dashboards
- Generate panel data and visualizations
- Handle dashboard templates
"""

import logging
from typing import Dict, List, Optional
from sqlalchemy.orm import Session
from sqlalchemy import and_, or_

from app.models.visualization_workspace import Dashboard, Panel

logger = logging.getLogger(__name__)


class VisualizationWorkspaceService:
    """Service for managing visualization dashboards and panels."""

    def __init__(self, db: Session):
        self.db = db

    # ==================== Dashboard Operations ====================

    def create_dashboard(
        self,
        project_id: int,
        name: str,
        description: str,
        owner_id: int,
        layout: Optional[Dict] = None,
        metadata: Optional[Dict] = None,
        is_public: bool = False,
        is_template: bool = False,
    ) -> Dashboard:
        """Create a new dashboard."""
        dashboard = Dashboard(
            project_id=project_id,
            name=name,
            description=description,
            owner_id=owner_id,
            layout=layout
            or {
                "cols": 12,
                "rowHeight": 100,
                "breakpoints": {"lg": 1200, "md": 996, "sm": 768},
            },
            metadata=metadata or {},
            is_public=is_public,
            is_template=is_template,
        )
        self.db.add(dashboard)
        self.db.commit()
        self.db.refresh(dashboard)
        logger.info(f"Created dashboard {dashboard.id}: {dashboard.name}")
        return dashboard

    def get_dashboard(self, dashboard_id: int) -> Optional[Dashboard]:
        """Get dashboard by ID."""
        return self.db.query(Dashboard).filter(Dashboard.id == dashboard_id).first()

    def list_dashboards(
        self,
        project_id: Optional[int] = None,
        owner_id: Optional[int] = None,
        is_public: Optional[bool] = None,
        is_template: Optional[bool] = None,
        search: Optional[str] = None,
    ) -> List[Dashboard]:
        """List dashboards with optional filters."""
        query = self.db.query(Dashboard)

        # Apply filters
        conditions = []
        if project_id is not None:
            conditions.append(Dashboard.project_id == project_id)
        if owner_id is not None:
            conditions.append(Dashboard.owner_id == owner_id)
        if is_public is not None:
            conditions.append(Dashboard.is_public == is_public)
        if is_template is not None:
            conditions.append(Dashboard.is_template == is_template)
        if search:
            search_pattern = f"%{search}%"
            conditions.append(
                or_(
                    Dashboard.name.ilike(search_pattern),
                    Dashboard.description.ilike(search_pattern),
                )
            )

        if conditions:
            query = query.filter(and_(*conditions))

        return query.order_by(Dashboard.updated_at.desc()).all()

    def update_dashboard(
        self,
        dashboard_id: int,
        name: Optional[str] = None,
        description: Optional[str] = None,
        layout: Optional[Dict] = None,
        metadata: Optional[Dict] = None,
        is_public: Optional[bool] = None,
    ) -> Optional[Dashboard]:
        """Update dashboard properties."""
        dashboard = self.get_dashboard(dashboard_id)
        if not dashboard:
            return None

        if name is not None:
            dashboard.name = name
        if description is not None:
            dashboard.description = description
        if layout is not None:
            dashboard.layout = layout
        if metadata is not None:
            dashboard.metadata = metadata
        if is_public is not None:
            dashboard.is_public = is_public

        self.db.commit()
        self.db.refresh(dashboard)
        logger.info(f"Updated dashboard {dashboard_id}")
        return dashboard

    def delete_dashboard(self, dashboard_id: int) -> bool:
        """Delete a dashboard and all its panels."""
        dashboard = self.get_dashboard(dashboard_id)
        if not dashboard:
            return False

        self.db.delete(dashboard)
        self.db.commit()
        logger.info(f"Deleted dashboard {dashboard_id}")
        return True

    def duplicate_dashboard(
        self, dashboard_id: int, new_name: str, owner_id: int
    ) -> Optional[Dashboard]:
        """Duplicate a dashboard with all its panels."""
        original = self.get_dashboard(dashboard_id)
        if not original:
            return None

        # Create new dashboard
        new_dashboard = Dashboard(
            project_id=original.project_id,
            name=new_name,
            description=f"Copy of {original.name}",
            owner_id=owner_id,
            layout=original.layout.copy() if original.layout else {},
            metadata=original.metadata.copy() if original.metadata else {},
            is_public=False,
            is_template=False,
        )
        self.db.add(new_dashboard)
        self.db.flush()

        # Duplicate panels
        for panel in original.panels:
            new_panel = Panel(
                dashboard_id=new_dashboard.id,
                panel_key=panel.panel_key,
                title=panel.title,
                description=panel.description,
                viz_type=panel.viz_type,
                data_source=panel.data_source.copy() if panel.data_source else {},
                viz_config=panel.viz_config.copy() if panel.viz_config else {},
                position=panel.position.copy() if panel.position else {},
                auto_refresh=panel.auto_refresh,
                refresh_interval=panel.refresh_interval,
            )
            self.db.add(new_panel)

        self.db.commit()
        self.db.refresh(new_dashboard)
        logger.info(f"Duplicated dashboard {dashboard_id} → {new_dashboard.id}")
        return new_dashboard

    # ==================== Panel Operations ====================

    def add_panel(
        self,
        dashboard_id: int,
        panel_key: str,
        title: str,
        viz_type: str,
        data_source: Dict,
        position: Dict,
        description: str = "",
        viz_config: Optional[Dict] = None,
        auto_refresh: bool = False,
        refresh_interval: Optional[int] = None,
    ) -> Optional[Panel]:
        """Add a panel to a dashboard."""
        # Check dashboard exists
        dashboard = self.get_dashboard(dashboard_id)
        if not dashboard:
            return None

        # Check panel_key uniqueness within dashboard
        existing = (
            self.db.query(Panel)
            .filter(
                and_(Panel.dashboard_id == dashboard_id, Panel.panel_key == panel_key)
            )
            .first()
        )
        if existing:
            logger.warning(
                f"Panel key {panel_key} already exists in dashboard {dashboard_id}"
            )
            return None

        panel = Panel(
            dashboard_id=dashboard_id,
            panel_key=panel_key,
            title=title,
            description=description,
            viz_type=viz_type,
            data_source=data_source,
            viz_config=viz_config or {},
            position=position,
            auto_refresh=auto_refresh,
            refresh_interval=refresh_interval,
        )
        self.db.add(panel)
        self.db.commit()
        self.db.refresh(panel)
        logger.info(f"Added panel {panel.id} to dashboard {dashboard_id}")
        return panel

    def get_panel(self, panel_id: int) -> Optional[Panel]:
        """Get panel by ID."""
        return self.db.query(Panel).filter(Panel.id == panel_id).first()

    def get_panel_by_key(self, dashboard_id: int, panel_key: str) -> Optional[Panel]:
        """Get panel by dashboard ID and key."""
        return (
            self.db.query(Panel)
            .filter(
                and_(Panel.dashboard_id == dashboard_id, Panel.panel_key == panel_key)
            )
            .first()
        )

    def list_panels(self, dashboard_id: int) -> List[Panel]:
        """List all panels in a dashboard."""
        return (
            self.db.query(Panel)
            .filter(Panel.dashboard_id == dashboard_id)
            .order_by(Panel.created_at)
            .all()
        )

    def update_panel(
        self,
        panel_id: int,
        title: Optional[str] = None,
        description: Optional[str] = None,
        viz_type: Optional[str] = None,
        data_source: Optional[Dict] = None,
        viz_config: Optional[Dict] = None,
        position: Optional[Dict] = None,
        auto_refresh: Optional[bool] = None,
        refresh_interval: Optional[int] = None,
    ) -> Optional[Panel]:
        """Update panel properties."""
        panel = self.get_panel(panel_id)
        if not panel:
            return None

        if title is not None:
            panel.title = title
        if description is not None:
            panel.description = description
        if viz_type is not None:
            panel.viz_type = viz_type
        if data_source is not None:
            panel.data_source = data_source
        if viz_config is not None:
            panel.viz_config = viz_config
        if position is not None:
            panel.position = position
        if auto_refresh is not None:
            panel.auto_refresh = auto_refresh
        if refresh_interval is not None:
            panel.refresh_interval = refresh_interval

        self.db.commit()
        self.db.refresh(panel)
        logger.info(f"Updated panel {panel_id}")
        return panel

    def delete_panel(self, panel_id: int) -> bool:
        """Delete a panel."""
        panel = self.get_panel(panel_id)
        if not panel:
            return False

        self.db.delete(panel)
        self.db.commit()
        logger.info(f"Deleted panel {panel_id}")
        return True

    def update_panel_positions(
        self, dashboard_id: int, positions: Dict[str, Dict]
    ) -> bool:
        """
        Batch update panel positions.

        Args:
            dashboard_id: Dashboard ID
            positions: Dict mapping panel_key to position dict {x, y, w, h}
        """
        panels = self.list_panels(dashboard_id)
        updated_count = 0

        for panel in panels:
            if panel.panel_key in positions:
                panel.position = positions[panel.panel_key]
                updated_count += 1

        self.db.commit()
        logger.info(
            f"Updated {updated_count} panel positions in dashboard {dashboard_id}"
        )
        return True

    # ==================== Template Operations ====================

    def create_template_from_dashboard(
        self, dashboard_id: int, template_name: str
    ) -> Optional[Dashboard]:
        """Convert a dashboard to a template."""
        dashboard = self.get_dashboard(dashboard_id)
        if not dashboard:
            return None

        # Create a duplicate as template
        template = self.duplicate_dashboard(
            dashboard_id, template_name, dashboard.owner_id
        )
        if template:
            template.is_template = True
            template.is_public = True
            self.db.commit()
            logger.info(f"Created template {template.id} from dashboard {dashboard_id}")
        return template

    def list_templates(self) -> List[Dashboard]:
        """List all available dashboard templates."""
        return (
            self.db.query(Dashboard)
            .filter(Dashboard.is_template == True)
            .order_by(Dashboard.name)
            .all()
        )

    def create_dashboard_from_template(
        self, template_id: int, name: str, project_id: int, owner_id: int
    ) -> Optional[Dashboard]:
        """Create a new dashboard from a template."""
        template = self.get_dashboard(template_id)
        if not template or not template.is_template:
            return None

        # Duplicate template
        dashboard = self.duplicate_dashboard(template_id, name, owner_id)
        if dashboard:
            dashboard.project_id = project_id
            dashboard.is_template = False
            self.db.commit()
            logger.info(f"Created dashboard {dashboard.id} from template {template_id}")
        return dashboard

    # ==================== Utility Methods ====================

    def get_dashboard_statistics(self, dashboard_id: int) -> Dict:
        """Get statistics for a dashboard."""
        dashboard = self.get_dashboard(dashboard_id)
        if not dashboard:
            return {}

        panels = self.list_panels(dashboard_id)
        viz_types = {}
        for panel in panels:
            viz_types[panel.viz_type] = viz_types.get(panel.viz_type, 0) + 1

        return {
            "dashboard_id": dashboard_id,
            "panel_count": len(panels),
            "viz_types": viz_types,
            "auto_refresh_panels": sum(1 for p in panels if p.auto_refresh),
            "created_at": dashboard.created_at.isoformat(),
            "updated_at": dashboard.updated_at.isoformat(),
        }


def get_visualization_workspace_service(db: Session) -> VisualizationWorkspaceService:
    """Get VisualizationWorkspaceService instance."""
    return VisualizationWorkspaceService(db)
