/**
 * Visualization Workspace Page
 * 
 * Dashboard builder for creating multi-panel visualizations.
 * Features:
 * - Drag-and-drop panel layout
 * - Multiple visualization types
 * - Save/load dashboards
 * - Template system
 * - Auto-refresh panels
 */

import React, { useState, useEffect } from 'react';
import {
  Box,
  Typography,
  Button,
  Card,
  CardContent,
  Grid,
  IconButton,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  TextField,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  Chip,
  AppBar,
  Toolbar,
  Menu,
  ListItemIcon,
  ListItemText,
  Divider,
  Alert,
  CircularProgress,
} from '@mui/material';
import {
  Add as AddIcon,
  Dashboard as DashboardIcon,
  Save as SaveIcon,
  FolderOpen as OpenIcon,
  ContentCopy as DuplicateIcon,
  Delete as DeleteIcon,
  Settings as SettingsIcon,
  GridOn as GridIcon,
  Refresh as RefreshIcon,
  Edit as EditIcon,
  Public as PublicIcon,
  Lock as PrivateIcon,
  Template as TemplateIcon,
} from '@mui/icons-material';
import axios from 'axios';

interface Dashboard {
  id: number;
  name: string;
  description: string;
  layout: any;
  metadata: any;
  is_public: boolean;
  is_template: boolean;
  panel_count: number;
  panels?: Panel[];
  created_at: string;
  updated_at: string;
}

interface Panel {
  id: number;
  dashboard_id: number;
  panel_key: string;
  title: string;
  description: string;
  viz_type: string;
  data_source: any;
  viz_config: any;
  position: {
    x: number;
    y: number;
    w: number;
    h: number;
  };
  auto_refresh: boolean;
  refresh_interval: number | null;
  created_at: string;
  updated_at: string;
}

const VIZ_TYPES = [
  { value: 'line', label: 'Line Chart', icon: '📈' },
  { value: 'bar', label: 'Bar Chart', icon: '📊' },
  { value: 'scatter', label: 'Scatter Plot', icon: '⚫' },
  { value: 'heatmap', label: 'Heatmap', icon: '🟥' },
  { value: 'table', label: 'Table', icon: '📋' },
  { value: 'metric', label: 'Metric', icon: '🔢' },
  { value: 'text', label: 'Text', icon: '📝' },
];

const VisualizationWorkspacePage: React.FC = () => {
  // State
  const [dashboards, setDashboards] = useState<Dashboard[]>([]);
  const [currentDashboard, setCurrentDashboard] = useState<Dashboard | null>(null);
  const [panels, setPanels] = useState<Panel[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Dialog states
  const [showCreateDialog, setShowCreateDialog] = useState(false);
  const [showLoadDialog, setShowLoadDialog] = useState(false);
  const [showPanelDialog, setShowPanelDialog] = useState(false);
  const [showSettingsDialog, setShowSettingsDialog] = useState(false);
  const [editingPanel, setEditingPanel] = useState<Panel | null>(null);

  // Menu state
  const [menuAnchor, setMenuAnchor] = useState<null | HTMLElement>(null);

  // Form states
  const [dashboardForm, setDashboardForm] = useState({
    name: '',
    description: '',
    project_id: 1,
    is_public: false,
  });

  const [panelForm, setPanelForm] = useState({
    panel_key: '',
    title: '',
    description: '',
    viz_type: 'line',
    data_source: {},
    viz_config: {},
    position: { x: 0, y: 0, w: 6, h: 4 },
    auto_refresh: false,
    refresh_interval: null as number | null,
  });

  // Load dashboards on mount
  useEffect(() => {
    loadDashboards();
  }, []);

  // Load panels when dashboard changes
  useEffect(() => {
    if (currentDashboard) {
      loadPanels(currentDashboard.id);
    }
  }, [currentDashboard]);

  // API Functions
  const loadDashboards = async () => {
    setLoading(true);
    setError(null);
    try {
      const response = await axios.get('/api/viz-workspace/dashboards');
      setDashboards(response.data);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to load dashboards');
    } finally {
      setLoading(false);
    }
  };

  const loadPanels = async (dashboardId: number) => {
    try {
      const response = await axios.get(`/api/viz-workspace/dashboards/${dashboardId}/panels`);
      setPanels(response.data);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to load panels');
    }
  };

  const createDashboard = async () => {
    try {
      const response = await axios.post('/api/viz-workspace/dashboards', dashboardForm);
      setCurrentDashboard(response.data);
      setPanels([]);
      setShowCreateDialog(false);
      loadDashboards();
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to create dashboard');
    }
  };

  const saveDashboard = async () => {
    if (!currentDashboard) return;

    try {
      await axios.put(`/api/viz-workspace/dashboards/${currentDashboard.id}`, {
        name: currentDashboard.name,
        description: currentDashboard.description,
        layout: currentDashboard.layout,
        is_public: currentDashboard.is_public,
      });
      alert('Dashboard saved successfully!');
      loadDashboards();
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to save dashboard');
    }
  };

  const deleteDashboard = async (dashboardId: number) => {
    if (!confirm('Are you sure you want to delete this dashboard?')) return;

    try {
      await axios.delete(`/api/viz-workspace/dashboards/${dashboardId}`);
      if (currentDashboard?.id === dashboardId) {
        setCurrentDashboard(null);
        setPanels([]);
      }
      loadDashboards();
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to delete dashboard');
    }
  };

  const duplicateDashboard = async (dashboardId: number) => {
    const newName = prompt('Enter name for duplicated dashboard:');
    if (!newName) return;

    try {
      const response = await axios.post(
        `/api/viz-workspace/dashboards/${dashboardId}/duplicate`,
        { new_name: newName }
      );
      setCurrentDashboard(response.data);
      loadDashboards();
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to duplicate dashboard');
    }
  };

  const addPanel = async () => {
    if (!currentDashboard) return;

    try {
      const response = await axios.post(
        `/api/viz-workspace/dashboards/${currentDashboard.id}/panels`,
        panelForm
      );
      setPanels([...panels, response.data]);
      setShowPanelDialog(false);
      resetPanelForm();
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to add panel');
    }
  };

  const updatePanel = async () => {
    if (!editingPanel) return;

    try {
      const response = await axios.put(
        `/api/viz-workspace/panels/${editingPanel.id}`,
        panelForm
      );
      setPanels(panels.map((p) => (p.id === editingPanel.id ? response.data : p)));
      setShowPanelDialog(false);
      setEditingPanel(null);
      resetPanelForm();
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to update panel');
    }
  };

  const deletePanel = async (panelId: number) => {
    if (!confirm('Are you sure you want to delete this panel?')) return;

    try {
      await axios.delete(`/api/viz-workspace/panels/${panelId}`);
      setPanels(panels.filter((p) => p.id !== panelId));
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to delete panel');
    }
  };

  const resetPanelForm = () => {
    setPanelForm({
      panel_key: `panel_${Date.now()}`,
      title: '',
      description: '',
      viz_type: 'line',
      data_source: {},
      viz_config: {},
      position: { x: 0, y: 0, w: 6, h: 4 },
      auto_refresh: false,
      refresh_interval: null,
    });
  };

  const openPanelDialog = (panel?: Panel) => {
    if (panel) {
      setEditingPanel(panel);
      setPanelForm({
        panel_key: panel.panel_key,
        title: panel.title,
        description: panel.description,
        viz_type: panel.viz_type,
        data_source: panel.data_source,
        viz_config: panel.viz_config,
        position: panel.position,
        auto_refresh: panel.auto_refresh,
        refresh_interval: panel.refresh_interval,
      });
    } else {
      setEditingPanel(null);
      resetPanelForm();
    }
    setShowPanelDialog(true);
  };

  return (
    <Box sx={{ height: '100vh', display: 'flex', flexDirection: 'column' }}>
      {/* AppBar */}
      <AppBar position="static" color="default" elevation={1}>
        <Toolbar>
          <DashboardIcon sx={{ mr: 2 }} />
          <Typography variant="h6" sx={{ flexGrow: 1 }}>
            Visualization Workspace
            {currentDashboard && (
              <Typography component="span" variant="subtitle2" sx={{ ml: 2, color: 'text.secondary' }}>
                {currentDashboard.name}
              </Typography>
            )}
          </Typography>

          <Button
            startIcon={<AddIcon />}
            onClick={() => setShowCreateDialog(true)}
            variant="outlined"
            sx={{ mr: 1 }}
          >
            New
          </Button>

          <Button
            startIcon={<OpenIcon />}
            onClick={() => setShowLoadDialog(true)}
            variant="outlined"
            sx={{ mr: 1 }}
          >
            Open
          </Button>

          {currentDashboard && (
            <>
              <Button
                startIcon={<SaveIcon />}
                onClick={saveDashboard}
                variant="contained"
                sx={{ mr: 1 }}
              >
                Save
              </Button>

              <Button
                startIcon={<AddIcon />}
                onClick={() => openPanelDialog()}
                variant="contained"
                color="secondary"
                sx={{ mr: 1 }}
              >
                Add Panel
              </Button>

              <IconButton
                onClick={(e) => setMenuAnchor(e.currentTarget)}
                sx={{ ml: 1 }}
              >
                <SettingsIcon />
              </IconButton>
            </>
          )}
        </Toolbar>
      </AppBar>

      {/* Error Alert */}
      {error && (
        <Alert severity="error" onClose={() => setError(null)} sx={{ m: 2 }}>
          {error}
        </Alert>
      )}

      {/* Main Content */}
      <Box sx={{ flex: 1, overflow: 'auto', p: 3 }}>
        {loading ? (
          <Box display="flex" justifyContent="center" alignItems="center" minHeight="400px">
            <CircularProgress />
          </Box>
        ) : currentDashboard ? (
          <Grid container spacing={2}>
            {panels.length === 0 ? (
              <Grid item xs={12}>
                <Card>
                  <CardContent sx={{ textAlign: 'center', py: 8 }}>
                    <GridIcon sx={{ fontSize: 64, color: 'text.secondary', mb: 2 }} />
                    <Typography variant="h6" gutterBottom>
                      No Panels Yet
                    </Typography>
                    <Typography variant="body2" color="text.secondary" paragraph>
                      Click "Add Panel" to create your first visualization
                    </Typography>
                    <Button
                      variant="contained"
                      startIcon={<AddIcon />}
                      onClick={() => openPanelDialog()}
                    >
                      Add Panel
                    </Button>
                  </CardContent>
                </Card>
              </Grid>
            ) : (
              panels.map((panel) => (
                <Grid item xs={12} md={panel.position.w || 6} key={panel.id}>
                  <Card>
                    <CardContent>
                      <Box display="flex" justifyContent="space-between" alignItems="center" mb={2}>
                        <Box>
                          <Typography variant="h6">{panel.title}</Typography>
                          <Chip label={panel.viz_type} size="small" sx={{ mt: 0.5 }} />
                        </Box>
                        <Box>
                          {panel.auto_refresh && (
                            <IconButton size="small" sx={{ mr: 1 }}>
                              <RefreshIcon fontSize="small" />
                            </IconButton>
                          )}
                          <IconButton size="small" onClick={() => openPanelDialog(panel)}>
                            <EditIcon fontSize="small" />
                          </IconButton>
                          <IconButton size="small" onClick={() => deletePanel(panel.id)}>
                            <DeleteIcon fontSize="small" />
                          </IconButton>
                        </Box>
                      </Box>
                      {panel.description && (
                        <Typography variant="body2" color="text.secondary" paragraph>
                          {panel.description}
                        </Typography>
                      )}
                      <Box
                        sx={{
                          height: 300,
                          bgcolor: 'grey.100',
                          borderRadius: 1,
                          display: 'flex',
                          alignItems: 'center',
                          justifyContent: 'center',
                        }}
                      >
                        <Typography variant="h4" color="text.secondary">
                          {VIZ_TYPES.find((v) => v.value === panel.viz_type)?.icon || '📊'}
                        </Typography>
                      </Box>
                    </CardContent>
                  </Card>
                </Grid>
              ))
            )}
          </Grid>
        ) : (
          <Card>
            <CardContent sx={{ textAlign: 'center', py: 8 }}>
              <DashboardIcon sx={{ fontSize: 64, color: 'text.secondary', mb: 2 }} />
              <Typography variant="h5" gutterBottom>
                Welcome to Visualization Workspace
              </Typography>
              <Typography variant="body1" color="text.secondary" paragraph>
                Create interactive dashboards with multiple visualization panels
              </Typography>
              <Button
                variant="contained"
                startIcon={<AddIcon />}
                onClick={() => setShowCreateDialog(true)}
                sx={{ mr: 2 }}
              >
                Create New Dashboard
              </Button>
              <Button
                variant="outlined"
                startIcon={<OpenIcon />}
                onClick={() => setShowLoadDialog(true)}
              >
                Open Existing
              </Button>
            </CardContent>
          </Card>
        )}
      </Box>

      {/* Create Dashboard Dialog */}
      <Dialog open={showCreateDialog} onClose={() => setShowCreateDialog(false)} maxWidth="sm" fullWidth>
        <DialogTitle>Create New Dashboard</DialogTitle>
        <DialogContent>
          <Box display="flex" flexDirection="column" gap={2} mt={1}>
            <TextField
              label="Dashboard Name"
              value={dashboardForm.name}
              onChange={(e: any) => setDashboardForm({ ...dashboardForm, name: e.target.value })}
              fullWidth
              required
            />
            <TextField
              label="Description"
              value={dashboardForm.description}
              onChange={(e: any) =>
                setDashboardForm({ ...dashboardForm, description: e.target.value })
              }
              fullWidth
              multiline
              rows={3}
            />
          </Box>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setShowCreateDialog(false)}>Cancel</Button>
          <Button onClick={createDashboard} variant="contained" disabled={!dashboardForm.name}>
            Create
          </Button>
        </DialogActions>
      </Dialog>

      {/* Load Dashboard Dialog */}
      <Dialog open={showLoadDialog} onClose={() => setShowLoadDialog(false)} maxWidth="md" fullWidth>
        <DialogTitle>Open Dashboard</DialogTitle>
        <DialogContent>
          <Grid container spacing={2} sx={{ mt: 1 }}>
            {dashboards.map((dashboard) => (
              <Grid item xs={12} sm={6} key={dashboard.id}>
                <Card
                  sx={{
                    cursor: 'pointer',
                    '&:hover': { boxShadow: 3 },
                    border: currentDashboard?.id === dashboard.id ? '2px solid #1976d2' : 'none',
                  }}
                  onClick={() => {
                    setCurrentDashboard(dashboard);
                    setShowLoadDialog(false);
                  }}
                >
                  <CardContent>
                    <Box display="flex" justifyContent="space-between" alignItems="start" mb={1}>
                      <Typography variant="h6">{dashboard.name}</Typography>
                      <Box>
                        {dashboard.is_public && <PublicIcon fontSize="small" color="action" />}
                        {dashboard.is_template && <TemplateIcon fontSize="small" color="primary" />}
                      </Box>
                    </Box>
                    <Typography variant="body2" color="text.secondary" paragraph>
                      {dashboard.description || 'No description'}
                    </Typography>
                    <Chip label={`${dashboard.panel_count} panels`} size="small" />
                    <Box mt={1} display="flex" gap={1}>
                      <IconButton size="small" onClick={(e) => {
                        e.stopPropagation();
                        duplicateDashboard(dashboard.id);
                      }}>
                        <DuplicateIcon fontSize="small" />
                      </IconButton>
                      <IconButton size="small" onClick={(e) => {
                        e.stopPropagation();
                        deleteDashboard(dashboard.id);
                      }}>
                        <DeleteIcon fontSize="small" />
                      </IconButton>
                    </Box>
                  </CardContent>
                </Card>
              </Grid>
            ))}
          </Grid>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setShowLoadDialog(false)}>Close</Button>
        </DialogActions>
      </Dialog>

      {/* Add/Edit Panel Dialog */}
      <Dialog open={showPanelDialog} onClose={() => setShowPanelDialog(false)} maxWidth="md" fullWidth>
        <DialogTitle>{editingPanel ? 'Edit Panel' : 'Add Panel'}</DialogTitle>
        <DialogContent>
          <Box display="flex" flexDirection="column" gap={2} mt={1}>
            <TextField
              label="Panel Title"
              value={panelForm.title}
              onChange={(e: any) => setPanelForm({ ...panelForm, title: e.target.value })}
              fullWidth
              required
            />
            <TextField
              label="Description"
              value={panelForm.description}
              onChange={(e: any) => setPanelForm({ ...panelForm, description: e.target.value })}
              fullWidth
              multiline
              rows={2}
            />
            <FormControl fullWidth>
              <InputLabel>Visualization Type</InputLabel>
              <Select
                value={panelForm.viz_type}
                onChange={(e: any) => setPanelForm({ ...panelForm, viz_type: e.target.value })}
                label="Visualization Type"
              >
                {VIZ_TYPES.map((type) => (
                  <MenuItem key={type.value} value={type.value}>
                    {type.icon} {type.label}
                  </MenuItem>
                ))}
              </Select>
            </FormControl>
            <Grid container spacing={2}>
              <Grid item xs={6}>
                <TextField
                  label="Width (grid units)"
                  type="number"
                  value={panelForm.position.w}
                  onChange={(e: any) =>
                    setPanelForm({
                      ...panelForm,
                      position: { ...panelForm.position, w: parseInt(e.target.value) },
                    })
                  }
                  fullWidth
                  InputProps={{ inputProps: { min: 1, max: 12 } }}
                />
              </Grid>
              <Grid item xs={6}>
                <TextField
                  label="Height (grid units)"
                  type="number"
                  value={panelForm.position.h}
                  onChange={(e: any) =>
                    setPanelForm({
                      ...panelForm,
                      position: { ...panelForm.position, h: parseInt(e.target.value) },
                    })
                  }
                  fullWidth
                  InputProps={{ inputProps: { min: 1, max: 12 } }}
                />
              </Grid>
            </Grid>
          </Box>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setShowPanelDialog(false)}>Cancel</Button>
          <Button
            onClick={editingPanel ? updatePanel : addPanel}
            variant="contained"
            disabled={!panelForm.title}
          >
            {editingPanel ? 'Update' : 'Add'}
          </Button>
        </DialogActions>
      </Dialog>

      {/* Settings Menu */}
      <Menu
        anchorEl={menuAnchor}
        open={Boolean(menuAnchor)}
        onClose={() => setMenuAnchor(null)}
      >
        <MenuItem onClick={() => {
          setMenuAnchor(null);
          setShowSettingsDialog(true);
        }}>
          <ListItemIcon>
            <SettingsIcon fontSize="small" />
          </ListItemIcon>
          <ListItemText>Dashboard Settings</ListItemText>
        </MenuItem>
        <MenuItem onClick={() => {
          setMenuAnchor(null);
          if (currentDashboard) duplicateDashboard(currentDashboard.id);
        }}>
          <ListItemIcon>
            <DuplicateIcon fontSize="small" />
          </ListItemIcon>
          <ListItemText>Duplicate</ListItemText>
        </MenuItem>
        <Divider />
        <MenuItem onClick={() => {
          setMenuAnchor(null);
          if (currentDashboard) deleteDashboard(currentDashboard.id);
        }}>
          <ListItemIcon>
            <DeleteIcon fontSize="small" color="error" />
          </ListItemIcon>
          <ListItemText>Delete Dashboard</ListItemText>
        </MenuItem>
      </Menu>
    </Box>
  );
};

export default VisualizationWorkspacePage;
