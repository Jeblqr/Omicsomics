import React, { useState, useEffect } from 'react';
import {
  Box,
  Button,
  Card,
  CardContent,
  Chip,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  FormControl,
  Grid,
  IconButton,
  InputLabel,
  LinearProgress,
  MenuItem,
  Paper,
  Select,
  Tab,
  Tabs,
  TextField,
  Typography,
  Alert,
  Stepper,
  Step,
  StepLabel,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Tooltip,
} from '@mui/material';
import {
  Add as AddIcon,
  PlayArrow as PlayIcon,
  Stop as StopIcon,
  Delete as DeleteIcon,
  Edit as EditIcon,
  Visibility as ViewIcon,
  Code as CodeIcon,
  ExpandMore as ExpandMoreIcon,
  Check as CheckIcon,
  Error as ErrorIcon,
  Schedule as ScheduleIcon,
  Description as DescriptionIcon,
} from '@mui/icons-material';
import { API_URL } from '../config';

interface Script {
  id: number;
  script_key: string;
  name: string;
  description: string;
  user_id: number;
  language: 'python' | 'r' | 'bash';
  script_content: string;
  entry_point?: string;
  parameters_schema?: any;
  requirements: string[];
  timeout: number;
  max_memory?: number;
  visibility: 'private' | 'project' | 'public';
  category?: string;
  tags: string[];
  is_verified: boolean;
  total_executions: number;
  successful_executions: number;
  failed_executions: number;
  avg_duration?: number;
  last_executed_at?: string;
  created_at: string;
  updated_at: string;
}

interface ScriptExecution {
  id: number;
  execution_key: string;
  script_id: number;
  user_id: number;
  parameters: any;
  input_file_ids: number[];
  output_file_ids: number[];
  status: 'pending' | 'running' | 'completed' | 'failed' | 'cancelled';
  description?: string;
  result_data?: any;
  output_text?: string;
  error_text?: string;
  exit_code?: number;
  duration?: number;
  memory_usage?: number;
  started_at?: string;
  completed_at?: string;
  created_at: string;
}

interface Template {
  language: string;
  name: string;
  description: string;
  script_content: string;
  parameters_schema: any;
}

const CustomScriptToolsPage: React.FC = () => {
  const [tab, setTab] = useState(0);
  const [scripts, setScripts] = useState<Script[]>([]);
  const [executions, setExecutions] = useState<ScriptExecution[]>([]);
  const [selectedScript, setSelectedScript] = useState<Script | null>(null);
  const [selectedExecution, setSelectedExecution] = useState<ScriptExecution | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Dialog states
  const [createDialogOpen, setCreateDialogOpen] = useState(false);
  const [executeDialogOpen, setExecuteDialogOpen] = useState(false);
  const [viewDialogOpen, setViewDialogOpen] = useState(false);
  const [executionDetailDialogOpen, setExecutionDetailDialogOpen] = useState(false);

  // Form states
  const [formData, setFormData] = useState({
    name: '',
    description: '',
    language: 'python' as 'python' | 'r' | 'bash',
    script_content: '',
    entry_point: '',
    requirements: '',
    timeout: 300,
    max_memory: 1024,
    visibility: 'private' as 'private' | 'project' | 'public',
    category: '',
    tags: '',
  });

  const [executeFormData, setExecuteFormData] = useState({
    parameters: '{}',
    input_file_ids: '',
    description: '',
  });

  // Filters
  const [filters, setFilters] = useState({
    language: '',
    category: '',
    search: '',
  });

  const [templates, setTemplates] = useState<Template[]>([]);

  useEffect(() => {
    loadScripts();
    loadExecutions();
    loadTemplates();
  }, []);

  useEffect(() => {
    loadScripts();
  }, [filters]);

  // Auto-refresh executions for active ones
  useEffect(() => {
    const hasActiveExecutions = executions.some(
      (exec) => exec.status === 'pending' || exec.status === 'running'
    );

    if (hasActiveExecutions) {
      const interval = setInterval(() => {
        loadExecutions();
      }, 3000);
      return () => clearInterval(interval);
    }
  }, [executions]);

  const loadScripts = async () => {
    try {
      setLoading(true);
      const params = new URLSearchParams();
      if (filters.language) params.append('language', filters.language);
      if (filters.category) params.append('category', filters.category);
      if (filters.search) params.append('search', filters.search);

      const response = await fetch(`${API_URL}/custom-scripts/scripts?${params}`, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });

      if (response.ok) {
        const data = await response.json();
        setScripts(data.items);
      } else {
        setError('Failed to load scripts');
      }
    } catch (err) {
      setError('Error loading scripts');
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  const loadExecutions = async () => {
    try {
      const response = await fetch(`${API_URL}/custom-scripts/executions`, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });

      if (response.ok) {
        const data = await response.json();
        setExecutions(data.items);
      }
    } catch (err) {
      console.error('Error loading executions:', err);
    }
  };

  const loadTemplates = async () => {
    try {
      const response = await fetch(`${API_URL}/custom-scripts/templates`, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });

      if (response.ok) {
        const data = await response.json();
        setTemplates(data);
      }
    } catch (err) {
      console.error('Error loading templates:', err);
    }
  };

  const handleCreateScript = async () => {
    try {
      setLoading(true);
      setError(null);

      const payload = {
        ...formData,
        requirements: formData.requirements.split(',').map((r) => r.trim()).filter(Boolean),
        tags: formData.tags.split(',').map((t) => t.trim()).filter(Boolean),
      };

      const response = await fetch(`${API_URL}/custom-scripts/scripts`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify(payload),
      });

      if (response.ok) {
        setCreateDialogOpen(false);
        resetFormData();
        loadScripts();
      } else {
        const errorData = await response.json();
        setError(errorData.detail || 'Failed to create script');
      }
    } catch (err) {
      setError('Error creating script');
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  const handleExecuteScript = async () => {
    if (!selectedScript) return;

    try {
      setLoading(true);
      setError(null);

      const payload = {
        parameters: JSON.parse(executeFormData.parameters || '{}'),
        input_file_ids: executeFormData.input_file_ids
          .split(',')
          .map((id) => parseInt(id.trim()))
          .filter((id) => !isNaN(id)),
        description: executeFormData.description || undefined,
      };

      const response = await fetch(
        `${API_URL}/custom-scripts/scripts/${selectedScript.id}/execute`,
        {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
            Authorization: `Bearer ${localStorage.getItem('token')}`,
          },
          body: JSON.stringify(payload),
        }
      );

      if (response.ok) {
        setExecuteDialogOpen(false);
        resetExecuteFormData();
        loadExecutions();
        setTab(1); // Switch to executions tab
      } else {
        const errorData = await response.json();
        setError(errorData.detail || 'Failed to execute script');
      }
    } catch (err) {
      setError('Error executing script');
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  const handleCancelExecution = async (executionId: number) => {
    try {
      const response = await fetch(
        `${API_URL}/custom-scripts/executions/${executionId}/cancel`,
        {
          method: 'POST',
          headers: {
            Authorization: `Bearer ${localStorage.getItem('token')}`,
          },
        }
      );

      if (response.ok) {
        loadExecutions();
      }
    } catch (err) {
      console.error('Error cancelling execution:', err);
    }
  };

  const handleDeleteScript = async (scriptId: number) => {
    if (!window.confirm('Are you sure you want to delete this script?')) return;

    try {
      const response = await fetch(`${API_URL}/custom-scripts/scripts/${scriptId}`, {
        method: 'DELETE',
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });

      if (response.ok) {
        loadScripts();
      }
    } catch (err) {
      console.error('Error deleting script:', err);
    }
  };

  const handleViewScript = async (scriptId: number) => {
    try {
      const response = await fetch(`${API_URL}/custom-scripts/scripts/${scriptId}`, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });

      if (response.ok) {
        const script = await response.json();
        setSelectedScript(script);
        setViewDialogOpen(true);
      }
    } catch (err) {
      console.error('Error loading script:', err);
    }
  };

  const handleViewExecution = async (executionId: number) => {
    try {
      const response = await fetch(`${API_URL}/custom-scripts/executions/${executionId}`, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });

      if (response.ok) {
        const execution = await response.json();
        setSelectedExecution(execution);
        setExecutionDetailDialogOpen(true);
      }
    } catch (err) {
      console.error('Error loading execution:', err);
    }
  };

  const handleUseTemplate = (template: Template) => {
    setFormData({
      ...formData,
      name: template.name,
      description: template.description,
      language: template.language as 'python' | 'r' | 'bash',
      script_content: template.script_content,
    });
    setCreateDialogOpen(true);
  };

  const resetFormData = () => {
    setFormData({
      name: '',
      description: '',
      language: 'python',
      script_content: '',
      entry_point: '',
      requirements: '',
      timeout: 300,
      max_memory: 1024,
      visibility: 'private',
      category: '',
      tags: '',
    });
  };

  const resetExecuteFormData = () => {
    setExecuteFormData({
      parameters: '{}',
      input_file_ids: '',
      description: '',
    });
  };

  const getStatusColor = (status: string) => {
    switch (status) {
      case 'completed':
        return 'success';
      case 'failed':
        return 'error';
      case 'running':
        return 'info';
      case 'pending':
        return 'warning';
      case 'cancelled':
        return 'default';
      default:
        return 'default';
    }
  };

  const getStatusIcon = (status: string) => {
    switch (status) {
      case 'completed':
        return <CheckIcon fontSize="small" />;
      case 'failed':
        return <ErrorIcon fontSize="small" />;
      case 'running':
        return <PlayIcon fontSize="small" />;
      case 'pending':
        return <ScheduleIcon fontSize="small" />;
      default:
        return null;
    }
  };

  return (
    <Box sx={{ p: 3 }}>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 3 }}>
        <Typography variant="h4">Custom Script Tools</Typography>
        <Button
          variant="contained"
          startIcon={<AddIcon />}
          onClick={() => setCreateDialogOpen(true)}
        >
          Create Script
        </Button>
      </Box>

      {error && (
        <Alert severity="error" sx={{ mb: 2 }} onClose={() => setError(null)}>
          {error}
        </Alert>
      )}

      <Tabs value={tab} onChange={(_, newValue) => setTab(newValue)} sx={{ mb: 3 }}>
        <Tab label="Scripts" />
        <Tab label="Executions" />
        <Tab label="Templates" />
      </Tabs>

      {/* Scripts Tab */}
      {tab === 0 && (
        <Box>
          {/* Filters */}
          <Paper sx={{ p: 2, mb: 3 }}>
            <Grid container spacing={2}>
              <Grid item xs={12} md={4}>
                <TextField
                  fullWidth
                  label="Search"
                  value={filters.search}
                  onChange={(e) => setFilters({ ...filters, search: e.target.value })}
                />
              </Grid>
              <Grid item xs={12} md={4}>
                <FormControl fullWidth>
                  <InputLabel>Language</InputLabel>
                  <Select
                    value={filters.language}
                    label="Language"
                    onChange={(e) => setFilters({ ...filters, language: e.target.value })}
                  >
                    <MenuItem value="">All</MenuItem>
                    <MenuItem value="python">Python</MenuItem>
                    <MenuItem value="r">R</MenuItem>
                    <MenuItem value="bash">Bash</MenuItem>
                  </Select>
                </FormControl>
              </Grid>
              <Grid item xs={12} md={4}>
                <TextField
                  fullWidth
                  label="Category"
                  value={filters.category}
                  onChange={(e) => setFilters({ ...filters, category: e.target.value })}
                />
              </Grid>
            </Grid>
          </Paper>

          {/* Scripts List */}
          {loading && scripts.length === 0 ? (
            <LinearProgress />
          ) : (
            <Grid container spacing={3}>
              {scripts.map((script) => (
                <Grid item xs={12} md={6} lg={4} key={script.id}>
                  <Card>
                    <CardContent>
                      <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 1 }}>
                        <Typography variant="h6" noWrap>
                          {script.name}
                        </Typography>
                        <Chip label={script.language} size="small" color="primary" />
                      </Box>

                      <Typography
                        variant="body2"
                        color="text.secondary"
                        sx={{
                          mb: 2,
                          overflow: 'hidden',
                          textOverflow: 'ellipsis',
                          display: '-webkit-box',
                          WebkitLineClamp: 2,
                          WebkitBoxOrient: 'vertical',
                        }}
                      >
                        {script.description}
                      </Typography>

                      <Box sx={{ display: 'flex', gap: 0.5, mb: 2, flexWrap: 'wrap' }}>
                        {script.tags.map((tag) => (
                          <Chip key={tag} label={tag} size="small" variant="outlined" />
                        ))}
                      </Box>

                      <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 2 }}>
                        <Typography variant="caption" color="text.secondary">
                          Executions: {script.total_executions} (
                          {script.successful_executions} success)
                        </Typography>
                        {script.is_verified && (
                          <Chip label="Verified" size="small" color="success" />
                        )}
                      </Box>

                      <Box sx={{ display: 'flex', gap: 1 }}>
                        <Button
                          size="small"
                          variant="contained"
                          startIcon={<PlayIcon />}
                          onClick={() => {
                            setSelectedScript(script);
                            setExecuteDialogOpen(true);
                          }}
                        >
                          Execute
                        </Button>
                        <IconButton
                          size="small"
                          onClick={() => handleViewScript(script.id)}
                        >
                          <ViewIcon />
                        </IconButton>
                        <IconButton
                          size="small"
                          color="error"
                          onClick={() => handleDeleteScript(script.id)}
                        >
                          <DeleteIcon />
                        </IconButton>
                      </Box>
                    </CardContent>
                  </Card>
                </Grid>
              ))}
            </Grid>
          )}
        </Box>
      )}

      {/* Executions Tab */}
      {tab === 1 && (
        <Box>
          <TableContainer component={Paper}>
            <Table>
              <TableHead>
                <TableRow>
                  <TableCell>Execution Key</TableCell>
                  <TableCell>Script</TableCell>
                  <TableCell>Status</TableCell>
                  <TableCell>Duration</TableCell>
                  <TableCell>Created</TableCell>
                  <TableCell>Actions</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {executions.map((execution) => {
                  const script = scripts.find((s) => s.id === execution.script_id);
                  return (
                    <TableRow key={execution.id}>
                      <TableCell>{execution.execution_key}</TableCell>
                      <TableCell>{script?.name || execution.script_id}</TableCell>
                      <TableCell>
                        <Chip
                          icon={getStatusIcon(execution.status)}
                          label={execution.status}
                          color={getStatusColor(execution.status) as any}
                          size="small"
                        />
                      </TableCell>
                      <TableCell>
                        {execution.duration ? `${execution.duration.toFixed(2)}s` : '-'}
                      </TableCell>
                      <TableCell>
                        {new Date(execution.created_at).toLocaleString()}
                      </TableCell>
                      <TableCell>
                        <IconButton
                          size="small"
                          onClick={() => handleViewExecution(execution.id)}
                        >
                          <ViewIcon />
                        </IconButton>
                        {(execution.status === 'pending' ||
                          execution.status === 'running') && (
                          <IconButton
                            size="small"
                            color="error"
                            onClick={() => handleCancelExecution(execution.id)}
                          >
                            <StopIcon />
                          </IconButton>
                        )}
                      </TableCell>
                    </TableRow>
                  );
                })}
              </TableBody>
            </Table>
          </TableContainer>
        </Box>
      )}

      {/* Templates Tab */}
      {tab === 2 && (
        <Box>
          <Grid container spacing={3}>
            {templates.map((template, index) => (
              <Grid item xs={12} md={6} key={index}>
                <Card>
                  <CardContent>
                    <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 1 }}>
                      <Typography variant="h6">{template.name}</Typography>
                      <Chip label={template.language} size="small" color="primary" />
                    </Box>
                    <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
                      {template.description}
                    </Typography>
                    <Button
                      variant="outlined"
                      startIcon={<CodeIcon />}
                      onClick={() => handleUseTemplate(template)}
                    >
                      Use Template
                    </Button>
                  </CardContent>
                </Card>
              </Grid>
            ))}
          </Grid>
        </Box>
      )}

      {/* Create Script Dialog */}
      <Dialog
        open={createDialogOpen}
        onClose={() => setCreateDialogOpen(false)}
        maxWidth="md"
        fullWidth
      >
        <DialogTitle>Create Custom Script</DialogTitle>
        <DialogContent>
          <Grid container spacing={2} sx={{ mt: 1 }}>
            <Grid item xs={12}>
              <TextField
                fullWidth
                label="Script Name"
                value={formData.name}
                onChange={(e) => setFormData({ ...formData, name: e.target.value })}
                required
              />
            </Grid>
            <Grid item xs={12}>
              <TextField
                fullWidth
                label="Description"
                multiline
                rows={2}
                value={formData.description}
                onChange={(e) => setFormData({ ...formData, description: e.target.value })}
                required
              />
            </Grid>
            <Grid item xs={12} md={6}>
              <FormControl fullWidth required>
                <InputLabel>Language</InputLabel>
                <Select
                  value={formData.language}
                  label="Language"
                  onChange={(e) =>
                    setFormData({
                      ...formData,
                      language: e.target.value as 'python' | 'r' | 'bash',
                    })
                  }
                >
                  <MenuItem value="python">Python</MenuItem>
                  <MenuItem value="r">R</MenuItem>
                  <MenuItem value="bash">Bash</MenuItem>
                </Select>
              </FormControl>
            </Grid>
            <Grid item xs={12} md={6}>
              <FormControl fullWidth>
                <InputLabel>Visibility</InputLabel>
                <Select
                  value={formData.visibility}
                  label="Visibility"
                  onChange={(e) =>
                    setFormData({
                      ...formData,
                      visibility: e.target.value as 'private' | 'project' | 'public',
                    })
                  }
                >
                  <MenuItem value="private">Private</MenuItem>
                  <MenuItem value="project">Project</MenuItem>
                  <MenuItem value="public">Public</MenuItem>
                </Select>
              </FormControl>
            </Grid>
            <Grid item xs={12}>
              <TextField
                fullWidth
                label="Script Content"
                multiline
                rows={10}
                value={formData.script_content}
                onChange={(e) =>
                  setFormData({ ...formData, script_content: e.target.value })
                }
                required
                sx={{ fontFamily: 'monospace' }}
              />
            </Grid>
            <Grid item xs={12} md={6}>
              <TextField
                fullWidth
                label="Requirements (comma-separated)"
                value={formData.requirements}
                onChange={(e) =>
                  setFormData({ ...formData, requirements: e.target.value })
                }
              />
            </Grid>
            <Grid item xs={12} md={6}>
              <TextField
                fullWidth
                label="Tags (comma-separated)"
                value={formData.tags}
                onChange={(e) => setFormData({ ...formData, tags: e.target.value })}
              />
            </Grid>
            <Grid item xs={12} md={4}>
              <TextField
                fullWidth
                label="Timeout (seconds)"
                type="number"
                value={formData.timeout}
                onChange={(e) =>
                  setFormData({ ...formData, timeout: parseInt(e.target.value) })
                }
              />
            </Grid>
            <Grid item xs={12} md={4}>
              <TextField
                fullWidth
                label="Max Memory (MB)"
                type="number"
                value={formData.max_memory}
                onChange={(e) =>
                  setFormData({ ...formData, max_memory: parseInt(e.target.value) })
                }
              />
            </Grid>
            <Grid item xs={12} md={4}>
              <TextField
                fullWidth
                label="Category"
                value={formData.category}
                onChange={(e) => setFormData({ ...formData, category: e.target.value })}
              />
            </Grid>
          </Grid>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setCreateDialogOpen(false)}>Cancel</Button>
          <Button onClick={handleCreateScript} variant="contained" disabled={loading}>
            Create
          </Button>
        </DialogActions>
      </Dialog>

      {/* Execute Script Dialog */}
      <Dialog
        open={executeDialogOpen}
        onClose={() => setExecuteDialogOpen(false)}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>Execute Script: {selectedScript?.name}</DialogTitle>
        <DialogContent>
          <Grid container spacing={2} sx={{ mt: 1 }}>
            <Grid item xs={12}>
              <TextField
                fullWidth
                label="Parameters (JSON)"
                multiline
                rows={4}
                value={executeFormData.parameters}
                onChange={(e) =>
                  setExecuteFormData({ ...executeFormData, parameters: e.target.value })
                }
                placeholder='{"param1": "value1"}'
                sx={{ fontFamily: 'monospace' }}
              />
            </Grid>
            <Grid item xs={12}>
              <TextField
                fullWidth
                label="Input File IDs (comma-separated)"
                value={executeFormData.input_file_ids}
                onChange={(e) =>
                  setExecuteFormData({
                    ...executeFormData,
                    input_file_ids: e.target.value,
                  })
                }
                placeholder="1, 2, 3"
              />
            </Grid>
            <Grid item xs={12}>
              <TextField
                fullWidth
                label="Description (optional)"
                value={executeFormData.description}
                onChange={(e) =>
                  setExecuteFormData({ ...executeFormData, description: e.target.value })
                }
              />
            </Grid>
          </Grid>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setExecuteDialogOpen(false)}>Cancel</Button>
          <Button onClick={handleExecuteScript} variant="contained" disabled={loading}>
            Execute
          </Button>
        </DialogActions>
      </Dialog>

      {/* View Script Dialog */}
      <Dialog
        open={viewDialogOpen}
        onClose={() => setViewDialogOpen(false)}
        maxWidth="md"
        fullWidth
      >
        <DialogTitle>{selectedScript?.name}</DialogTitle>
        <DialogContent>
          {selectedScript && (
            <Box>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
                {selectedScript.description}
              </Typography>
              <Box sx={{ mb: 2 }}>
                <Chip label={selectedScript.language} size="small" sx={{ mr: 1 }} />
                <Chip label={selectedScript.visibility} size="small" sx={{ mr: 1 }} />
                {selectedScript.is_verified && (
                  <Chip label="Verified" size="small" color="success" />
                )}
              </Box>
              <Typography variant="subtitle2" sx={{ mb: 1 }}>
                Script Content:
              </Typography>
              <Paper sx={{ p: 2, bgcolor: 'grey.100', mb: 2 }}>
                <pre style={{ margin: 0, fontSize: '0.875rem', overflow: 'auto' }}>
                  {selectedScript.script_content}
                </pre>
              </Paper>
              <Typography variant="subtitle2" sx={{ mb: 1 }}>
                Statistics:
              </Typography>
              <Typography variant="body2">
                Total Executions: {selectedScript.total_executions}
                <br />
                Successful: {selectedScript.successful_executions}
                <br />
                Failed: {selectedScript.failed_executions}
                {selectedScript.avg_duration && (
                  <>
                    <br />
                    Average Duration: {selectedScript.avg_duration.toFixed(2)}s
                  </>
                )}
              </Typography>
            </Box>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setViewDialogOpen(false)}>Close</Button>
        </DialogActions>
      </Dialog>

      {/* Execution Detail Dialog */}
      <Dialog
        open={executionDetailDialogOpen}
        onClose={() => setExecutionDetailDialogOpen(false)}
        maxWidth="md"
        fullWidth
      >
        <DialogTitle>Execution Details: {selectedExecution?.execution_key}</DialogTitle>
        <DialogContent>
          {selectedExecution && (
            <Box>
              <Box sx={{ mb: 2 }}>
                <Chip
                  icon={getStatusIcon(selectedExecution.status)}
                  label={selectedExecution.status}
                  color={getStatusColor(selectedExecution.status) as any}
                />
              </Box>

              <Accordion defaultExpanded>
                <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                  <Typography>Parameters</Typography>
                </AccordionSummary>
                <AccordionDetails>
                  <Paper sx={{ p: 2, bgcolor: 'grey.100' }}>
                    <pre style={{ margin: 0, fontSize: '0.875rem' }}>
                      {JSON.stringify(selectedExecution.parameters, null, 2)}
                    </pre>
                  </Paper>
                </AccordionDetails>
              </Accordion>

              {selectedExecution.output_text && (
                <Accordion>
                  <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                    <Typography>Output</Typography>
                  </AccordionSummary>
                  <AccordionDetails>
                    <Paper sx={{ p: 2, bgcolor: 'grey.100' }}>
                      <pre style={{ margin: 0, fontSize: '0.875rem' }}>
                        {selectedExecution.output_text}
                      </pre>
                    </Paper>
                  </AccordionDetails>
                </Accordion>
              )}

              {selectedExecution.error_text && (
                <Accordion defaultExpanded={selectedExecution.status === 'failed'}>
                  <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                    <Typography color="error">Error</Typography>
                  </AccordionSummary>
                  <AccordionDetails>
                    <Paper sx={{ p: 2, bgcolor: 'error.light' }}>
                      <pre style={{ margin: 0, fontSize: '0.875rem' }}>
                        {selectedExecution.error_text}
                      </pre>
                    </Paper>
                  </AccordionDetails>
                </Accordion>
              )}

              <Box sx={{ mt: 2 }}>
                <Typography variant="body2">
                  Duration: {selectedExecution.duration?.toFixed(2) || '-'}s
                  <br />
                  Exit Code: {selectedExecution.exit_code ?? '-'}
                  <br />
                  Started: {selectedExecution.started_at ? new Date(selectedExecution.started_at).toLocaleString() : '-'}
                  <br />
                  Completed: {selectedExecution.completed_at ? new Date(selectedExecution.completed_at).toLocaleString() : '-'}
                </Typography>
              </Box>
            </Box>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setExecutionDetailDialogOpen(false)}>Close</Button>
        </DialogActions>
      </Dialog>
    </Box>
  );
};

export default CustomScriptToolsPage;
