/**
 * Data Editing Page
 * 
 * In-place data editing with:
 * - Operation builder
 * - Live preview
 * - Validation
 * - Undo/redo
 * - Apply/revert
 */

import React, { useState, useEffect } from 'react';
import {
  Box,
  Container,
  Typography,
  Button,
  Paper,
  Stepper,
  Step,
  StepLabel,
  Card,
  CardContent,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  TextField,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  List,
  ListItem,
  ListItemText,
  ListItemSecondaryAction,
  IconButton,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Alert,
  Snackbar,
  Chip,
  Grid,
} from '@mui/material';
import {
  Add as AddIcon,
  Delete as DeleteIcon,
  Edit as EditIcon,
  Preview as PreviewIcon,
  CheckCircle as ApplyIcon,
  Undo as UndoIcon,
  Refresh as RefreshIcon,
  Save as SaveIcon,
} from '@mui/icons-material';
import axios from 'axios';

interface EditSession {
  id: number;
  session_key: string;
  name: string;
  description?: string;
  file_id: number;
  status: string;
  operations: any[];
  current_operation_index: number;
  preview_data?: any;
  preview_summary?: any;
  validation_errors?: any[];
  validation_warnings?: any[];
  is_valid: boolean;
  output_file_id?: number;
  created_at: string;
}

interface Operation {
  type: string;
  parameters: any;
  description?: string;
}

const DataEditingPage: React.FC = () => {
  // State
  const [sessions, setSessions] = useState<EditSession[]>([]);
  const [currentSession, setCurrentSession] = useState<EditSession | null>(null);
  const [activeStep, setActiveStep] = useState(0);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  // Dialogs
  const [showCreateDialog, setShowCreateDialog] = useState(false);
  const [showAddOperationDialog, setShowAddOperationDialog] = useState(false);
  
  // Forms
  const [sessionForm, setSessionForm] = useState({
    name: '',
    description: '',
    file_id: '',
  });
  
  const [operationForm, setOperationForm] = useState({
    operation_type: 'filter_rows',
    parameters: {},
    description: '',
  });
  
  const [snackbar, setSnackbar] = useState({ open: false, message: '', severity: 'success' as 'success' | 'error' });

  const steps = ['Select File', 'Add Operations', 'Preview', 'Apply'];

  // Load sessions
  const loadSessions = async () => {
    setLoading(true);
    try {
      const response = await axios.get('/api/data-editing/sessions');
      setSessions(response.data);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to load sessions');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadSessions();
  }, []);

  // Create session
  const handleCreateSession = async () => {
    if (!sessionForm.name || !sessionForm.file_id) {
      setSnackbar({ open: true, message: 'Please fill required fields', severity: 'error' });
      return;
    }

    try {
      const response = await axios.post('/api/data-editing/sessions', {
        name: sessionForm.name,
        description: sessionForm.description,
        file_id: parseInt(sessionForm.file_id),
      });
      
      setCurrentSession(response.data);
      setActiveStep(1);
      setShowCreateDialog(false);
      setSnackbar({ open: true, message: 'Session created', severity: 'success' });
      loadSessions();
    } catch (err: any) {
      setSnackbar({ open: true, message: err.response?.data?.detail || 'Failed to create session', severity: 'error' });
    }
  };

  // Add operation
  const handleAddOperation = async () => {
    if (!currentSession) return;

    try {
      const response = await axios.post(
        `/api/data-editing/sessions/${currentSession.id}/operations`,
        operationForm
      );
      
      setCurrentSession(response.data);
      setShowAddOperationDialog(false);
      setOperationForm({ operation_type: 'filter_rows', parameters: {}, description: '' });
      setSnackbar({ open: true, message: 'Operation added', severity: 'success' });
    } catch (err: any) {
      setSnackbar({ open: true, message: err.response?.data?.detail || 'Failed to add operation', severity: 'error' });
    }
  };

  // Remove operation
  const handleRemoveOperation = async (index: number) => {
    if (!currentSession) return;

    try {
      const response = await axios.delete(
        `/api/data-editing/sessions/${currentSession.id}/operations/${index}`
      );
      setCurrentSession(response.data);
      setSnackbar({ open: true, message: 'Operation removed', severity: 'success' });
    } catch (err: any) {
      setSnackbar({ open: true, message: 'Failed to remove operation', severity: 'error' });
    }
  };

  // Preview
  const handlePreview = async () => {
    if (!currentSession) return;

    setLoading(true);
    try {
      const response = await axios.post(
        `/api/data-editing/sessions/${currentSession.id}/preview`,
        { sample_size: 100 }
      );
      
      // Reload session to get preview data
      const sessionResponse = await axios.get(`/api/data-editing/sessions/${currentSession.id}`);
      setCurrentSession(sessionResponse.data);
      setActiveStep(2);
      setSnackbar({ open: true, message: 'Preview generated', severity: 'success' });
    } catch (err: any) {
      setSnackbar({ open: true, message: err.response?.data?.detail || 'Failed to preview', severity: 'error' });
    } finally {
      setLoading(false);
    }
  };

  // Apply
  const handleApply = async () => {
    if (!currentSession) return;

    if (!confirm('Apply operations? This will create a new file with the edited data.')) {
      return;
    }

    setLoading(true);
    try {
      const response = await axios.post(`/api/data-editing/sessions/${currentSession.id}/apply`);
      setCurrentSession(response.data);
      setActiveStep(3);
      setSnackbar({ open: true, message: 'Operations applied successfully', severity: 'success' });
      loadSessions();
    } catch (err: any) {
      setSnackbar({ open: true, message: err.response?.data?.detail || 'Failed to apply', severity: 'error' });
    } finally {
      setLoading(false);
    }
  };

  // Delete session
  const handleDeleteSession = async (sessionId: number) => {
    if (!confirm('Delete this edit session?')) return;

    try {
      await axios.delete(`/api/data-editing/sessions/${sessionId}`);
      setSnackbar({ open: true, message: 'Session deleted', severity: 'success' });
      if (currentSession?.id === sessionId) {
        setCurrentSession(null);
        setActiveStep(0);
      }
      loadSessions();
    } catch (err: any) {
      setSnackbar({ open: true, message: 'Failed to delete session', severity: 'error' });
    }
  };

  // Load session
  const handleLoadSession = async (sessionId: number) => {
    try {
      const response = await axios.get(`/api/data-editing/sessions/${sessionId}`);
      setCurrentSession(response.data);
      
      // Determine step based on status
      const session = response.data;
      if (session.status === 'applied') {
        setActiveStep(3);
      } else if (session.preview_data) {
        setActiveStep(2);
      } else if (session.operations.length > 0) {
        setActiveStep(1);
      } else {
        setActiveStep(1);
      }
    } catch (err: any) {
      setSnackbar({ open: true, message: 'Failed to load session', severity: 'error' });
    }
  };

  // Render parameter inputs based on operation type
  const renderParameterInputs = () => {
    const { operation_type } = operationForm;

    switch (operation_type) {
      case 'filter_rows':
        return (
          <>
            <TextField
              fullWidth
              label="Column"
              value={operationForm.parameters.column || ''}
              onChange={(e) => setOperationForm({
                ...operationForm,
                parameters: { ...operationForm.parameters, column: e.target.value }
              })}
              sx={{ mb: 2 }}
            />
            <FormControl fullWidth sx={{ mb: 2 }}>
              <InputLabel>Operator</InputLabel>
              <Select
                value={operationForm.parameters.operator || 'equals'}
                label="Operator"
                onChange={(e) => setOperationForm({
                  ...operationForm,
                  parameters: { ...operationForm.parameters, operator: e.target.value }
                })}
              >
                <MenuItem value="equals">Equals</MenuItem>
                <MenuItem value="not_equals">Not Equals</MenuItem>
                <MenuItem value="greater_than">Greater Than</MenuItem>
                <MenuItem value="less_than">Less Than</MenuItem>
                <MenuItem value="contains">Contains</MenuItem>
              </Select>
            </FormControl>
            <TextField
              fullWidth
              label="Value"
              value={operationForm.parameters.value || ''}
              onChange={(e) => setOperationForm({
                ...operationForm,
                parameters: { ...operationForm.parameters, value: e.target.value }
              })}
            />
          </>
        );

      case 'delete_column':
        return (
          <TextField
            fullWidth
            label="Columns (comma-separated)"
            value={operationForm.parameters.columns?.join(', ') || ''}
            onChange={(e) => setOperationForm({
              ...operationForm,
              parameters: { 
                ...operationForm.parameters, 
                columns: e.target.value.split(',').map(c => c.trim()) 
              }
            })}
            helperText="Enter column names to delete, separated by commas"
          />
        );

      case 'rename_column':
        return (
          <>
            <TextField
              fullWidth
              label="Old Name"
              value={operationForm.parameters.old_name || ''}
              onChange={(e) => setOperationForm({
                ...operationForm,
                parameters: { ...operationForm.parameters, old_name: e.target.value }
              })}
              sx={{ mb: 2 }}
            />
            <TextField
              fullWidth
              label="New Name"
              value={operationForm.parameters.new_name || ''}
              onChange={(e) => setOperationForm({
                ...operationForm,
                parameters: { ...operationForm.parameters, new_name: e.target.value }
              })}
            />
          </>
        );

      case 'fill_missing':
        return (
          <>
            <TextField
              fullWidth
              label="Column"
              value={operationForm.parameters.column || ''}
              onChange={(e) => setOperationForm({
                ...operationForm,
                parameters: { ...operationForm.parameters, column: e.target.value }
              })}
              sx={{ mb: 2 }}
            />
            <FormControl fullWidth>
              <InputLabel>Method</InputLabel>
              <Select
                value={operationForm.parameters.method || 'mean'}
                label="Method"
                onChange={(e) => setOperationForm({
                  ...operationForm,
                  parameters: { ...operationForm.parameters, method: e.target.value }
                })}
              >
                <MenuItem value="value">Specific Value</MenuItem>
                <MenuItem value="mean">Mean</MenuItem>
                <MenuItem value="median">Median</MenuItem>
                <MenuItem value="forward_fill">Forward Fill</MenuItem>
                <MenuItem value="backward_fill">Backward Fill</MenuItem>
              </Select>
            </FormControl>
          </>
        );

      case 'convert_type':
        return (
          <>
            <TextField
              fullWidth
              label="Column"
              value={operationForm.parameters.column || ''}
              onChange={(e) => setOperationForm({
                ...operationForm,
                parameters: { ...operationForm.parameters, column: e.target.value }
              })}
              sx={{ mb: 2 }}
            />
            <FormControl fullWidth>
              <InputLabel>Target Type</InputLabel>
              <Select
                value={operationForm.parameters.target_type || 'string'}
                label="Target Type"
                onChange={(e) => setOperationForm({
                  ...operationForm,
                  parameters: { ...operationForm.parameters, target_type: e.target.value }
                })}
              >
                <MenuItem value="int">Integer</MenuItem>
                <MenuItem value="float">Float</MenuItem>
                <MenuItem value="string">String</MenuItem>
                <MenuItem value="datetime">Datetime</MenuItem>
                <MenuItem value="bool">Boolean</MenuItem>
              </Select>
            </FormControl>
          </>
        );

      default:
        return (
          <Alert severity="info">
            Select an operation type to see parameter inputs
          </Alert>
        );
    }
  };

  return (
    <Container maxWidth="xl">
      <Box sx={{ mt: 4, mb: 4 }}>
        {/* Header */}
        <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 3 }}>
          <Typography variant="h4" component="h1">
            Data Editing
          </Typography>
          <Button
            variant="contained"
            startIcon={<AddIcon />}
            onClick={() => setShowCreateDialog(true)}
          >
            New Edit Session
          </Button>
        </Box>

        <Grid container spacing={3}>
          {/* Session List */}
          <Grid item xs={12} md={3}>
            <Paper sx={{ p: 2 }}>
              <Typography variant="h6" gutterBottom>
                Sessions
              </Typography>
              <List>
                {sessions.map((session) => (
                  <ListItem
                    key={session.id}
                    button
                    selected={currentSession?.id === session.id}
                    onClick={() => handleLoadSession(session.id)}
                  >
                    <ListItemText
                      primary={session.name}
                      secondary={
                        <>
                          <Chip label={session.status} size="small" sx={{ mr: 1 }} />
                          {session.operations.length} ops
                        </>
                      }
                    />
                    <ListItemSecondaryAction>
                      <IconButton
                        edge="end"
                        onClick={(e) => {
                          e.stopPropagation();
                          handleDeleteSession(session.id);
                        }}
                      >
                        <DeleteIcon />
                      </IconButton>
                    </ListItemSecondaryAction>
                  </ListItem>
                ))}
              </List>
            </Paper>
          </Grid>

          {/* Main Content */}
          <Grid item xs={12} md={9}>
            {currentSession ? (
              <Paper sx={{ p: 3 }}>
                {/* Stepper */}
                <Stepper activeStep={activeStep} sx={{ mb: 4 }}>
                  {steps.map((label) => (
                    <Step key={label}>
                      <StepLabel>{label}</StepLabel>
                    </Step>
                  ))}
                </Stepper>

                {/* Session Info */}
                <Card sx={{ mb: 3 }}>
                  <CardContent>
                    <Typography variant="h6">{currentSession.name}</Typography>
                    {currentSession.description && (
                      <Typography color="text.secondary">{currentSession.description}</Typography>
                    )}
                    <Box sx={{ mt: 2 }}>
                      <Chip label={`Status: ${currentSession.status}`} sx={{ mr: 1 }} />
                      <Chip label={`Operations: ${currentSession.operations.length}`} />
                      {currentSession.is_valid && <Chip label="Valid" color="success" sx={{ ml: 1 }} />}
                    </Box>
                  </CardContent>
                </Card>

                {/* Operations List */}
                {activeStep >= 1 && (
                  <Card sx={{ mb: 3 }}>
                    <CardContent>
                      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
                        <Typography variant="h6">Operations</Typography>
                        <Button
                          variant="outlined"
                          startIcon={<AddIcon />}
                          onClick={() => setShowAddOperationDialog(true)}
                          disabled={currentSession.status === 'applied'}
                        >
                          Add Operation
                        </Button>
                      </Box>
                      
                      {currentSession.operations.length === 0 ? (
                        <Alert severity="info">No operations added yet</Alert>
                      ) : (
                        <List>
                          {currentSession.operations.map((op: any, index: number) => (
                            <ListItem key={index}>
                              <ListItemText
                                primary={`${index + 1}. ${op.type}`}
                                secondary={op.description || JSON.stringify(op.parameters)}
                              />
                              <ListItemSecondaryAction>
                                <IconButton
                                  edge="end"
                                  onClick={() => handleRemoveOperation(index)}
                                  disabled={currentSession.status === 'applied'}
                                >
                                  <DeleteIcon />
                                </IconButton>
                              </ListItemSecondaryAction>
                            </ListItem>
                          ))}
                        </List>
                      )}
                    </CardContent>
                  </Card>
                )}

                {/* Preview */}
                {activeStep >= 2 && currentSession.preview_data && (
                  <Card sx={{ mb: 3 }}>
                    <CardContent>
                      <Typography variant="h6" gutterBottom>Preview</Typography>
                      
                      {/* Summary */}
                      {currentSession.preview_summary && (
                        <Box sx={{ mb: 2 }}>
                          <Typography variant="subtitle2">Summary:</Typography>
                          <Typography variant="body2">
                            Original: {currentSession.preview_summary.original_shape[0]} rows × {currentSession.preview_summary.original_shape[1]} cols
                          </Typography>
                          <Typography variant="body2">
                            Result: {currentSession.preview_summary.result_shape[0]} rows × {currentSession.preview_summary.result_shape[1]} cols
                          </Typography>
                        </Box>
                      )}

                      {/* Validation */}
                      {currentSession.validation_errors && currentSession.validation_errors.length > 0 && (
                        <Alert severity="error" sx={{ mb: 2 }}>
                          Validation errors: {currentSession.validation_errors.length}
                        </Alert>
                      )}

                      {/* Data Table */}
                      <TableContainer sx={{ maxHeight: 400 }}>
                        <Table size="small" stickyHeader>
                          <TableHead>
                            <TableRow>
                              {currentSession.preview_data.columns?.map((col: string) => (
                                <TableCell key={col}><strong>{col}</strong></TableCell>
                              ))}
                            </TableRow>
                          </TableHead>
                          <TableBody>
                            {currentSession.preview_data.data?.slice(0, 10).map((row: any, idx: number) => (
                              <TableRow key={idx}>
                                {currentSession.preview_data.columns?.map((col: string) => (
                                  <TableCell key={col}>{String(row[col])}</TableCell>
                                ))}
                              </TableRow>
                            ))}
                          </TableBody>
                        </Table>
                      </TableContainer>
                    </CardContent>
                  </Card>
                )}

                {/* Actions */}
                <Box sx={{ display: 'flex', gap: 2 }}>
                  {activeStep === 1 && currentSession.operations.length > 0 && (
                    <Button
                      variant="contained"
                      startIcon={<PreviewIcon />}
                      onClick={handlePreview}
                      disabled={loading}
                    >
                      Preview
                    </Button>
                  )}
                  
                  {activeStep === 2 && currentSession.is_valid && (
                    <Button
                      variant="contained"
                      color="success"
                      startIcon={<ApplyIcon />}
                      onClick={handleApply}
                      disabled={loading}
                    >
                      Apply Operations
                    </Button>
                  )}

                  {activeStep === 3 && (
                    <Alert severity="success" sx={{ flex: 1 }}>
                      Operations applied successfully! Output file ID: {currentSession.output_file_id}
                    </Alert>
                  )}
                </Box>
              </Paper>
            ) : (
              <Paper sx={{ p: 8, textAlign: 'center' }}>
                <Typography variant="h6" color="text.secondary" gutterBottom>
                  No Session Selected
                </Typography>
                <Typography color="text.secondary" sx={{ mb: 3 }}>
                  Create a new edit session or select an existing one
                </Typography>
                <Button
                  variant="contained"
                  startIcon={<AddIcon />}
                  onClick={() => setShowCreateDialog(true)}
                >
                  Create Edit Session
                </Button>
              </Paper>
            )}
          </Grid>
        </Grid>

        {/* Create Session Dialog */}
        <Dialog open={showCreateDialog} onClose={() => setShowCreateDialog(false)} maxWidth="sm" fullWidth>
          <DialogTitle>Create Edit Session</DialogTitle>
          <DialogContent>
            <Box sx={{ pt: 2 }}>
              <TextField
                fullWidth
                label="Session Name"
                value={sessionForm.name}
                onChange={(e) => setSessionForm({ ...sessionForm, name: e.target.value })}
                sx={{ mb: 2 }}
                required
              />
              <TextField
                fullWidth
                label="Description"
                value={sessionForm.description}
                onChange={(e) => setSessionForm({ ...sessionForm, description: e.target.value })}
                multiline
                rows={2}
                sx={{ mb: 2 }}
              />
              <TextField
                fullWidth
                label="File ID"
                type="number"
                value={sessionForm.file_id}
                onChange={(e) => setSessionForm({ ...sessionForm, file_id: e.target.value })}
                helperText="ID of the file to edit"
                required
              />
            </Box>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setShowCreateDialog(false)}>Cancel</Button>
            <Button onClick={handleCreateSession} variant="contained">Create</Button>
          </DialogActions>
        </Dialog>

        {/* Add Operation Dialog */}
        <Dialog open={showAddOperationDialog} onClose={() => setShowAddOperationDialog(false)} maxWidth="sm" fullWidth>
          <DialogTitle>Add Operation</DialogTitle>
          <DialogContent>
            <Box sx={{ pt: 2 }}>
              <FormControl fullWidth sx={{ mb: 2 }}>
                <InputLabel>Operation Type</InputLabel>
                <Select
                  value={operationForm.operation_type}
                  label="Operation Type"
                  onChange={(e) => setOperationForm({ 
                    ...operationForm, 
                    operation_type: e.target.value,
                    parameters: {} 
                  })}
                >
                  <MenuItem value="filter_rows">Filter Rows</MenuItem>
                  <MenuItem value="delete_column">Delete Column</MenuItem>
                  <MenuItem value="rename_column">Rename Column</MenuItem>
                  <MenuItem value="fill_missing">Fill Missing Values</MenuItem>
                  <MenuItem value="convert_type">Convert Type</MenuItem>
                  <MenuItem value="sort_rows">Sort Rows</MenuItem>
                  <MenuItem value="deduplicate">Remove Duplicates</MenuItem>
                </Select>
              </FormControl>

              {renderParameterInputs()}

              <TextField
                fullWidth
                label="Description (optional)"
                value={operationForm.description}
                onChange={(e) => setOperationForm({ ...operationForm, description: e.target.value })}
                sx={{ mt: 2 }}
              />
            </Box>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setShowAddOperationDialog(false)}>Cancel</Button>
            <Button onClick={handleAddOperation} variant="contained">Add</Button>
          </DialogActions>
        </Dialog>

        {/* Snackbar */}
        <Snackbar
          open={snackbar.open}
          autoHideDuration={6000}
          onClose={() => setSnackbar({ ...snackbar, open: false })}
        >
          <Alert
            onClose={() => setSnackbar({ ...snackbar, open: false })}
            severity={snackbar.severity}
            sx={{ width: '100%' }}
          >
            {snackbar.message}
          </Alert>
        </Snackbar>
      </Box>
    </Container>
  );
};

export default DataEditingPage;
