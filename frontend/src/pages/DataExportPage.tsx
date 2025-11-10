/**
 * Data Export Page
 * 
 * Provides batch export functionality with:
 * - File selection from Data Browser
 * - Format conversion options
 * - Metadata and lineage inclusion
 * - Progress tracking
 * - Download management
 */

import React, { useState, useEffect } from 'react';
import {
  Box,
  Container,
  Typography,
  Button,
  Paper,
  Grid,
  Card,
  CardContent,
  CardActions,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  TextField,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  FormControlLabel,
  Checkbox,
  LinearProgress,
  Chip,
  IconButton,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Alert,
  Snackbar,
} from '@mui/material';
import {
  Add as AddIcon,
  Download as DownloadIcon,
  Delete as DeleteIcon,
  Cancel as CancelIcon,
  Refresh as RefreshIcon,
  CheckCircle as CheckCircleIcon,
  Error as ErrorIcon,
  HourglassEmpty as PendingIcon,
  PlayArrow as ProcessingIcon,
} from '@mui/icons-material';
import axios from 'axios';

interface ExportJob {
  id: number;
  job_key: string;
  name: string;
  description?: string;
  user_id: number;
  project_id?: number;
  export_format: string;
  include_metadata: boolean;
  include_lineage: boolean;
  compress: boolean;
  file_ids: number[];
  dataset_ids?: number[];
  export_options?: any;
  output_path?: string;
  output_size?: number;
  download_url?: string;
  status: 'pending' | 'processing' | 'completed' | 'failed' | 'cancelled';
  progress: number;
  error_message?: string;
  total_files: number;
  processed_files: number;
  failed_files: number;
  created_at: string;
  started_at?: string;
  completed_at?: string;
  expires_at?: string;
}

interface ExportJobForm {
  name: string;
  description: string;
  file_ids: number[];
  export_format: string;
  include_metadata: boolean;
  include_lineage: boolean;
  compress: boolean;
  ttl_hours: number;
}

const DataExportPage: React.FC = () => {
  // State
  const [jobs, setJobs] = useState<ExportJob[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [showCreateDialog, setShowCreateDialog] = useState(false);
  const [snackbar, setSnackbar] = useState({ open: false, message: '', severity: 'success' as 'success' | 'error' });

  // Form state
  const [formData, setFormData] = useState<ExportJobForm>({
    name: '',
    description: '',
    file_ids: [],
    export_format: 'csv',
    include_metadata: true,
    include_lineage: false,
    compress: true,
    ttl_hours: 48,
  });

  // Load export jobs
  const loadJobs = async () => {
    setLoading(true);
    setError(null);
    try {
      const response = await axios.get('/api/data-export/jobs');
      setJobs(response.data);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to load export jobs');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadJobs();
    
    // Auto-refresh every 5 seconds for active jobs
    const interval = setInterval(() => {
      const hasActiveJobs = jobs.some(job => 
        job.status === 'pending' || job.status === 'processing'
      );
      if (hasActiveJobs) {
        loadJobs();
      }
    }, 5000);

    return () => clearInterval(interval);
  }, [jobs]);

  // Create export job
  const handleCreateJob = async () => {
    if (!formData.name || formData.file_ids.length === 0) {
      setSnackbar({
        open: true,
        message: 'Please provide a name and select files',
        severity: 'error'
      });
      return;
    }

    try {
      await axios.post('/api/data-export/jobs', formData);
      setSnackbar({
        open: true,
        message: 'Export job created successfully',
        severity: 'success'
      });
      setShowCreateDialog(false);
      resetForm();
      loadJobs();
    } catch (err: any) {
      setSnackbar({
        open: true,
        message: err.response?.data?.detail || 'Failed to create export job',
        severity: 'error'
      });
    }
  };

  // Cancel export job
  const handleCancelJob = async (jobId: number) => {
    try {
      await axios.post(`/api/data-export/jobs/${jobId}/cancel`);
      setSnackbar({
        open: true,
        message: 'Export job cancelled',
        severity: 'success'
      });
      loadJobs();
    } catch (err: any) {
      setSnackbar({
        open: true,
        message: err.response?.data?.detail || 'Failed to cancel job',
        severity: 'error'
      });
    }
  };

  // Delete export job
  const handleDeleteJob = async (jobId: number) => {
    if (!confirm('Are you sure you want to delete this export job? This will also delete the exported files.')) {
      return;
    }

    try {
      await axios.delete(`/api/data-export/jobs/${jobId}`);
      setSnackbar({
        open: true,
        message: 'Export job deleted',
        severity: 'success'
      });
      loadJobs();
    } catch (err: any) {
      setSnackbar({
        open: true,
        message: err.response?.data?.detail || 'Failed to delete job',
        severity: 'error'
      });
    }
  };

  // Download export
  const handleDownload = (job: ExportJob) => {
    if (job.download_url) {
      window.open(job.download_url, '_blank');
    }
  };

  // Reset form
  const resetForm = () => {
    setFormData({
      name: '',
      description: '',
      file_ids: [],
      export_format: 'csv',
      include_metadata: true,
      include_lineage: false,
      compress: true,
      ttl_hours: 48,
    });
  };

  // Get status icon
  const getStatusIcon = (status: string) => {
    switch (status) {
      case 'pending':
        return <PendingIcon />;
      case 'processing':
        return <ProcessingIcon />;
      case 'completed':
        return <CheckCircleIcon color="success" />;
      case 'failed':
        return <ErrorIcon color="error" />;
      case 'cancelled':
        return <CancelIcon color="disabled" />;
      default:
        return null;
    }
  };

  // Get status color
  const getStatusColor = (status: string): "default" | "primary" | "secondary" | "error" | "info" | "success" | "warning" => {
    switch (status) {
      case 'pending':
        return 'default';
      case 'processing':
        return 'info';
      case 'completed':
        return 'success';
      case 'failed':
        return 'error';
      case 'cancelled':
        return 'default';
      default:
        return 'default';
    }
  };

  // Format file size
  const formatSize = (bytes?: number): string => {
    if (!bytes) return 'N/A';
    const units = ['B', 'KB', 'MB', 'GB'];
    let size = bytes;
    let unitIndex = 0;
    while (size >= 1024 && unitIndex < units.length - 1) {
      size /= 1024;
      unitIndex++;
    }
    return `${size.toFixed(2)} ${units[unitIndex]}`;
  };

  // Format date
  const formatDate = (dateString?: string): string => {
    if (!dateString) return 'N/A';
    return new Date(dateString).toLocaleString();
  };

  return (
    <Container maxWidth="xl">
      <Box sx={{ mt: 4, mb: 4 }}>
        {/* Header */}
        <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 3 }}>
          <Typography variant="h4" component="h1">
            Data Export
          </Typography>
          <Box>
            <Button
              variant="outlined"
              startIcon={<RefreshIcon />}
              onClick={loadJobs}
              sx={{ mr: 2 }}
            >
              Refresh
            </Button>
            <Button
              variant="contained"
              startIcon={<AddIcon />}
              onClick={() => setShowCreateDialog(true)}
            >
              New Export
            </Button>
          </Box>
        </Box>

        {/* Error Alert */}
        {error && (
          <Alert severity="error" sx={{ mb: 3 }} onClose={() => setError(null)}>
            {error}
          </Alert>
        )}

        {/* Export Jobs */}
        {loading && jobs.length === 0 ? (
          <Box sx={{ textAlign: 'center', py: 8 }}>
            <Typography>Loading export jobs...</Typography>
          </Box>
        ) : jobs.length === 0 ? (
          <Paper sx={{ p: 8, textAlign: 'center' }}>
            <Typography variant="h6" color="text.secondary" gutterBottom>
              No Export Jobs
            </Typography>
            <Typography color="text.secondary" sx={{ mb: 3 }}>
              Create your first export job to batch export files with format conversion
            </Typography>
            <Button
              variant="contained"
              startIcon={<AddIcon />}
              onClick={() => setShowCreateDialog(true)}
            >
              Create Export Job
            </Button>
          </Paper>
        ) : (
          <Grid container spacing={3}>
            {jobs.map((job) => (
              <Grid item xs={12} md={6} lg={4} key={job.id}>
                <Card>
                  <CardContent>
                    {/* Job Name and Status */}
                    <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', mb: 2 }}>
                      <Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
                        {job.name}
                      </Typography>
                      <Chip
                        label={job.status}
                        color={getStatusColor(job.status)}
                        size="small"
                        icon={getStatusIcon(job.status)}
                      />
                    </Box>

                    {/* Description */}
                    {job.description && (
                      <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
                        {job.description}
                      </Typography>
                    )}

                    {/* Details */}
                    <Box sx={{ mb: 2 }}>
                      <Typography variant="body2">
                        <strong>Format:</strong> {job.export_format.toUpperCase()}
                      </Typography>
                      <Typography variant="body2">
                        <strong>Files:</strong> {job.processed_files}/{job.total_files}
                        {job.failed_files > 0 && ` (${job.failed_files} failed)`}
                      </Typography>
                      {job.output_size && (
                        <Typography variant="body2">
                          <strong>Size:</strong> {formatSize(job.output_size)}
                        </Typography>
                      )}
                      <Typography variant="body2">
                        <strong>Created:</strong> {formatDate(job.created_at)}
                      </Typography>
                      {job.expires_at && (
                        <Typography variant="body2" color="warning.main">
                          <strong>Expires:</strong> {formatDate(job.expires_at)}
                        </Typography>
                      )}
                    </Box>

                    {/* Progress Bar */}
                    {(job.status === 'pending' || job.status === 'processing') && (
                      <Box sx={{ mb: 2 }}>
                        <LinearProgress variant="determinate" value={job.progress} />
                        <Typography variant="caption" color="text.secondary">
                          {job.progress}% complete
                        </Typography>
                      </Box>
                    )}

                    {/* Error Message */}
                    {job.error_message && (
                      <Alert severity="error" sx={{ mb: 2 }}>
                        {job.error_message}
                      </Alert>
                    )}
                  </CardContent>

                  <CardActions>
                    {job.status === 'completed' && job.download_url && (
                      <Button
                        size="small"
                        startIcon={<DownloadIcon />}
                        onClick={() => handleDownload(job)}
                      >
                        Download
                      </Button>
                    )}
                    {(job.status === 'pending' || job.status === 'processing') && (
                      <Button
                        size="small"
                        startIcon={<CancelIcon />}
                        onClick={() => handleCancelJob(job.id)}
                      >
                        Cancel
                      </Button>
                    )}
                    <Button
                      size="small"
                      color="error"
                      startIcon={<DeleteIcon />}
                      onClick={() => handleDeleteJob(job.id)}
                    >
                      Delete
                    </Button>
                  </CardActions>
                </Card>
              </Grid>
            ))}
          </Grid>
        )}

        {/* Create Export Dialog */}
        <Dialog
          open={showCreateDialog}
          onClose={() => setShowCreateDialog(false)}
          maxWidth="md"
          fullWidth
        >
          <DialogTitle>Create Export Job</DialogTitle>
          <DialogContent>
            <Box sx={{ pt: 2 }}>
              <TextField
                fullWidth
                label="Export Name"
                value={formData.name}
                onChange={(e) => setFormData({ ...formData, name: e.target.value })}
                sx={{ mb: 2 }}
                required
              />

              <TextField
                fullWidth
                label="Description"
                value={formData.description}
                onChange={(e) => setFormData({ ...formData, description: e.target.value })}
                multiline
                rows={2}
                sx={{ mb: 2 }}
              />

              <TextField
                fullWidth
                label="File IDs (comma-separated)"
                value={formData.file_ids.join(', ')}
                onChange={(e) =>
                  setFormData({
                    ...formData,
                    file_ids: e.target.value.split(',').map(id => parseInt(id.trim())).filter(id => !isNaN(id))
                  })
                }
                helperText="Enter file IDs to export, separated by commas"
                sx={{ mb: 2 }}
                required
              />

              <FormControl fullWidth sx={{ mb: 2 }}>
                <InputLabel>Export Format</InputLabel>
                <Select
                  value={formData.export_format}
                  label="Export Format"
                  onChange={(e) => setFormData({ ...formData, export_format: e.target.value })}
                >
                  <MenuItem value="csv">CSV</MenuItem>
                  <MenuItem value="tsv">TSV</MenuItem>
                  <MenuItem value="json">JSON</MenuItem>
                  <MenuItem value="excel">Excel</MenuItem>
                  <MenuItem value="parquet">Parquet</MenuItem>
                  <MenuItem value="hdf5">HDF5</MenuItem>
                  <MenuItem value="zip">ZIP</MenuItem>
                </Select>
              </FormControl>

              <TextField
                fullWidth
                type="number"
                label="Time to Live (hours)"
                value={formData.ttl_hours}
                onChange={(e) => setFormData({ ...formData, ttl_hours: parseInt(e.target.value) })}
                helperText="Export will be automatically deleted after this time"
                inputProps={{ min: 1, max: 168 }}
                sx={{ mb: 2 }}
              />

              <FormControlLabel
                control={
                  <Checkbox
                    checked={formData.include_metadata}
                    onChange={(e) => setFormData({ ...formData, include_metadata: e.target.checked })}
                  />
                }
                label="Include Metadata"
              />

              <FormControlLabel
                control={
                  <Checkbox
                    checked={formData.include_lineage}
                    onChange={(e) => setFormData({ ...formData, include_lineage: e.target.checked })}
                  />
                }
                label="Include Lineage"
              />

              <FormControlLabel
                control={
                  <Checkbox
                    checked={formData.compress}
                    onChange={(e) => setFormData({ ...formData, compress: e.target.checked })}
                  />
                }
                label="Compress (ZIP)"
              />
            </Box>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setShowCreateDialog(false)}>Cancel</Button>
            <Button onClick={handleCreateJob} variant="contained">
              Create Export
            </Button>
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

export default DataExportPage;
