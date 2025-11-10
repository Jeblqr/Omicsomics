/**
 * Dataset Input Selector Component
 * 
 * Allows users to select datasets as input sources for pipeline nodes.
 * Used in Pipeline Builder to replace file selection with dataset selection.
 */

import React, { useState, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Box,
  Typography,
  Card,
  CardContent,
  CardActions,
  Grid,
  Chip,
  TextField,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  CircularProgress,
  Alert,
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  Checkbox,
  IconButton,
  Divider,
} from '@mui/material';
import {
  DatasetOutlined,
  FolderOutlined,
  CheckCircleOutline,
  WarningOutlined,
  InfoOutlined,
  RefreshOutlined,
} from '@mui/icons-material';
import axios from 'axios';

interface Dataset {
  id: number;
  name: string;
  description: string;
  data_type: string;
  status: string;
  file_count: number;
  created_at: string;
  tags?: Array<{ name: string; color: string }>;
}

interface DatasetFile {
  id: number;
  path: string;
  name: string;
  type: string;
  role: string;
  size: number;
}

interface DatasetInputSelectorProps {
  open: boolean;
  onClose: () => void;
  onSelect: (datasetId: number, files: DatasetFile[]) => void;
  nodeId: string;
  inputName: string;
  requiredFileTypes?: string[];
  projectId?: number;
}

const DatasetInputSelector: React.FC<DatasetInputSelectorProps> = ({
  open,
  onClose,
  onSelect,
  nodeId,
  inputName,
  requiredFileTypes,
  projectId,
}) => {
  const [datasets, setDatasets] = useState<Dataset[]>([]);
  const [selectedDataset, setSelectedDataset] = useState<Dataset | null>(null);
  const [datasetFiles, setDatasetFiles] = useState<DatasetFile[]>([]);
  const [selectedFiles, setSelectedFiles] = useState<Set<number>>(new Set());
  const [loading, setLoading] = useState(false);
  const [validating, setValidating] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [validationResult, setValidationResult] = useState<{
    is_valid: boolean;
    errors: string[];
  } | null>(null);

  // Filter states
  const [filterDataType, setFilterDataType] = useState<string>('');
  const [searchQuery, setSearchQuery] = useState('');

  // Load datasets
  useEffect(() => {
    if (open) {
      loadDatasets();
    }
  }, [open, projectId]);

  const loadDatasets = async () => {
    setLoading(true);
    setError(null);
    try {
      const params: any = { status: 'active' };
      if (projectId) {
        params.project_id = projectId;
      }
      if (filterDataType) {
        params.data_type = filterDataType;
      }
      if (searchQuery) {
        params.search = searchQuery;
      }

      const response = await axios.get('/api/datasets/', { params });
      setDatasets(response.data);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to load datasets');
    } finally {
      setLoading(false);
    }
  };

  // Load dataset files when dataset is selected
  useEffect(() => {
    if (selectedDataset) {
      loadDatasetFiles(selectedDataset.id);
      validateDataset(selectedDataset.id);
    } else {
      setDatasetFiles([]);
      setValidationResult(null);
    }
  }, [selectedDataset]);

  const loadDatasetFiles = async (datasetId: number) => {
    try {
      const response = await axios.post('/api/pipeline-datasets/input-files', {
        dataset_id: datasetId,
      });
      setDatasetFiles(response.data.files);
      // Auto-select all files initially
      setSelectedFiles(new Set(response.data.files.map((f: DatasetFile) => f.id)));
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to load dataset files');
    }
  };

  const validateDataset = async (datasetId: number) => {
    setValidating(true);
    try {
      const response = await axios.post('/api/pipeline-datasets/validate', {
        dataset_id: datasetId,
        required_file_types: requiredFileTypes,
      });
      setValidationResult(response.data);
    } catch (err: any) {
      console.error('Validation error:', err);
    } finally {
      setValidating(false);
    }
  };

  const handleDatasetSelect = (dataset: Dataset) => {
    setSelectedDataset(dataset);
  };

  const handleFileToggle = (fileId: number) => {
    const newSelected = new Set(selectedFiles);
    if (newSelected.has(fileId)) {
      newSelected.delete(fileId);
    } else {
      newSelected.add(fileId);
    }
    setSelectedFiles(newSelected);
  };

  const handleConfirm = () => {
    if (!selectedDataset) return;

    const filesToPass = datasetFiles.filter((f) => selectedFiles.has(f.id));
    onSelect(selectedDataset.id, filesToPass);
    handleClose();
  };

  const handleClose = () => {
    setSelectedDataset(null);
    setDatasetFiles([]);
    setSelectedFiles(new Set());
    setValidationResult(null);
    setError(null);
    onClose();
  };

  const formatFileSize = (bytes: number): string => {
    if (bytes < 1024) return bytes + ' B';
    if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(2) + ' KB';
    if (bytes < 1024 * 1024 * 1024)
      return (bytes / (1024 * 1024)).toFixed(2) + ' MB';
    return (bytes / (1024 * 1024 * 1024)).toFixed(2) + ' GB';
  };

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="lg" fullWidth>
      <DialogTitle>
        <Box display="flex" alignItems="center" justifyContent="space-between">
          <Box display="flex" alignItems="center" gap={1}>
            <DatasetOutlined />
            <Typography variant="h6">Select Dataset Input</Typography>
          </Box>
          <IconButton onClick={loadDatasets} size="small">
            <RefreshOutlined />
          </IconButton>
        </Box>
        <Typography variant="body2" color="text.secondary">
          Node: {nodeId} • Input: {inputName}
        </Typography>
      </DialogTitle>

      <DialogContent dividers>
        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        {/* Filters */}
        <Box display="flex" gap={2} mb={3}>
          <TextField
            label="Search datasets"
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            size="small"
            sx={{ flex: 1 }}
          />
          <FormControl size="small" sx={{ minWidth: 200 }}>
            <InputLabel>Data Type</InputLabel>
            <Select
              value={filterDataType}
              onChange={(e) => setFilterDataType(e.target.value)}
              label="Data Type"
            >
              <MenuItem value="">All Types</MenuItem>
              <MenuItem value="genomics">Genomics</MenuItem>
              <MenuItem value="transcriptomics">Transcriptomics</MenuItem>
              <MenuItem value="proteomics">Proteomics</MenuItem>
              <MenuItem value="metabolomics">Metabolomics</MenuItem>
              <MenuItem value="single-cell">Single-cell</MenuItem>
              <MenuItem value="general">General</MenuItem>
            </Select>
          </FormControl>
          <Button onClick={loadDatasets} variant="outlined">
            Apply
          </Button>
        </Box>

        <Grid container spacing={3}>
          {/* Dataset List */}
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2" gutterBottom>
              Available Datasets ({datasets.length})
            </Typography>
            <Box sx={{ maxHeight: 400, overflowY: 'auto' }}>
              {loading ? (
                <Box display="flex" justifyContent="center" py={4}>
                  <CircularProgress />
                </Box>
              ) : datasets.length === 0 ? (
                <Alert severity="info">No active datasets found</Alert>
              ) : (
                datasets.map((dataset) => (
                  <Card
                    key={dataset.id}
                    sx={{
                      mb: 2,
                      cursor: 'pointer',
                      border:
                        selectedDataset?.id === dataset.id
                          ? '2px solid #1976d2'
                          : '1px solid #e0e0e0',
                    }}
                    onClick={() => handleDatasetSelect(dataset)}
                  >
                    <CardContent>
                      <Box display="flex" alignItems="center" gap={1} mb={1}>
                        <DatasetOutlined color="primary" />
                        <Typography variant="subtitle1">{dataset.name}</Typography>
                        {selectedDataset?.id === dataset.id && (
                          <CheckCircleOutline color="primary" fontSize="small" />
                        )}
                      </Box>
                      <Typography variant="body2" color="text.secondary" mb={1}>
                        {dataset.description}
                      </Typography>
                      <Box display="flex" gap={1} flexWrap="wrap">
                        <Chip
                          label={dataset.data_type}
                          size="small"
                          color="primary"
                          variant="outlined"
                        />
                        <Chip
                          label={`${dataset.file_count} files`}
                          size="small"
                          icon={<FolderOutlined />}
                        />
                        {dataset.tags?.map((tag) => (
                          <Chip
                            key={tag.name}
                            label={tag.name}
                            size="small"
                            sx={{
                              backgroundColor: tag.color + '20',
                              color: tag.color,
                            }}
                          />
                        ))}
                      </Box>
                    </CardContent>
                  </Card>
                ))
              )}
            </Box>
          </Grid>

          {/* Dataset Files */}
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2" gutterBottom>
              Dataset Files
            </Typography>

            {validating && (
              <Box display="flex" alignItems="center" gap={1} mb={2}>
                <CircularProgress size={16} />
                <Typography variant="body2">Validating...</Typography>
              </Box>
            )}

            {validationResult && (
              <Alert
                severity={validationResult.is_valid ? 'success' : 'warning'}
                sx={{ mb: 2 }}
              >
                {validationResult.is_valid
                  ? 'Dataset is valid for pipeline input'
                  : 'Validation issues: ' + validationResult.errors.join(', ')}
              </Alert>
            )}

            {selectedDataset ? (
              <Box sx={{ maxHeight: 400, overflowY: 'auto' }}>
                <List>
                  {datasetFiles.map((file) => (
                    <ListItem
                      key={file.id}
                      dense
                      button
                      onClick={() => handleFileToggle(file.id)}
                    >
                      <ListItemIcon>
                        <Checkbox
                          edge="start"
                          checked={selectedFiles.has(file.id)}
                          tabIndex={-1}
                          disableRipple
                        />
                      </ListItemIcon>
                      <ListItemText
                        primary={file.name}
                        secondary={
                          <Box display="flex" gap={1} alignItems="center">
                            <Chip label={file.type} size="small" />
                            <Chip label={file.role} size="small" variant="outlined" />
                            <Typography variant="caption">
                              {formatFileSize(file.size)}
                            </Typography>
                          </Box>
                        }
                      />
                    </ListItem>
                  ))}
                </List>
              </Box>
            ) : (
              <Alert severity="info" icon={<InfoOutlined />}>
                Select a dataset to view its files
              </Alert>
            )}
          </Grid>
        </Grid>

        {/* Selection Summary */}
        {selectedDataset && (
          <Box mt={3}>
            <Divider sx={{ mb: 2 }} />
            <Typography variant="body2" color="text.secondary">
              Selected: <strong>{selectedDataset.name}</strong> •{' '}
              {selectedFiles.size} of {datasetFiles.length} files
            </Typography>
          </Box>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={handleClose}>Cancel</Button>
        <Button
          onClick={handleConfirm}
          variant="contained"
          disabled={!selectedDataset || selectedFiles.size === 0}
        >
          Use Selected Files ({selectedFiles.size})
        </Button>
      </DialogActions>
    </Dialog>
  );
};

export default DatasetInputSelector;
