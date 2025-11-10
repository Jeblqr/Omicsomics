import React, { useState, useEffect } from 'react';
import {
  Box,
  Button,
  Card,
  CardContent,
  Typography,
  Grid,
  Chip,
  TextField,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  IconButton,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
  Tabs,
  Tab,
} from '@mui/material';
import {
  Add as AddIcon,
  Edit as EditIcon,
  Delete as DeleteIcon,
  Archive as ArchiveIcon,
  Visibility as ViewIcon,
  Timeline as TimelineIcon,
  Label as LabelIcon,
  FolderOpen as FolderIcon,
} from '@mui/icons-material';
import axios from 'axios';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:8000';

interface Dataset {
  id: string;
  name: string;
  description?: string;
  version: number;
  data_type?: string;
  file_format?: string;
  status: string;
  created_at: string;
  updated_at: string;
  file_count: number;
  tag_names: string[];
}

interface DatasetFile {
  id: string;
  file_name: string;
  file_path: string;
  file_type?: string;
  file_size?: number;
  role?: string;
  description?: string;
  added_at: string;
}

interface LineageRecord {
  id: string;
  operation_type: string;
  tool_name?: string;
  tool_version?: string;
  executed_at: string;
}

const DatasetManagerPage: React.FC = () => {
  const [datasets, setDatasets] = useState<Dataset[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  // Filters
  const [searchQuery, setSearchQuery] = useState('');
  const [filterStatus, setFilterStatus] = useState<string>('all');
  const [filterDataType, setFilterDataType] = useState<string>('all');
  
  // Modals
  const [showCreateDialog, setShowCreateDialog] = useState(false);
  const [showDetailDialog, setShowDetailDialog] = useState(false);
  const [selectedDataset, setSelectedDataset] = useState<Dataset | null>(null);
  const [selectedTab, setSelectedTab] = useState(0);
  
  // Dataset details
  const [datasetFiles, setDatasetFiles] = useState<DatasetFile[]>([]);
  const [datasetLineage, setDatasetLineage] = useState<LineageRecord[]>([]);
  const [datasetVersions, setDatasetVersions] = useState<Dataset[]>([]);
  
  // Create form
  const [newDataset, setNewDataset] = useState({
    name: '',
    description: '',
    data_type: '',
    project_id: 'default-project', // TODO: Get from context
  });

  useEffect(() => {
    loadDatasets();
  }, [filterStatus, filterDataType, searchQuery]);

  const loadDatasets = async () => {
    setLoading(true);
    setError(null);
    
    try {
      const params: any = {};
      if (filterStatus !== 'all') params.status = filterStatus;
      if (filterDataType !== 'all') params.data_type = filterDataType;
      if (searchQuery) params.search = searchQuery;
      
      const response = await axios.get(`${API_BASE_URL}/api/datasets/`, { params });
      setDatasets(response.data);
    } catch (err: any) {
      console.error('Failed to load datasets:', err);
      setError(err.message || 'Failed to load datasets');
    } finally {
      setLoading(false);
    }
  };

  const handleCreateDataset = async () => {
    try {
      await axios.post(`${API_BASE_URL}/api/datasets/`, newDataset);
      setShowCreateDialog(false);
      setNewDataset({
        name: '',
        description: '',
        data_type: '',
        project_id: 'default-project',
      });
      loadDatasets();
    } catch (err: any) {
      console.error('Failed to create dataset:', err);
      alert('Failed to create dataset: ' + (err.response?.data?.detail || err.message));
    }
  };

  const handleViewDataset = async (dataset: Dataset) => {
    setSelectedDataset(dataset);
    setShowDetailDialog(true);
    setSelectedTab(0);
    
    // Load dataset details
    try {
      const [filesRes, lineageRes, versionsRes] = await Promise.all([
        axios.get(`${API_BASE_URL}/api/datasets/${dataset.id}/files`),
        axios.get(`${API_BASE_URL}/api/datasets/${dataset.id}/lineage`),
        axios.get(`${API_BASE_URL}/api/datasets/${dataset.id}/versions`),
      ]);
      
      setDatasetFiles(filesRes.data);
      setDatasetLineage(lineageRes.data);
      setDatasetVersions(versionsRes.data);
    } catch (err) {
      console.error('Failed to load dataset details:', err);
    }
  };

  const handleArchiveDataset = async (datasetId: string) => {
    if (!confirm('Archive this dataset?')) return;
    
    try {
      await axios.post(`${API_BASE_URL}/api/datasets/${datasetId}/archive`);
      loadDatasets();
    } catch (err: any) {
      alert('Failed to archive dataset: ' + (err.response?.data?.detail || err.message));
    }
  };

  const handleDeleteDataset = async (datasetId: string) => {
    if (!confirm('Delete this dataset? This action cannot be undone.')) return;
    
    try {
      await axios.delete(`${API_BASE_URL}/api/datasets/${datasetId}`);
      loadDatasets();
    } catch (err: any) {
      alert('Failed to delete dataset: ' + (err.response?.data?.detail || err.message));
    }
  };

  const formatFileSize = (bytes?: number) => {
    if (!bytes) return 'N/A';
    if (bytes < 1024) return `${bytes} B`;
    if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
    if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
    return `${(bytes / (1024 * 1024 * 1024)).toFixed(1)} GB`;
  };

  const formatDate = (dateString: string) => {
    return new Date(dateString).toLocaleString();
  };

  return (
    <Box sx={{ p: 3, bgcolor: '#0f172a', minHeight: '100vh' }}>
      {/* Header */}
      <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
        <Box>
          <Typography variant="h4" sx={{ color: '#f3f4f6', fontWeight: 700 }}>
            📦 Dataset Manager
          </Typography>
          <Typography variant="body2" sx={{ color: '#9ca3af', mt: 0.5 }}>
            Organize and manage your data with version control and lineage tracking
          </Typography>
        </Box>
        <Button
          variant="contained"
          startIcon={<AddIcon />}
          onClick={() => setShowCreateDialog(true)}
        >
          Create Dataset
        </Button>
      </Box>

      {/* Filters */}
      <Grid container spacing={2} sx={{ mb: 3 }}>
        <Grid item xs={12} md={6}>
          <TextField
            fullWidth
            placeholder="Search datasets..."
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            sx={{ bgcolor: '#1f2937' }}
          />
        </Grid>
        <Grid item xs={12} md={3}>
          <FormControl fullWidth>
            <InputLabel>Status</InputLabel>
            <Select
              value={filterStatus}
              onChange={(e) => setFilterStatus(e.target.value)}
              sx={{ bgcolor: '#1f2937' }}
            >
              <MenuItem value="all">All</MenuItem>
              <MenuItem value="active">Active</MenuItem>
              <MenuItem value="archived">Archived</MenuItem>
            </Select>
          </FormControl>
        </Grid>
        <Grid item xs={12} md={3}>
          <FormControl fullWidth>
            <InputLabel>Data Type</InputLabel>
            <Select
              value={filterDataType}
              onChange={(e) => setFilterDataType(e.target.value)}
              sx={{ bgcolor: '#1f2937' }}
            >
              <MenuItem value="all">All</MenuItem>
              <MenuItem value="genomics">Genomics</MenuItem>
              <MenuItem value="proteomics">Proteomics</MenuItem>
              <MenuItem value="transcriptomics">Transcriptomics</MenuItem>
              <MenuItem value="metabolomics">Metabolomics</MenuItem>
            </Select>
          </FormControl>
        </Grid>
      </Grid>

      {/* Error message */}
      {error && (
        <Typography color="error" sx={{ mb: 2 }}>
          {error}
        </Typography>
      )}

      {/* Dataset Grid */}
      <Grid container spacing={2}>
        {loading ? (
          <Grid item xs={12}>
            <Typography sx={{ color: '#9ca3af', textAlign: 'center', py: 4 }}>
              Loading datasets...
            </Typography>
          </Grid>
        ) : datasets.length === 0 ? (
          <Grid item xs={12}>
            <Typography sx={{ color: '#9ca3af', textAlign: 'center', py: 4 }}>
              No datasets found. Create one to get started!
            </Typography>
          </Grid>
        ) : (
          datasets.map((dataset) => (
            <Grid item xs={12} md={6} lg={4} key={dataset.id}>
              <Card sx={{ bgcolor: '#1f2937', border: '1px solid #374151' }}>
                <CardContent>
                  <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 1 }}>
                    <Typography variant="h6" sx={{ color: '#f3f4f6' }}>
                      {dataset.name}
                    </Typography>
                    <Chip
                      label={`v${dataset.version}`}
                      size="small"
                      color="primary"
                    />
                  </Box>
                  
                  <Typography variant="body2" sx={{ color: '#9ca3af', mb: 2, minHeight: 40 }}>
                    {dataset.description || 'No description'}
                  </Typography>
                  
                  <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap', mb: 2 }}>
                    {dataset.data_type && (
                      <Chip label={dataset.data_type} size="small" />
                    )}
                    <Chip
                      label={`${dataset.file_count} files`}
                      size="small"
                      icon={<FolderIcon />}
                    />
                    {dataset.tag_names.map((tag) => (
                      <Chip
                        key={tag}
                        label={tag}
                        size="small"
                        color="secondary"
                      />
                    ))}
                  </Box>
                  
                  <Typography variant="caption" sx={{ color: '#6b7280', display: 'block', mb: 2 }}>
                    Updated: {formatDate(dataset.updated_at)}
                  </Typography>
                  
                  <Box sx={{ display: 'flex', gap: 1 }}>
                    <IconButton
                      size="small"
                      onClick={() => handleViewDataset(dataset)}
                      sx={{ color: '#3b82f6' }}
                    >
                      <ViewIcon />
                    </IconButton>
                    <IconButton
                      size="small"
                      onClick={() => handleArchiveDataset(dataset.id)}
                      sx={{ color: '#f59e0b' }}
                    >
                      <ArchiveIcon />
                    </IconButton>
                    <IconButton
                      size="small"
                      onClick={() => handleDeleteDataset(dataset.id)}
                      sx={{ color: '#ef4444' }}
                    >
                      <DeleteIcon />
                    </IconButton>
                  </Box>
                </CardContent>
              </Card>
            </Grid>
          ))
        )}
      </Grid>

      {/* Create Dataset Dialog */}
      <Dialog
        open={showCreateDialog}
        onClose={() => setShowCreateDialog(false)}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>Create New Dataset</DialogTitle>
        <DialogContent>
          <TextField
            fullWidth
            label="Name"
            value={newDataset.name}
            onChange={(e) => setNewDataset({ ...newDataset, name: e.target.value })}
            margin="normal"
            required
          />
          <TextField
            fullWidth
            label="Description"
            value={newDataset.description}
            onChange={(e) => setNewDataset({ ...newDataset, description: e.target.value })}
            margin="normal"
            multiline
            rows={3}
          />
          <FormControl fullWidth margin="normal">
            <InputLabel>Data Type</InputLabel>
            <Select
              value={newDataset.data_type}
              onChange={(e) => setNewDataset({ ...newDataset, data_type: e.target.value })}
            >
              <MenuItem value="">None</MenuItem>
              <MenuItem value="genomics">Genomics</MenuItem>
              <MenuItem value="proteomics">Proteomics</MenuItem>
              <MenuItem value="transcriptomics">Transcriptomics</MenuItem>
              <MenuItem value="metabolomics">Metabolomics</MenuItem>
              <MenuItem value="other">Other</MenuItem>
            </Select>
          </FormControl>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setShowCreateDialog(false)}>Cancel</Button>
          <Button
            variant="contained"
            onClick={handleCreateDataset}
            disabled={!newDataset.name}
          >
            Create
          </Button>
        </DialogActions>
      </Dialog>

      {/* Dataset Detail Dialog */}
      <Dialog
        open={showDetailDialog}
        onClose={() => setShowDetailDialog(false)}
        maxWidth="lg"
        fullWidth
      >
        {selectedDataset && (
          <>
            <DialogTitle>
              <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                <Typography variant="h6">{selectedDataset.name}</Typography>
                <Chip label={`v${selectedDataset.version}`} color="primary" />
              </Box>
              <Typography variant="body2" color="text.secondary">
                {selectedDataset.description}
              </Typography>
            </DialogTitle>
            <DialogContent>
              <Tabs value={selectedTab} onChange={(_, v) => setSelectedTab(v)}>
                <Tab label="Files" />
                <Tab label="Lineage" />
                <Tab label="Versions" />
              </Tabs>
              
              {selectedTab === 0 && (
                <TableContainer component={Paper} sx={{ mt: 2 }}>
                  <Table>
                    <TableHead>
                      <TableRow>
                        <TableCell>File Name</TableCell>
                        <TableCell>Type</TableCell>
                        <TableCell>Size</TableCell>
                        <TableCell>Role</TableCell>
                        <TableCell>Added</TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {datasetFiles.map((file) => (
                        <TableRow key={file.id}>
                          <TableCell>{file.file_name}</TableCell>
                          <TableCell>{file.file_type || 'N/A'}</TableCell>
                          <TableCell>{formatFileSize(file.file_size)}</TableCell>
                          <TableCell>{file.role || 'N/A'}</TableCell>
                          <TableCell>{formatDate(file.added_at)}</TableCell>
                        </TableRow>
                      ))}
                      {datasetFiles.length === 0 && (
                        <TableRow>
                          <TableCell colSpan={5} align="center">
                            No files in this dataset
                          </TableCell>
                        </TableRow>
                      )}
                    </TableBody>
                  </Table>
                </TableContainer>
              )}
              
              {selectedTab === 1 && (
                <TableContainer component={Paper} sx={{ mt: 2 }}>
                  <Table>
                    <TableHead>
                      <TableRow>
                        <TableCell>Operation</TableCell>
                        <TableCell>Tool</TableCell>
                        <TableCell>Version</TableCell>
                        <TableCell>Executed</TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {datasetLineage.map((record) => (
                        <TableRow key={record.id}>
                          <TableCell>{record.operation_type}</TableCell>
                          <TableCell>{record.tool_name || 'N/A'}</TableCell>
                          <TableCell>{record.tool_version || 'N/A'}</TableCell>
                          <TableCell>{formatDate(record.executed_at)}</TableCell>
                        </TableRow>
                      ))}
                      {datasetLineage.length === 0 && (
                        <TableRow>
                          <TableCell colSpan={4} align="center">
                            No lineage records
                          </TableCell>
                        </TableRow>
                      )}
                    </TableBody>
                  </Table>
                </TableContainer>
              )}
              
              {selectedTab === 2 && (
                <TableContainer component={Paper} sx={{ mt: 2 }}>
                  <Table>
                    <TableHead>
                      <TableRow>
                        <TableCell>Version</TableCell>
                        <TableCell>Description</TableCell>
                        <TableCell>Files</TableCell>
                        <TableCell>Status</TableCell>
                        <TableCell>Created</TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {datasetVersions.map((version) => (
                        <TableRow key={version.id}>
                          <TableCell>v{version.version}</TableCell>
                          <TableCell>{version.description || 'N/A'}</TableCell>
                          <TableCell>{version.file_count}</TableCell>
                          <TableCell>
                            <Chip label={version.status} size="small" />
                          </TableCell>
                          <TableCell>{formatDate(version.created_at)}</TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                </TableContainer>
              )}
            </DialogContent>
            <DialogActions>
              <Button onClick={() => setShowDetailDialog(false)}>Close</Button>
            </DialogActions>
          </>
        )}
      </Dialog>
    </Box>
  );
};

export default DatasetManagerPage;
