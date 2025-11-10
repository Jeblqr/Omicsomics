import React, { useState, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  TextField,
  Checkbox,
  FormControlLabel,
  Typography,
  Box,
  Paper,
  Table,
  TableHead,
  TableRow,
  TableCell,
  TableBody,
  Chip,
  LinearProgress,
  Alert,
  IconButton,
  List,
  ListItem,
  ListItemText,
  ListItemSecondaryAction,
} from '@mui/material';
import DeleteIcon from '@mui/icons-material/Delete';
import PreviewIcon from '@mui/icons-material/Preview';
import MergeIcon from '@mui/icons-material/MergeType';
import VerticalAlignCenterIcon from '@mui/icons-material/VerticalAlignCenter';
import CompareArrowsIcon from '@mui/icons-material/CompareArrows';
import AutoFixHighIcon from '@mui/icons-material/AutoFixHigh';
import axios from 'axios';

interface DataMergerModalProps {
  open: boolean;
  onClose: () => void;
  selectedFiles?: string[];
  onSuccess?: () => void;
}

interface MergeMode {
  value: string;
  label: string;
  description: string;
  icon: React.ReactNode;
}

interface JoinType {
  value: string;
  label: string;
  description: string;
}

interface PreviewData {
  columns: string[];
  rows: any[][];
}

const DataMergerModal: React.FC<DataMergerModalProps> = ({
  open,
  onClose,
  selectedFiles = [],
  onSuccess,
}) => {
  const [files, setFiles] = useState<string[]>(selectedFiles);
  const [mode, setMode] = useState<string>('smart');
  const [joinType, setJoinType] = useState<string>('outer');
  const [keyColumns, setKeyColumns] = useState<string>('');
  const [duplicateHandling, setDuplicateHandling] = useState<string>('keep_all');
  const [addSourceColumn, setAddSourceColumn] = useState<boolean>(true);
  const [ignoreIndex, setIgnoreIndex] = useState<boolean>(false);
  const [sortKeys, setSortKeys] = useState<boolean>(false);
  const [outputFilename, setOutputFilename] = useState<string>('merged_data.csv');
  
  const [mergeModes, setMergeModes] = useState<MergeMode[]>([]);
  const [joinTypes, setJoinTypes] = useState<JoinType[]>([]);
  
  const [previewData, setPreviewData] = useState<PreviewData | null>(null);
  const [previewLoading, setPreviewLoading] = useState(false);
  const [previewError, setPreviewError] = useState<string>('');
  const [detectedMode, setDetectedMode] = useState<string>('');
  const [estimatedDimensions, setEstimatedDimensions] = useState<{rows: number; columns: number}>({rows: 0, columns: 0});
  
  const [merging, setMerging] = useState(false);
  const [jobId, setJobId] = useState<string>('');
  const [progress, setProgress] = useState(0);
  const [error, setError] = useState<string>('');
  const [success, setSuccess] = useState(false);

  // Load merge modes and join types
  useEffect(() => {
    if (open) {
      loadMergeModes();
      loadJoinTypes();
    }
  }, [open]);

  const loadMergeModes = async () => {
    try {
      const response = await axios.get('/api/data-merger/merge-modes');
      const modes = response.data.modes.map((m: any) => ({
        ...m,
        icon: getModeIcon(m.value),
      }));
      setMergeModes(modes);
    } catch (err) {
      console.error('Failed to load merge modes:', err);
    }
  };

  const loadJoinTypes = async () => {
    try {
      const response = await axios.get('/api/data-merger/join-types');
      setJoinTypes(response.data.join_types);
    } catch (err) {
      console.error('Failed to load join types:', err);
    }
  };

  const getModeIcon = (value: string) => {
    switch (value) {
      case 'vertical':
        return <VerticalAlignCenterIcon />;
      case 'horizontal':
        return <CompareArrowsIcon />;
      case 'smart':
        return <AutoFixHighIcon />;
      default:
        return <MergeIcon />;
    }
  };

  const handleRemoveFile = (index: number) => {
    setFiles(files.filter((_, i) => i !== index));
  };

  const handlePreview = async () => {
    if (files.length < 2) {
      setPreviewError('Need at least 2 files to merge');
      return;
    }

    setPreviewLoading(true);
    setPreviewError('');

    try {
      const response = await axios.post('/api/data-merger/preview', {
        file_paths: files,
        mode,
        join_type: joinType,
        key_columns: keyColumns ? keyColumns.split(',').map(k => k.trim()) : null,
        add_source_column: addSourceColumn,
        n_rows: 10,
      });

      if (response.data.success) {
        setPreviewData(response.data.preview_data);
        setDetectedMode(response.data.detected_mode || mode);
        setEstimatedDimensions({
          rows: response.data.estimated_rows,
          columns: response.data.estimated_columns,
        });
      } else {
        setPreviewError(response.data.error || 'Preview failed');
      }
    } catch (err: any) {
      setPreviewError(err.response?.data?.detail || 'Failed to generate preview');
    } finally {
      setPreviewLoading(false);
    }
  };

  const handleMerge = async () => {
    if (files.length < 2) {
      setError('Need at least 2 files to merge');
      return;
    }

    setMerging(true);
    setError('');
    setProgress(0);

    try {
      // Start merge job
      const response = await axios.post('/api/data-merger/merge', {
        file_paths: files,
        output_filename: outputFilename,
        mode,
        join_type: joinType,
        key_columns: keyColumns ? keyColumns.split(',').map(k => k.trim()) : null,
        duplicate_handling: duplicateHandling,
        add_source_column: addSourceColumn,
        ignore_index: ignoreIndex,
        sort_keys: sortKeys,
      });

      const jobId = response.data.job_id;
      setJobId(jobId);

      // Poll for status
      pollJobStatus(jobId);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to start merge');
      setMerging(false);
    }
  };

  const pollJobStatus = async (jobId: string) => {
    const interval = setInterval(async () => {
      try {
        const response = await axios.get(`/api/data-merger/jobs/${jobId}`);
        const { status, progress: jobProgress, error: jobError } = response.data;

        setProgress(jobProgress);

        if (status === 'completed') {
          clearInterval(interval);
          setSuccess(true);
          setMerging(false);
          if (onSuccess) onSuccess();
        } else if (status === 'failed') {
          clearInterval(interval);
          setError(jobError || 'Merge failed');
          setMerging(false);
        }
      } catch (err) {
        clearInterval(interval);
        setError('Failed to check job status');
        setMerging(false);
      }
    }, 1000);
  };

  const handleDownload = async () => {
    if (!jobId) return;

    try {
      const response = await axios.get(`/api/data-merger/jobs/${jobId}/download`, {
        responseType: 'blob',
      });

      const url = window.URL.createObjectURL(new Blob([response.data]));
      const link = document.createElement('a');
      link.href = url;
      link.setAttribute('download', outputFilename);
      document.body.appendChild(link);
      link.click();
      link.remove();
    } catch (err) {
      setError('Failed to download merged file');
    }
  };

  const handleClose = () => {
    if (!merging) {
      onClose();
      // Reset state
      setPreviewData(null);
      setError('');
      setSuccess(false);
      setJobId('');
    }
  };

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="lg" fullWidth>
      <DialogTitle>
        <Box display="flex" alignItems="center" gap={1}>
          <MergeIcon />
          Data Merger
        </Box>
      </DialogTitle>

      <DialogContent>
        {/* File List */}
        <Box mb={3}>
          <Typography variant="subtitle2" gutterBottom>
            Files to Merge ({files.length})
          </Typography>
          <Paper variant="outlined" sx={{ maxHeight: 150, overflow: 'auto' }}>
            <List dense>
              {files.map((file, index) => (
                <ListItem key={index}>
                  <ListItemText
                    primary={file.split('/').pop()}
                    secondary={file}
                  />
                  <ListItemSecondaryAction>
                    <IconButton
                      edge="end"
                      onClick={() => handleRemoveFile(index)}
                      disabled={merging}
                    >
                      <DeleteIcon />
                    </IconButton>
                  </ListItemSecondaryAction>
                </ListItem>
              ))}
            </List>
          </Paper>
        </Box>

        {/* Merge Configuration */}
        <Box display="flex" flexDirection="column" gap={2}>
          {/* Merge Mode */}
          <FormControl fullWidth>
            <InputLabel>Merge Mode</InputLabel>
            <Select
              value={mode}
              onChange={(e) => setMode(e.target.value)}
              disabled={merging}
            >
              {mergeModes.map((m) => (
                <MenuItem key={m.value} value={m.value}>
                  <Box display="flex" alignItems="center" gap={1}>
                    {m.icon}
                    <Box>
                      <Typography variant="body2">{m.label}</Typography>
                      <Typography variant="caption" color="text.secondary">
                        {m.description}
                      </Typography>
                    </Box>
                  </Box>
                </MenuItem>
              ))}
            </Select>
          </FormControl>

          {/* Horizontal Merge Options */}
          {mode === 'horizontal' && (
            <>
              <FormControl fullWidth>
                <InputLabel>Join Type</InputLabel>
                <Select
                  value={joinType}
                  onChange={(e) => setJoinType(e.target.value)}
                  disabled={merging}
                >
                  {joinTypes.map((jt) => (
                    <MenuItem key={jt.value} value={jt.value}>
                      <Box>
                        <Typography variant="body2">{jt.label}</Typography>
                        <Typography variant="caption" color="text.secondary">
                          {jt.description}
                        </Typography>
                      </Box>
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>

              <TextField
                label="Key Columns (comma-separated)"
                value={keyColumns}
                onChange={(e) => setKeyColumns(e.target.value)}
                disabled={merging}
                placeholder="e.g., id, sample_id"
                helperText="Columns used to join the files"
              />
            </>
          )}

          {/* Additional Options */}
          <Box display="flex" flexDirection="column">
            <FormControlLabel
              control={
                <Checkbox
                  checked={addSourceColumn}
                  onChange={(e) => setAddSourceColumn(e.target.checked)}
                  disabled={merging}
                />
              }
              label="Add source file column"
            />
            <FormControlLabel
              control={
                <Checkbox
                  checked={ignoreIndex}
                  onChange={(e) => setIgnoreIndex(e.target.checked)}
                  disabled={merging}
                />
              }
              label="Reset row index"
            />
            {mode === 'horizontal' && (
              <FormControlLabel
                control={
                  <Checkbox
                    checked={sortKeys}
                    onChange={(e) => setSortKeys(e.target.checked)}
                    disabled={merging}
                  />
                }
                label="Sort by key columns"
              />
            )}
          </Box>

          {/* Output Filename */}
          <TextField
            label="Output Filename"
            value={outputFilename}
            onChange={(e) => setOutputFilename(e.target.value)}
            disabled={merging}
          />
        </Box>

        {/* Preview Button */}
        <Box mt={2} mb={2}>
          <Button
            variant="outlined"
            startIcon={<PreviewIcon />}
            onClick={handlePreview}
            disabled={previewLoading || merging || files.length < 2}
            fullWidth
          >
            {previewLoading ? 'Loading Preview...' : 'Preview Merge'}
          </Button>
        </Box>

        {/* Preview Data */}
        {previewData && (
          <Box mt={2}>
            <Box display="flex" justifyContent="space-between" alignItems="center" mb={1}>
              <Typography variant="subtitle2">
                Preview (First 10 rows)
              </Typography>
              {detectedMode && (
                <Chip
                  label={`Mode: ${detectedMode}`}
                  size="small"
                  color="primary"
                />
              )}
            </Box>
            <Typography variant="caption" color="text.secondary" gutterBottom>
              Estimated: {estimatedDimensions.rows} rows × {estimatedDimensions.columns} columns
            </Typography>
            <Paper variant="outlined" sx={{ maxHeight: 300, overflow: 'auto', mt: 1 }}>
              <Table size="small">
                <TableHead>
                  <TableRow>
                    {previewData.columns.map((col, idx) => (
                      <TableCell key={idx}>{col}</TableCell>
                    ))}
                  </TableRow>
                </TableHead>
                <TableBody>
                  {previewData.rows.map((row, idx) => (
                    <TableRow key={idx}>
                      {row.map((cell, cellIdx) => (
                        <TableCell key={cellIdx}>
                          {cell === null || cell === undefined ? 'NaN' : String(cell)}
                        </TableCell>
                      ))}
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </Paper>
          </Box>
        )}

        {/* Preview Error */}
        {previewError && (
          <Alert severity="error" sx={{ mt: 2 }}>
            {previewError}
          </Alert>
        )}

        {/* Progress */}
        {merging && (
          <Box mt={2}>
            <Typography variant="body2" gutterBottom>
              Merging files... {progress}%
            </Typography>
            <LinearProgress variant="determinate" value={progress} />
          </Box>
        )}

        {/* Error */}
        {error && (
          <Alert severity="error" sx={{ mt: 2 }}>
            {error}
          </Alert>
        )}

        {/* Success */}
        {success && (
          <Alert severity="success" sx={{ mt: 2 }}>
            Files merged successfully! You can download the result.
          </Alert>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={handleClose} disabled={merging}>
          Cancel
        </Button>
        {success && (
          <Button
            variant="contained"
            color="primary"
            onClick={handleDownload}
          >
            Download Result
          </Button>
        )}
        {!success && (
          <Button
            variant="contained"
            color="primary"
            onClick={handleMerge}
            disabled={merging || files.length < 2}
            startIcon={<MergeIcon />}
          >
            {merging ? 'Merging...' : 'Merge Files'}
          </Button>
        )}
      </DialogActions>
    </Dialog>
  );
};

export default DataMergerModal;
