/**
 * Create Dataset from Run Dialog
 * 
 * Dialog for creating a dataset from a completed pipeline run's outputs.
 * Automatically adds all output files and records lineage.
 */

import React, { useState, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  TextField,
  Box,
  Typography,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Chip,
  Alert,
  CircularProgress,
  Autocomplete,
  List,
  ListItem,
  ListItemText,
  Divider,
} from '@mui/material';
import {
  SaveOutlined,
  AutoAwesomeOutlined,
  FolderOutlined,
} from '@mui/icons-material';
import axios from 'axios';

interface Run {
  id: number;
  name: string;
  status: string;
  pipeline_type: string;
  output_files: number[];
}

interface CreateDatasetFromRunDialogProps {
  open: boolean;
  onClose: () => void;
  run: Run | null;
  onSuccess?: (datasetId: number) => void;
}

const DATA_TYPES = [
  'genomics',
  'transcriptomics',
  'proteomics',
  'metabolomics',
  'single-cell',
  'general',
];

const COMMON_TAGS = [
  'pipeline-output',
  'processed',
  'analysis-ready',
  'qc-passed',
  'normalized',
  'filtered',
];

const CreateDatasetFromRunDialog: React.FC<CreateDatasetFromRunDialogProps> = ({
  open,
  onClose,
  run,
  onSuccess,
}) => {
  const [datasetName, setDatasetName] = useState('');
  const [description, setDescription] = useState('');
  const [dataType, setDataType] = useState('general');
  const [selectedTags, setSelectedTags] = useState<string[]>(['pipeline-output']);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [autoMode, setAutoMode] = useState(false);

  // Initialize form with run data
  useEffect(() => {
    if (run && open) {
      setDatasetName(`${run.name}_output`);
      setDescription(`Output dataset from pipeline run: ${run.name}`);
      
      // Infer data type from pipeline type
      const pipelineType = run.pipeline_type?.toLowerCase() || '';
      if (pipelineType.includes('rna') || pipelineType.includes('transcriptomics')) {
        setDataType('transcriptomics');
      } else if (pipelineType.includes('dna') || pipelineType.includes('genomics')) {
        setDataType('genomics');
      } else if (pipelineType.includes('protein') || pipelineType.includes('proteomics')) {
        setDataType('proteomics');
      } else if (pipelineType.includes('metabol')) {
        setDataType('metabolomics');
      } else if (pipelineType.includes('single') && pipelineType.includes('cell')) {
        setDataType('single-cell');
      } else {
        setDataType('general');
      }

      // Add pipeline type tag
      if (run.pipeline_type) {
        setSelectedTags(['pipeline-output', `pipeline-${run.pipeline_type}`]);
      }
    }
  }, [run, open]);

  const handleCreateDataset = async () => {
    if (!run) return;

    setLoading(true);
    setError(null);

    try {
      const response = await axios.post('/api/pipeline-datasets/create-from-run', {
        run_id: run.id,
        dataset_name: datasetName,
        dataset_description: description,
        data_type: dataType,
        tags: selectedTags,
      });

      if (onSuccess) {
        onSuccess(response.data.dataset.id);
      }

      handleClose();
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to create dataset');
    } finally {
      setLoading(false);
    }
  };

  const handleAutoCreate = async () => {
    if (!run) return;

    setLoading(true);
    setError(null);
    setAutoMode(true);

    try {
      const response = await axios.post('/api/pipeline-datasets/auto-create', {
        run_id: run.id,
        auto_tags: true,
      });

      if (onSuccess) {
        onSuccess(response.data.dataset.id);
      }

      handleClose();
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to auto-create dataset');
    } finally {
      setLoading(false);
      setAutoMode(false);
    }
  };

  const handleClose = () => {
    setDatasetName('');
    setDescription('');
    setDataType('general');
    setSelectedTags(['pipeline-output']);
    setError(null);
    onClose();
  };

  if (!run) return null;

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="md" fullWidth>
      <DialogTitle>
        <Box display="flex" alignItems="center" gap={1}>
          <SaveOutlined />
          <Typography variant="h6">Create Dataset from Run</Typography>
        </Box>
      </DialogTitle>

      <DialogContent dividers>
        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        {/* Run Information */}
        <Box mb={3} p={2} bgcolor="grey.50" borderRadius={1}>
          <Typography variant="subtitle2" gutterBottom>
            Pipeline Run Information
          </Typography>
          <List dense>
            <ListItem>
              <ListItemText primary="Run Name" secondary={run.name} />
            </ListItem>
            <ListItem>
              <ListItemText primary="Pipeline Type" secondary={run.pipeline_type} />
            </ListItem>
            <ListItem>
              <ListItemText
                primary="Output Files"
                secondary={`${run.output_files?.length || 0} files`}
              />
            </ListItem>
          </List>
        </Box>

        <Divider sx={{ my: 2 }} />

        {/* Dataset Configuration */}
        <Box display="flex" flexDirection="column" gap={2}>
          <TextField
            label="Dataset Name"
            value={datasetName}
            onChange={(e) => setDatasetName(e.target.value)}
            fullWidth
            required
            helperText="Unique name for the dataset"
          />

          <TextField
            label="Description"
            value={description}
            onChange={(e) => setDescription(e.target.value)}
            fullWidth
            multiline
            rows={3}
            helperText="Describe the contents and purpose of this dataset"
          />

          <FormControl fullWidth>
            <InputLabel>Data Type</InputLabel>
            <Select
              value={dataType}
              onChange={(e) => setDataType(e.target.value)}
              label="Data Type"
            >
              {DATA_TYPES.map((type) => (
                <MenuItem key={type} value={type}>
                  {type.charAt(0).toUpperCase() + type.slice(1)}
                </MenuItem>
              ))}
            </Select>
          </FormControl>

          <Autocomplete
            multiple
            options={COMMON_TAGS}
            value={selectedTags}
            onChange={(_, newValue) => setSelectedTags(newValue)}
            freeSolo
            renderTags={(value, getTagProps) =>
              value.map((option, index) => (
                <Chip label={option} {...getTagProps({ index })} />
              ))
            }
            renderInput={(params) => (
              <TextField
                {...params}
                label="Tags"
                helperText="Add tags to organize and categorize this dataset"
              />
            )}
          />
        </Box>

        <Alert severity="info" sx={{ mt: 3 }} icon={<FolderOutlined />}>
          All output files from this run will be automatically added to the dataset with
          hash computation for integrity verification.
        </Alert>
      </DialogContent>

      <DialogActions>
        <Button onClick={handleClose} disabled={loading}>
          Cancel
        </Button>
        <Button
          onClick={handleAutoCreate}
          startIcon={<AutoAwesomeOutlined />}
          disabled={loading}
          color="secondary"
        >
          {autoMode && loading ? (
            <CircularProgress size={20} />
          ) : (
            'Auto Create'
          )}
        </Button>
        <Button
          onClick={handleCreateDataset}
          variant="contained"
          disabled={!datasetName || loading}
          startIcon={<SaveOutlined />}
        >
          {!autoMode && loading ? (
            <CircularProgress size={20} color="inherit" />
          ) : (
            'Create Dataset'
          )}
        </Button>
      </DialogActions>
    </Dialog>
  );
};

export default CreateDatasetFromRunDialog;
