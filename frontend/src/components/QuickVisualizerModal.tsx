import React, { useState, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Box,
  Typography,
  Alert,
  CircularProgress,
  Card,
  CardContent,
  Grid,
  Chip,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  TextField,
  Divider,
  Paper,
} from '@mui/material';
import {
  BarChart as BarChartIcon,
  ShowChart as LineChartIcon,
  ScatterPlot as ScatterIcon,
  Dashboard as HeatmapIcon,
  Timeline as HistogramIcon,
  ViewModule as BoxPlotIcon,
} from '@mui/icons-material';
import axios from 'axios';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:8000';

interface QuickVisualizerModalProps {
  open: boolean;
  onClose: () => void;
  filePath: string;
  fileName: string;
}

interface ColumnInfo {
  name: string;
  dtype: string;
  n_unique: number;
  n_missing: number;
  is_numeric: boolean;
  is_categorical: boolean;
  min?: number;
  max?: number;
  mean?: number;
  median?: number;
}

interface Suggestion {
  chart_type: string;
  title: string;
  description: string;
  options: Record<string, any>;
  priority: number;
}

interface AnalyzeResult {
  success: boolean;
  error?: string;
  data_type?: string;
  n_rows?: number;
  n_columns?: number;
  limited?: boolean;
  column_info?: ColumnInfo[];
  suggestions?: Suggestion[];
}

const CHART_ICONS: Record<string, React.ReactElement> = {
  histogram: <HistogramIcon />,
  scatter: <ScatterIcon />,
  line: <LineChartIcon />,
  bar: <BarChartIcon />,
  boxplot: <BoxPlotIcon />,
  heatmap: <HeatmapIcon />,
  violin: <BoxPlotIcon />,
  correlation: <HeatmapIcon />,
};

export const QuickVisualizerModal: React.FC<QuickVisualizerModalProps> = ({
  open,
  onClose,
  filePath,
  fileName,
}) => {
  const [loading, setLoading] = useState(false);
  const [analyzing, setAnalyzing] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  const [analyzeResult, setAnalyzeResult] = useState<AnalyzeResult | null>(null);
  const [selectedChartType, setSelectedChartType] = useState<string>('');
  const [chartOptions, setChartOptions] = useState<Record<string, any>>({});
  const [generatedImage, setGeneratedImage] = useState<string | null>(null);

  // Auto-analyze when modal opens
  useEffect(() => {
    if (open && filePath) {
      handleAnalyze();
    } else {
      // Reset state when modal closes
      setAnalyzeResult(null);
      setSelectedChartType('');
      setChartOptions({});
      setGeneratedImage(null);
      setError(null);
    }
  }, [open, filePath]);

  const handleAnalyze = async () => {
    setAnalyzing(true);
    setError(null);
    
    try {
      const response = await axios.post<AnalyzeResult>(
        `${API_BASE_URL}/api/quick-visualizer/analyze`,
        { file_path: filePath }
      );
      
      if (response.data.success) {
        setAnalyzeResult(response.data);
        
        // Auto-select first suggestion
        if (response.data.suggestions && response.data.suggestions.length > 0) {
          const firstSuggestion = response.data.suggestions[0];
          setSelectedChartType(firstSuggestion.chart_type);
          setChartOptions(firstSuggestion.options || {});
        }
      } else {
        setError(response.data.error || 'Analysis failed');
      }
    } catch (err: any) {
      console.error('Analysis error:', err);
      setError(err.response?.data?.detail || err.message || 'Failed to analyze file');
    } finally {
      setAnalyzing(false);
    }
  };

  const handleGenerateVisualization = async () => {
    if (!selectedChartType) {
      setError('Please select a chart type');
      return;
    }

    setLoading(true);
    setError(null);
    setGeneratedImage(null);
    
    try {
      const response = await axios.post(
        `${API_BASE_URL}/api/quick-visualizer/visualize`,
        {
          file_path: filePath,
          chart_type: selectedChartType,
          options: chartOptions,
        }
      );
      
      if (response.data.success) {
        setGeneratedImage(response.data.image);
      } else {
        setError(response.data.error || 'Visualization failed');
      }
    } catch (err: any) {
      console.error('Visualization error:', err);
      setError(err.response?.data?.detail || err.message || 'Failed to generate visualization');
    } finally {
      setLoading(false);
    }
  };

  const handleOptionChange = (key: string, value: any) => {
    setChartOptions(prev => ({
      ...prev,
      [key]: value,
    }));
  };

  const getNumericColumns = (): string[] => {
    if (!analyzeResult?.column_info) return [];
    return analyzeResult.column_info
      .filter(col => col.is_numeric)
      .map(col => col.name);
  };

  const getCategoricalColumns = (): string[] => {
    if (!analyzeResult?.column_info) return [];
    return analyzeResult.column_info
      .filter(col => col.is_categorical)
      .map(col => col.name);
  };

  const renderChartOptions = () => {
    if (!selectedChartType) return null;

    const numericCols = getNumericColumns();
    const categoricalCols = getCategoricalColumns();

    switch (selectedChartType) {
      case 'histogram':
        return (
          <FormControl fullWidth margin="normal">
            <InputLabel>Column</InputLabel>
            <Select
              value={chartOptions.column || ''}
              onChange={(e) => handleOptionChange('column', e.target.value)}
            >
              {numericCols.map(col => (
                <MenuItem key={col} value={col}>{col}</MenuItem>
              ))}
            </Select>
          </FormControl>
        );

      case 'scatter':
        return (
          <>
            <FormControl fullWidth margin="normal">
              <InputLabel>X Column</InputLabel>
              <Select
                value={chartOptions.x_column || ''}
                onChange={(e) => handleOptionChange('x_column', e.target.value)}
              >
                {numericCols.map(col => (
                  <MenuItem key={col} value={col}>{col}</MenuItem>
                ))}
              </Select>
            </FormControl>
            <FormControl fullWidth margin="normal">
              <InputLabel>Y Column</InputLabel>
              <Select
                value={chartOptions.y_column || ''}
                onChange={(e) => handleOptionChange('y_column', e.target.value)}
              >
                {numericCols.map(col => (
                  <MenuItem key={col} value={col}>{col}</MenuItem>
                ))}
              </Select>
            </FormControl>
          </>
        );

      case 'bar':
        return (
          <FormControl fullWidth margin="normal">
            <InputLabel>Column</InputLabel>
            <Select
              value={chartOptions.column || ''}
              onChange={(e) => handleOptionChange('column', e.target.value)}
            >
              {categoricalCols.map(col => (
                <MenuItem key={col} value={col}>{col}</MenuItem>
              ))}
            </Select>
          </FormControl>
        );

      case 'boxplot':
      case 'violin':
        return (
          <TextField
            fullWidth
            margin="normal"
            label="Columns (comma-separated)"
            placeholder="col1,col2,col3"
            value={Array.isArray(chartOptions.columns) ? chartOptions.columns.join(',') : chartOptions.columns || ''}
            onChange={(e) => {
              const cols = e.target.value.split(',').map(c => c.trim()).filter(c => c);
              handleOptionChange('columns', cols);
            }}
            helperText={`Available: ${numericCols.slice(0, 5).join(', ')}...`}
          />
        );

      default:
        return null;
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="lg" fullWidth>
      <DialogTitle>
        <Box display="flex" alignItems="center" gap={1}>
          <ScatterIcon />
          <Typography variant="h6">Quick Visualization</Typography>
        </Box>
        <Typography variant="body2" color="text.secondary">
          {fileName}
        </Typography>
      </DialogTitle>

      <DialogContent>
        {error && (
          <Alert severity="error" sx={{ mb: 2 }} onClose={() => setError(null)}>
            {error}
          </Alert>
        )}

        {analyzing && (
          <Box display="flex" justifyContent="center" alignItems="center" py={4}>
            <CircularProgress />
            <Typography variant="body2" sx={{ ml: 2 }}>
              Analyzing file...
            </Typography>
          </Box>
        )}

        {analyzeResult && !analyzing && (
          <Grid container spacing={2}>
            {/* File Info */}
            <Grid item xs={12}>
              <Paper variant="outlined" sx={{ p: 2 }}>
                <Typography variant="subtitle2" gutterBottom>
                  File Information
                </Typography>
                <Box display="flex" gap={2} flexWrap="wrap">
                  <Chip
                    label={`Type: ${analyzeResult.data_type}`}
                    color="primary"
                    size="small"
                  />
                  <Chip
                    label={`${analyzeResult.n_rows} rows × ${analyzeResult.n_columns} columns`}
                    size="small"
                  />
                  {analyzeResult.limited && (
                    <Chip
                      label="Preview only (first 10,000 rows)"
                      color="warning"
                      size="small"
                    />
                  )}
                </Box>
              </Paper>
            </Grid>

            {/* Suggestions */}
            {analyzeResult.suggestions && analyzeResult.suggestions.length > 0 && (
              <Grid item xs={12}>
                <Typography variant="subtitle2" gutterBottom>
                  Suggested Visualizations
                </Typography>
                <Box display="flex" gap={1} flexWrap="wrap">
                  {analyzeResult.suggestions.map((suggestion) => (
                    <Card
                      key={suggestion.chart_type}
                      variant="outlined"
                      sx={{
                        cursor: 'pointer',
                        border: selectedChartType === suggestion.chart_type ? 2 : 1,
                        borderColor: selectedChartType === suggestion.chart_type 
                          ? 'primary.main' 
                          : 'divider',
                        '&:hover': {
                          borderColor: 'primary.main',
                          bgcolor: 'action.hover',
                        },
                      }}
                      onClick={() => {
                        setSelectedChartType(suggestion.chart_type);
                        setChartOptions(suggestion.options || {});
                        setGeneratedImage(null);
                      }}
                    >
                      <CardContent sx={{ p: 1.5, '&:last-child': { pb: 1.5 } }}>
                        <Box display="flex" alignItems="center" gap={1}>
                          {CHART_ICONS[suggestion.chart_type]}
                          <Box>
                            <Typography variant="body2" fontWeight="medium">
                              {suggestion.title}
                            </Typography>
                            <Typography variant="caption" color="text.secondary">
                              {suggestion.description}
                            </Typography>
                          </Box>
                        </Box>
                      </CardContent>
                    </Card>
                  ))}
                </Box>
              </Grid>
            )}

            {/* Chart Options */}
            {selectedChartType && (
              <Grid item xs={12}>
                <Divider sx={{ my: 1 }} />
                <Typography variant="subtitle2" gutterBottom>
                  Chart Options
                </Typography>
                {renderChartOptions()}
              </Grid>
            )}

            {/* Generated Visualization */}
            {generatedImage && (
              <Grid item xs={12}>
                <Divider sx={{ my: 1 }} />
                <Typography variant="subtitle2" gutterBottom>
                  Visualization
                </Typography>
                <Box
                  component="img"
                  src={generatedImage}
                  alt="Generated visualization"
                  sx={{
                    width: '100%',
                    maxHeight: '600px',
                    objectFit: 'contain',
                    border: '1px solid',
                    borderColor: 'divider',
                    borderRadius: 1,
                  }}
                />
              </Grid>
            )}
          </Grid>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose}>Close</Button>
        {analyzeResult && (
          <Button
            variant="contained"
            onClick={handleGenerateVisualization}
            disabled={loading || !selectedChartType}
            startIcon={loading ? <CircularProgress size={20} /> : null}
          >
            {loading ? 'Generating...' : 'Generate Visualization'}
          </Button>
        )}
      </DialogActions>
    </Dialog>
  );
};

export default QuickVisualizerModal;
