import { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Loader2, PlayCircle, CheckCircle, AlertCircle } from 'lucide-react';

interface WorkflowResult {
  workflow_id: number;
  status: string;
  message: string;
}

export default function MetabolomicsAnalysis() {
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<WorkflowResult | null>(null);
  const [error, setError] = useState<string | null>(null);

  // Feature Detection state
  const [detectionForm, setDetectionForm] = useState({
    mzml_files: '',
    output_dir: '',
    tool: 'xcms',
    ppm: 25,
    peakwidth_min: 10,
    peakwidth_max: 60,
  });

  // Annotation state
  const [annotationForm, setAnnotationForm] = useState({
    feature_table: '',
    output_dir: '',
    tool: 'gnps',
    library_file: '',
    mz_tolerance: 0.02,
    rt_tolerance: 10,
  });

  // Quantification state
  const [quantForm, setQuantForm] = useState({
    feature_table: '',
    sample_metadata: '',
    output_dir: '',
    normalization: 'median',
    imputation: 'knn',
    scaling: 'pareto',
  });

  const runFeatureDetection = async () => {
    setLoading(true);
    setError(null);
    try {
      const response = await fetch('/api/metabolomics/feature-detection', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify({
          sample_id: 1,
          ...detectionForm,
          mzml_files: detectionForm.mzml_files.split(',').map(f => f.trim()),
        }),
      });

      if (!response.ok) {
        throw new Error('Feature detection request failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error');
    } finally {
      setLoading(false);
    }
  };

  const runAnnotation = async () => {
    setLoading(true);
    setError(null);
    try {
      const response = await fetch('/api/metabolomics/annotation', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify({
          sample_id: 1,
          ...annotationForm,
        }),
      });

      if (!response.ok) {
        throw new Error('Annotation request failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error');
    } finally {
      setLoading(false);
    }
  };

  const runQuantification = async () => {
    setLoading(true);
    setError(null);
    try {
      const response = await fetch('/api/metabolomics/quantification', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify({
          sample_id: 1,
          ...quantForm,
        }),
      });

      if (!response.ok) {
        throw new Error('Quantification request failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container mx-auto py-6 space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-3xl font-bold">Metabolomics Analysis</h1>
          <p className="text-muted-foreground">
            Feature detection, annotation, and quantification for metabolomics data
          </p>
        </div>
      </div>

      {result && (
        <Alert>
          <CheckCircle className="h-4 w-4" />
          <AlertDescription>
            Workflow #{result.workflow_id} {result.status}: {result.message}
          </AlertDescription>
        </Alert>
      )}

      {error && (
        <Alert variant="destructive">
          <AlertCircle className="h-4 w-4" />
          <AlertDescription>{error}</AlertDescription>
        </Alert>
      )}

      <Tabs defaultValue="detection" className="w-full">
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="detection">Feature Detection</TabsTrigger>
          <TabsTrigger value="annotation">Annotation</TabsTrigger>
          <TabsTrigger value="quantification">Quantification</TabsTrigger>
        </TabsList>

        <TabsContent value="detection">
          <Card>
            <CardHeader>
              <CardTitle>Feature Detection</CardTitle>
              <CardDescription>
                Detect metabolite features from MS data using XCMS or MZmine
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-2">
                <Label htmlFor="mzml_files">mzML Files (comma-separated)</Label>
                <Input
                  id="mzml_files"
                  value={detectionForm.mzml_files}
                  onChange={(e) => setDetectionForm({ ...detectionForm, mzml_files: e.target.value })}
                  placeholder="/data/sample1.mzML, /data/sample2.mzML"
                />
              </div>

              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="detection_output">Output Directory</Label>
                  <Input
                    id="detection_output"
                    value={detectionForm.output_dir}
                    onChange={(e) => setDetectionForm({ ...detectionForm, output_dir: e.target.value })}
                    placeholder="/output/features"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="tool">Detection Tool</Label>
                  <Select
                    value={detectionForm.tool}
                    onValueChange={(value) => setDetectionForm({ ...detectionForm, tool: value })}
                  >
                    <SelectTrigger id="tool">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="xcms">XCMS</SelectItem>
                      <SelectItem value="mzmine">MZmine</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
              </div>

              <div className="grid grid-cols-3 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="ppm">PPM Tolerance</Label>
                  <Input
                    id="ppm"
                    type="number"
                    value={detectionForm.ppm}
                    onChange={(e) => setDetectionForm({ ...detectionForm, ppm: parseFloat(e.target.value) })}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="peakwidth_min">Peak Width Min (s)</Label>
                  <Input
                    id="peakwidth_min"
                    type="number"
                    value={detectionForm.peakwidth_min}
                    onChange={(e) => setDetectionForm({ ...detectionForm, peakwidth_min: parseFloat(e.target.value) })}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="peakwidth_max">Peak Width Max (s)</Label>
                  <Input
                    id="peakwidth_max"
                    type="number"
                    value={detectionForm.peakwidth_max}
                    onChange={(e) => setDetectionForm({ ...detectionForm, peakwidth_max: parseFloat(e.target.value) })}
                  />
                </div>
              </div>

              <Button onClick={runFeatureDetection} disabled={loading} className="w-full">
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <PlayCircle className="mr-2 h-4 w-4" />
                    Run Feature Detection
                  </>
                )}
              </Button>
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="annotation">
          <Card>
            <CardHeader>
              <CardTitle>Metabolite Annotation</CardTitle>
              <CardDescription>
                Annotate detected features using spectral libraries (GNPS or MS-DIAL)
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-2">
                <Label htmlFor="feature_table">Feature Table</Label>
                <Input
                  id="feature_table"
                  value={annotationForm.feature_table}
                  onChange={(e) => setAnnotationForm({ ...annotationForm, feature_table: e.target.value })}
                  placeholder="/output/features/feature_table.csv"
                />
              </div>

              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="annotation_output">Output Directory</Label>
                  <Input
                    id="annotation_output"
                    value={annotationForm.output_dir}
                    onChange={(e) => setAnnotationForm({ ...annotationForm, output_dir: e.target.value })}
                    placeholder="/output/annotation"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="annotation_tool">Annotation Tool</Label>
                  <Select
                    value={annotationForm.tool}
                    onValueChange={(value) => setAnnotationForm({ ...annotationForm, tool: value })}
                  >
                    <SelectTrigger id="annotation_tool">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="gnps">GNPS</SelectItem>
                      <SelectItem value="msdial">MS-DIAL</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
              </div>

              <div className="space-y-2">
                <Label htmlFor="library_file">Spectral Library (optional for GNPS)</Label>
                <Input
                  id="library_file"
                  value={annotationForm.library_file}
                  onChange={(e) => setAnnotationForm({ ...annotationForm, library_file: e.target.value })}
                  placeholder="/data/library.mgf"
                />
              </div>

              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="mz_tolerance">m/z Tolerance</Label>
                  <Input
                    id="mz_tolerance"
                    type="number"
                    step="0.001"
                    value={annotationForm.mz_tolerance}
                    onChange={(e) => setAnnotationForm({ ...annotationForm, mz_tolerance: parseFloat(e.target.value) })}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="rt_tolerance">RT Tolerance (s)</Label>
                  <Input
                    id="rt_tolerance"
                    type="number"
                    value={annotationForm.rt_tolerance}
                    onChange={(e) => setAnnotationForm({ ...annotationForm, rt_tolerance: parseFloat(e.target.value) })}
                  />
                </div>
              </div>

              <Button onClick={runAnnotation} disabled={loading} className="w-full">
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <PlayCircle className="mr-2 h-4 w-4" />
                    Run Annotation
                  </>
                )}
              </Button>
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="quantification">
          <Card>
            <CardHeader>
              <CardTitle>Feature Quantification</CardTitle>
              <CardDescription>
                Normalize, impute, and scale metabolite quantification data
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-2">
                <Label htmlFor="quant_feature_table">Feature Table</Label>
                <Input
                  id="quant_feature_table"
                  value={quantForm.feature_table}
                  onChange={(e) => setQuantForm({ ...quantForm, feature_table: e.target.value })}
                  placeholder="/output/features/feature_table.csv"
                />
              </div>

              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="sample_metadata">Sample Metadata</Label>
                  <Input
                    id="sample_metadata"
                    value={quantForm.sample_metadata}
                    onChange={(e) => setQuantForm({ ...quantForm, sample_metadata: e.target.value })}
                    placeholder="/data/sample_metadata.csv"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="quant_output">Output Directory</Label>
                  <Input
                    id="quant_output"
                    value={quantForm.output_dir}
                    onChange={(e) => setQuantForm({ ...quantForm, output_dir: e.target.value })}
                    placeholder="/output/quantification"
                  />
                </div>
              </div>

              <div className="grid grid-cols-3 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="normalization">Normalization</Label>
                  <Select
                    value={quantForm.normalization}
                    onValueChange={(value) => setQuantForm({ ...quantForm, normalization: value })}
                  >
                    <SelectTrigger id="normalization">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="median">Median</SelectItem>
                      <SelectItem value="quantile">Quantile</SelectItem>
                      <SelectItem value="total">Total Sum</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
                <div className="space-y-2">
                  <Label htmlFor="imputation">Imputation</Label>
                  <Select
                    value={quantForm.imputation}
                    onValueChange={(value) => setQuantForm({ ...quantForm, imputation: value })}
                  >
                    <SelectTrigger id="imputation">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="knn">KNN</SelectItem>
                      <SelectItem value="zero">Zero</SelectItem>
                      <SelectItem value="min">Minimum</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
                <div className="space-y-2">
                  <Label htmlFor="scaling">Scaling</Label>
                  <Select
                    value={quantForm.scaling}
                    onValueChange={(value) => setQuantForm({ ...quantForm, scaling: value })}
                  >
                    <SelectTrigger id="scaling">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="pareto">Pareto</SelectItem>
                      <SelectItem value="auto">Auto (UV)</SelectItem>
                      <SelectItem value="range">Range</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
              </div>

              <Button onClick={runQuantification} disabled={loading} className="w-full">
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <PlayCircle className="mr-2 h-4 w-4" />
                    Run Quantification
                  </>
                )}
              </Button>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </div>
  );
}
