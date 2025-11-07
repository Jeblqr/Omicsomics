import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Checkbox } from '@/components/ui/checkbox';
import { Loader2, Grid3x3, Users, TrendingUp } from 'lucide-react';

interface QCParams {
  sampleId: string;
  minGenes: number;
  maxGenes: number;
  maxMtPercent: number;
  minCells: number;
}

interface ClusterParams {
  sampleId: string;
  nPcs: number;
  resolution: number;
  algorithm: 'louvain' | 'leiden';
  nNeighbors: number;
}

interface DEParams {
  sampleId: string;
  cluster1: string;
  cluster2: string;
  method: 'wilcoxon' | 'ttest' | 'deseq2';
  minPct: number;
  logfcThreshold: number;
}

interface TrajectoryParams {
  sampleId: string;
  rootCell: string;
  method: 'monocle3' | 'slingshot' | 'paga';
  nComponents: number;
}

const SingleCellAnalysis: React.FC = () => {
  const [activeTab, setActiveTab] = useState<string>('qc');
  const [loading, setLoading] = useState<boolean>(false);
  const [result, setResult] = useState<any>(null);
  const [error, setError] = useState<string>('');

  // Quality Control State
  const [qcParams, setQcParams] = useState<QCParams>({
    sampleId: '',
    minGenes: 200,
    maxGenes: 5000,
    maxMtPercent: 10,
    minCells: 3
  });

  // Clustering State
  const [clusterParams, setClusterParams] = useState<ClusterParams>({
    sampleId: '',
    nPcs: 20,
    resolution: 0.5,
    algorithm: 'louvain',
    nNeighbors: 15
  });

  // Differential Expression State
  const [deParams, setDeParams] = useState<DEParams>({
    sampleId: '',
    cluster1: '',
    cluster2: '',
    method: 'wilcoxon',
    minPct: 0.25,
    logfcThreshold: 0.5
  });

  // Trajectory Analysis State
  const [trajectoryParams, setTrajectoryParams] = useState<TrajectoryParams>({
    sampleId: '',
    rootCell: '',
    method: 'monocle3',
    nComponents: 2
  });

  const handleQualityControl = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/single-cell/qc', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(qcParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Quality control failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during quality control');
    } finally {
      setLoading(false);
    }
  };

  const handleClustering = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/single-cell/clustering', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(clusterParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Clustering failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during clustering');
    } finally {
      setLoading(false);
    }
  };

  const handleDifferentialExpression = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/single-cell/differential', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(deParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Differential expression analysis failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during differential expression analysis');
    } finally {
      setLoading(false);
    }
  };

  const handleTrajectoryAnalysis = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/single-cell/trajectory', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(trajectoryParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Trajectory analysis failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during trajectory analysis');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container mx-auto py-8 px-4">
      <div className="mb-8">
        <h1 className="text-4xl font-bold mb-2">Single-Cell RNA-seq Analysis</h1>
        <p className="text-muted-foreground">
          Perform quality control, clustering, differential expression, and trajectory analysis on scRNA-seq data
        </p>
      </div>

      <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
        <TabsList className="grid w-full grid-cols-4">
          <TabsTrigger value="qc" className="flex items-center gap-2">
            <Grid3x3 className="h-4 w-4" />
            Quality Control
          </TabsTrigger>
          <TabsTrigger value="clustering" className="flex items-center gap-2">
            <Users className="h-4 w-4" />
            Clustering
          </TabsTrigger>
          <TabsTrigger value="de" className="flex items-center gap-2">
            <TrendingUp className="h-4 w-4" />
            Differential Expression
          </TabsTrigger>
          <TabsTrigger value="trajectory" className="flex items-center gap-2">
            <TrendingUp className="h-4 w-4" />
            Trajectory
          </TabsTrigger>
        </TabsList>

        {/* Quality Control Tab */}
        <TabsContent value="qc" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Quality Control & Filtering</CardTitle>
              <CardDescription>
                Filter cells and genes based on quality metrics (gene count, UMI count, mitochondrial percentage)
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="qc-sample">Sample ID</Label>
                  <Input
                    id="qc-sample"
                    placeholder="Enter sample ID"
                    value={qcParams.sampleId}
                    onChange={(e) => setQcParams({ ...qcParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="min-genes">Min Genes per Cell</Label>
                    <Input
                      id="min-genes"
                      type="number"
                      value={qcParams.minGenes}
                      onChange={(e) => setQcParams({ ...qcParams, minGenes: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Remove low-quality cells</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="max-genes">Max Genes per Cell</Label>
                    <Input
                      id="max-genes"
                      type="number"
                      value={qcParams.maxGenes}
                      onChange={(e) => setQcParams({ ...qcParams, maxGenes: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Remove potential doublets</p>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="max-mt">Max Mitochondrial % </Label>
                    <Input
                      id="max-mt"
                      type="number"
                      step="0.5"
                      value={qcParams.maxMtPercent}
                      onChange={(e) => setQcParams({ ...qcParams, maxMtPercent: parseFloat(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Remove dying/stressed cells</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="min-cells">Min Cells per Gene</Label>
                    <Input
                      id="min-cells"
                      type="number"
                      value={qcParams.minCells}
                      onChange={(e) => setQcParams({ ...qcParams, minCells: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Remove rare genes</p>
                  </div>
                </div>

                <Button
                  onClick={handleQualityControl}
                  disabled={loading || !qcParams.sampleId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Running QC...
                    </>
                  ) : (
                    'Run Quality Control'
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Clustering Tab */}
        <TabsContent value="clustering" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Cell Clustering</CardTitle>
              <CardDescription>
                Identify cell populations using graph-based clustering (Louvain or Leiden algorithm)
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="cluster-sample">Sample ID</Label>
                  <Input
                    id="cluster-sample"
                    placeholder="Enter sample ID"
                    value={clusterParams.sampleId}
                    onChange={(e) => setClusterParams({ ...clusterParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="n-pcs">Number of PCs</Label>
                    <Input
                      id="n-pcs"
                      type="number"
                      value={clusterParams.nPcs}
                      onChange={(e) => setClusterParams({ ...clusterParams, nPcs: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Principal components for clustering</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="resolution">Resolution</Label>
                    <Input
                      id="resolution"
                      type="number"
                      step="0.1"
                      value={clusterParams.resolution}
                      onChange={(e) => setClusterParams({ ...clusterParams, resolution: parseFloat(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Higher = more clusters</p>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="algorithm">Clustering Algorithm</Label>
                    <Select
                      value={clusterParams.algorithm}
                      onValueChange={(value: any) => setClusterParams({ ...clusterParams, algorithm: value })}
                    >
                      <SelectTrigger id="algorithm">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="louvain">Louvain (faster)</SelectItem>
                        <SelectItem value="leiden">Leiden (more accurate)</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="n-neighbors">Number of Neighbors</Label>
                    <Input
                      id="n-neighbors"
                      type="number"
                      value={clusterParams.nNeighbors}
                      onChange={(e) => setClusterParams({ ...clusterParams, nNeighbors: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">For k-NN graph construction</p>
                  </div>
                </div>

                <Button
                  onClick={handleClustering}
                  disabled={loading || !clusterParams.sampleId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Clustering...
                    </>
                  ) : (
                    'Run Clustering'
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Differential Expression Tab */}
        <TabsContent value="de" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Differential Expression Analysis</CardTitle>
              <CardDescription>
                Find marker genes that distinguish between cell clusters or conditions
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="de-sample">Sample ID</Label>
                  <Input
                    id="de-sample"
                    placeholder="Enter sample ID"
                    value={deParams.sampleId}
                    onChange={(e) => setDeParams({ ...deParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="cluster1">Cluster 1</Label>
                    <Input
                      id="cluster1"
                      placeholder="e.g., '0' or 'all'"
                      value={deParams.cluster1}
                      onChange={(e) => setDeParams({ ...deParams, cluster1: e.target.value })}
                    />
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="cluster2">Cluster 2 (optional)</Label>
                    <Input
                      id="cluster2"
                      placeholder="Leave empty for all others"
                      value={deParams.cluster2}
                      onChange={(e) => setDeParams({ ...deParams, cluster2: e.target.value })}
                    />
                  </div>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="de-method">Statistical Test</Label>
                  <Select
                    value={deParams.method}
                    onValueChange={(value: any) => setDeParams({ ...deParams, method: value })}
                  >
                    <SelectTrigger id="de-method">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="wilcoxon">Wilcoxon Rank-Sum (default)</SelectItem>
                      <SelectItem value="ttest">Student's t-test</SelectItem>
                      <SelectItem value="deseq2">DESeq2 (count-based)</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="min-pct">Min % Expressed</Label>
                    <Input
                      id="min-pct"
                      type="number"
                      step="0.05"
                      value={deParams.minPct}
                      onChange={(e) => setDeParams({ ...deParams, minPct: parseFloat(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Min fraction of cells expressing gene</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="logfc">Log Fold Change Threshold</Label>
                    <Input
                      id="logfc"
                      type="number"
                      step="0.1"
                      value={deParams.logfcThreshold}
                      onChange={(e) => setDeParams({ ...deParams, logfcThreshold: parseFloat(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Min absolute log2 fold change</p>
                  </div>
                </div>

                <Button
                  onClick={handleDifferentialExpression}
                  disabled={loading || !deParams.sampleId || !deParams.cluster1}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Analyzing...
                    </>
                  ) : (
                    'Run Differential Expression'
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Trajectory Analysis Tab */}
        <TabsContent value="trajectory" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Trajectory Analysis</CardTitle>
              <CardDescription>
                Infer cell developmental trajectories and pseudotime ordering
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="traj-sample">Sample ID</Label>
                  <Input
                    id="traj-sample"
                    placeholder="Enter sample ID"
                    value={trajectoryParams.sampleId}
                    onChange={(e) => setTrajectoryParams({ ...trajectoryParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="space-y-2">
                  <Label htmlFor="root-cell">Root Cell/Cluster</Label>
                  <Input
                    id="root-cell"
                    placeholder="Starting point for pseudotime (optional)"
                    value={trajectoryParams.rootCell}
                    onChange={(e) => setTrajectoryParams({ ...trajectoryParams, rootCell: e.target.value })}
                  />
                  <p className="text-xs text-muted-foreground">Specify initial state or leave blank for auto-detection</p>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="traj-method">Trajectory Method</Label>
                    <Select
                      value={trajectoryParams.method}
                      onValueChange={(value: any) => setTrajectoryParams({ ...trajectoryParams, method: value })}
                    >
                      <SelectTrigger id="traj-method">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="monocle3">Monocle 3 (UMAP-based)</SelectItem>
                        <SelectItem value="slingshot">Slingshot (smooth curves)</SelectItem>
                        <SelectItem value="paga">PAGA (graph abstraction)</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="n-components">Number of Components</Label>
                    <Input
                      id="n-components"
                      type="number"
                      value={trajectoryParams.nComponents}
                      onChange={(e) => setTrajectoryParams({ ...trajectoryParams, nComponents: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Dimensions for trajectory embedding</p>
                  </div>
                </div>

                <Button
                  onClick={handleTrajectoryAnalysis}
                  disabled={loading || !trajectoryParams.sampleId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Computing Trajectory...
                    </>
                  ) : (
                    'Run Trajectory Analysis'
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>

      {/* Error Display */}
      {error && (
        <Alert variant="destructive" className="mt-4">
          <AlertDescription>{error}</AlertDescription>
        </Alert>
      )}

      {/* Results Display */}
      {result && (
        <Card className="mt-4">
          <CardHeader>
            <CardTitle>Analysis Results</CardTitle>
          </CardHeader>
          <CardContent>
            <pre className="bg-muted p-4 rounded-lg overflow-auto max-h-96">
              {JSON.stringify(result, null, 2)}
            </pre>
          </CardContent>
        </Card>
      )}
    </div>
  );
};

export default SingleCellAnalysis;
