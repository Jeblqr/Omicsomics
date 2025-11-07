import { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Badge } from '@/components/ui/badge';
import { Loader2, PlayCircle, CheckCircle, AlertCircle, Plus, X } from 'lucide-react';

interface WorkflowResult {
  workflow_id: number;
  status: string;
  message: string;
}

interface OmicsLayer {
  id: string;
  name: string;
  file: string;
}

export default function MultiOmicsAnalysis() {
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<WorkflowResult | null>(null);
  const [error, setError] = useState<string | null>(null);

  // MOFA2 state
  const [mofa2Layers, setMofa2Layers] = useState<OmicsLayer[]>([
    { id: '1', name: 'transcriptomics', file: '' },
    { id: '2', name: 'proteomics', file: '' },
  ]);
  const [mofa2Output, setMofa2Output] = useState('');
  const [nFactors, setNFactors] = useState(10);

  // DIABLO state
  const [diabloLayers, setDiabloLayers] = useState<OmicsLayer[]>([
    { id: '1', name: 'transcriptomics', file: '' },
    { id: '2', name: 'proteomics', file: '' },
  ]);
  const [phenotypeFile, setPhenotypeFile] = useState('');
  const [diabloOutput, setDiabloOutput] = useState('');
  const [nComponents, setNComponents] = useState(2);

  const addMofa2Layer = () => {
    const newId = (mofa2Layers.length + 1).toString();
    setMofa2Layers([...mofa2Layers, { id: newId, name: '', file: '' }]);
  };

  const removeMofa2Layer = (id: string) => {
    setMofa2Layers(mofa2Layers.filter((layer) => layer.id !== id));
  };

  const updateMofa2Layer = (id: string, field: 'name' | 'file', value: string) => {
    setMofa2Layers(
      mofa2Layers.map((layer) =>
        layer.id === id ? { ...layer, [field]: value } : layer
      )
    );
  };

  const addDiabloLayer = () => {
    const newId = (diabloLayers.length + 1).toString();
    setDiabloLayers([...diabloLayers, { id: newId, name: '', file: '' }]);
  };

  const removeDiabloLayer = (id: string) => {
    setDiabloLayers(diabloLayers.filter((layer) => layer.id !== id));
  };

  const updateDiabloLayer = (id: string, field: 'name' | 'file', value: string) => {
    setDiabloLayers(
      diabloLayers.map((layer) =>
        layer.id === id ? { ...layer, [field]: value } : layer
      )
    );
  };

  const runMOFA2 = async () => {
    setLoading(true);
    setError(null);
    try {
      const data_matrices: Record<string, string> = {};
      mofa2Layers.forEach((layer) => {
        if (layer.name && layer.file) {
          data_matrices[layer.name] = layer.file;
        }
      });

      const response = await fetch('/api/multiomics/mofa2', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify({
          sample_id: 1,
          data_matrices,
          output_dir: mofa2Output,
          n_factors: nFactors,
          convergence_mode: 'fast',
        }),
      });

      if (!response.ok) {
        throw new Error('MOFA2 request failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error');
    } finally {
      setLoading(false);
    }
  };

  const runDIABLO = async () => {
    setLoading(true);
    setError(null);
    try {
      const data_matrices: Record<string, string> = {};
      diabloLayers.forEach((layer) => {
        if (layer.name && layer.file) {
          data_matrices[layer.name] = layer.file;
        }
      });

      const response = await fetch('/api/multiomics/diablo', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify({
          sample_id: 1,
          data_matrices,
          phenotype_file: phenotypeFile,
          output_dir: diabloOutput,
          n_components: nComponents,
          design_correlation: 0.1,
        }),
      });

      if (!response.ok) {
        throw new Error('DIABLO request failed');
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
          <h1 className="text-3xl font-bold">Multi-Omics Integration</h1>
          <p className="text-muted-foreground">
            Integrate multiple omics layers with MOFA2 and DIABLO
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

      <Tabs defaultValue="mofa2" className="w-full">
        <TabsList className="grid w-full grid-cols-2">
          <TabsTrigger value="mofa2">MOFA2 (Unsupervised)</TabsTrigger>
          <TabsTrigger value="diablo">DIABLO (Supervised)</TabsTrigger>
        </TabsList>

        <TabsContent value="mofa2">
          <Card>
            <CardHeader>
              <CardTitle>MOFA2 - Multi-Omics Factor Analysis</CardTitle>
              <CardDescription>
                Unsupervised integration to discover latent factors across omics layers
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="flex items-center justify-between">
                  <Label>Omics Data Layers</Label>
                  <Button onClick={addMofa2Layer} variant="outline" size="sm">
                    <Plus className="h-4 w-4 mr-2" />
                    Add Layer
                  </Button>
                </div>

                {mofa2Layers.map((layer) => (
                  <div key={layer.id} className="flex items-end gap-4">
                    <div className="flex-1 space-y-2">
                      <Label>Omics Type</Label>
                      <Select
                        value={layer.name}
                        onValueChange={(value) => updateMofa2Layer(layer.id, 'name', value)}
                      >
                        <SelectTrigger>
                          <SelectValue placeholder="Select type..." />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="transcriptomics">Transcriptomics</SelectItem>
                          <SelectItem value="proteomics">Proteomics</SelectItem>
                          <SelectItem value="metabolomics">Metabolomics</SelectItem>
                          <SelectItem value="epigenomics">Epigenomics</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                    <div className="flex-[2] space-y-2">
                      <Label>Data Matrix File (CSV)</Label>
                      <Input
                        value={layer.file}
                        onChange={(e) => updateMofa2Layer(layer.id, 'file', e.target.value)}
                        placeholder="/data/rna_counts.csv"
                      />
                    </div>
                    <Button
                      onClick={() => removeMofa2Layer(layer.id)}
                      variant="ghost"
                      size="icon"
                      disabled={mofa2Layers.length <= 2}
                    >
                      <X className="h-4 w-4" />
                    </Button>
                  </div>
                ))}
              </div>

              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="mofa2_output">Output Directory</Label>
                  <Input
                    id="mofa2_output"
                    value={mofa2Output}
                    onChange={(e) => setMofa2Output(e.target.value)}
                    placeholder="/output/mofa2_results"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="n_factors">Number of Factors</Label>
                  <Input
                    id="n_factors"
                    type="number"
                    value={nFactors}
                    onChange={(e) => setNFactors(parseInt(e.target.value))}
                    min={1}
                    max={50}
                  />
                </div>
              </div>

              <div className="flex items-center gap-2 text-sm text-muted-foreground">
                <Badge variant="secondary">{mofa2Layers.length} layers</Badge>
                <span>MOFA2 will identify {nFactors} latent factors across omics</span>
              </div>

              <Button onClick={runMOFA2} disabled={loading} className="w-full">
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <PlayCircle className="mr-2 h-4 w-4" />
                    Run MOFA2 Integration
                  </>
                )}
              </Button>
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="diablo">
          <Card>
            <CardHeader>
              <CardTitle>DIABLO - Data Integration for Biomarker Discovery</CardTitle>
              <CardDescription>
                Supervised integration to identify multi-omics biomarkers for phenotype prediction
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="flex items-center justify-between">
                  <Label>Omics Data Layers</Label>
                  <Button onClick={addDiabloLayer} variant="outline" size="sm">
                    <Plus className="h-4 w-4 mr-2" />
                    Add Layer
                  </Button>
                </div>

                {diabloLayers.map((layer) => (
                  <div key={layer.id} className="flex items-end gap-4">
                    <div className="flex-1 space-y-2">
                      <Label>Omics Type</Label>
                      <Select
                        value={layer.name}
                        onValueChange={(value) => updateDiabloLayer(layer.id, 'name', value)}
                      >
                        <SelectTrigger>
                          <SelectValue placeholder="Select type..." />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="transcriptomics">Transcriptomics</SelectItem>
                          <SelectItem value="proteomics">Proteomics</SelectItem>
                          <SelectItem value="metabolomics">Metabolomics</SelectItem>
                          <SelectItem value="epigenomics">Epigenomics</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                    <div className="flex-[2] space-y-2">
                      <Label>Data Matrix File (CSV)</Label>
                      <Input
                        value={layer.file}
                        onChange={(e) => updateDiabloLayer(layer.id, 'file', e.target.value)}
                        placeholder="/data/rna_counts.csv"
                      />
                    </div>
                    <Button
                      onClick={() => removeDiabloLayer(layer.id)}
                      variant="ghost"
                      size="icon"
                      disabled={diabloLayers.length <= 2}
                    >
                      <X className="h-4 w-4" />
                    </Button>
                  </div>
                ))}
              </div>

              <div className="space-y-2">
                <Label htmlFor="phenotype">Phenotype File (CSV with class labels)</Label>
                <Input
                  id="phenotype"
                  value={phenotypeFile}
                  onChange={(e) => setPhenotypeFile(e.target.value)}
                  placeholder="/data/phenotype.csv"
                />
              </div>

              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="diablo_output">Output Directory</Label>
                  <Input
                    id="diablo_output"
                    value={diabloOutput}
                    onChange={(e) => setDiabloOutput(e.target.value)}
                    placeholder="/output/diablo_results"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="n_components">Number of Components</Label>
                  <Input
                    id="n_components"
                    type="number"
                    value={nComponents}
                    onChange={(e) => setNComponents(parseInt(e.target.value))}
                    min={1}
                    max={10}
                  />
                </div>
              </div>

              <div className="flex items-center gap-2 text-sm text-muted-foreground">
                <Badge variant="secondary">{diabloLayers.length} layers</Badge>
                <span>DIABLO will identify biomarkers across omics for phenotype classification</span>
              </div>

              <Button onClick={runDIABLO} disabled={loading} className="w-full">
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <PlayCircle className="mr-2 h-4 w-4" />
                    Run DIABLO Integration
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
