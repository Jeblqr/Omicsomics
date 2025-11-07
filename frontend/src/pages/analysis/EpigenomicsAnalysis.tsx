import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Checkbox } from '@/components/ui/checkbox';
import { Loader2, Dna, Target, Map } from 'lucide-react';

interface ChIPSeqParams {
  sampleId: string;
  controlId: string;
  peakCaller: 'macs2' | 'homer' | 'sicer';
  pvalueThreshold: number;
  foldEnrichment: number;
  genomeSize: string;
  fragmentLength: number;
}

interface ATACSeqParams {
  sampleId: string;
  peakCaller: 'macs2' | 'genrich';
  minFragmentLength: number;
  maxFragmentLength: number;
  tssEnrichment: boolean;
  nucleosomeFree: boolean;
}

interface MethylationParams {
  sampleId: string;
  analysisType: 'wgbs' | 'rrbs' | 'array';
  minCoverage: number;
  minSamples: number;
  diffMethod: 'dss' | 'methylkit' | 'dmrfinder';
  qvalueThreshold: number;
  methylDiffThreshold: number;
}

const EpigenomicsAnalysis: React.FC = () => {
  const [activeTab, setActiveTab] = useState<string>('chipseq');
  const [loading, setLoading] = useState<boolean>(false);
  const [result, setResult] = useState<any>(null);
  const [error, setError] = useState<string>('');

  // ChIP-seq State
  const [chipParams, setChipParams] = useState<ChIPSeqParams>({
    sampleId: '',
    controlId: '',
    peakCaller: 'macs2',
    pvalueThreshold: 0.05,
    foldEnrichment: 2,
    genomeSize: 'hs',
    fragmentLength: 200
  });

  // ATAC-seq State
  const [atacParams, setAtacParams] = useState<ATACSeqParams>({
    sampleId: '',
    peakCaller: 'macs2',
    minFragmentLength: 38,
    maxFragmentLength: 2000,
    tssEnrichment: true,
    nucleosomeFree: true
  });

  // DNA Methylation State
  const [methylParams, setMethylParams] = useState<MethylationParams>({
    sampleId: '',
    analysisType: 'wgbs',
    minCoverage: 10,
    minSamples: 3,
    diffMethod: 'dss',
    qvalueThreshold: 0.01,
    methylDiffThreshold: 0.25
  });

  const handleChIPSeq = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/epigenomics/chipseq', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(chipParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'ChIP-seq analysis failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during ChIP-seq analysis');
    } finally {
      setLoading(false);
    }
  };

  const handleATACSeq = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/epigenomics/atacseq', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(atacParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'ATAC-seq analysis failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during ATAC-seq analysis');
    } finally {
      setLoading(false);
    }
  };

  const handleMethylation = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/epigenomics/methylation', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(methylParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Methylation analysis failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during methylation analysis');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container mx-auto py-8 px-4">
      <div className="mb-8">
        <h1 className="text-4xl font-bold mb-2">Epigenomics Analysis</h1>
        <p className="text-muted-foreground">
          Analyze ChIP-seq, ATAC-seq, and DNA methylation data to study epigenetic modifications
        </p>
      </div>

      <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="chipseq" className="flex items-center gap-2">
            <Target className="h-4 w-4" />
            ChIP-seq
          </TabsTrigger>
          <TabsTrigger value="atacseq" className="flex items-center gap-2">
            <Map className="h-4 w-4" />
            ATAC-seq
          </TabsTrigger>
          <TabsTrigger value="methylation" className="flex items-center gap-2">
            <Dna className="h-4 w-4" />
            DNA Methylation
          </TabsTrigger>
        </TabsList>

        {/* ChIP-seq Tab */}
        <TabsContent value="chipseq" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>ChIP-seq Peak Calling</CardTitle>
              <CardDescription>
                Identify protein-DNA binding sites and histone modifications using ChIP-seq data
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="chip-sample">Treatment Sample ID</Label>
                    <Input
                      id="chip-sample"
                      placeholder="ChIP sample ID"
                      value={chipParams.sampleId}
                      onChange={(e) => setChipParams({ ...chipParams, sampleId: e.target.value })}
                    />
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="chip-control">Control/Input Sample ID</Label>
                    <Input
                      id="chip-control"
                      placeholder="Control sample ID"
                      value={chipParams.controlId}
                      onChange={(e) => setChipParams({ ...chipParams, controlId: e.target.value })}
                    />
                  </div>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="peak-caller">Peak Caller</Label>
                  <Select
                    value={chipParams.peakCaller}
                    onValueChange={(value: any) => setChipParams({ ...chipParams, peakCaller: value })}
                  >
                    <SelectTrigger id="peak-caller">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="macs2">MACS2 (Model-based Analysis)</SelectItem>
                      <SelectItem value="homer">HOMER</SelectItem>
                      <SelectItem value="sicer">SICER (Broad peaks)</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="pvalue">P-value Threshold</Label>
                    <Input
                      id="pvalue"
                      type="number"
                      step="0.001"
                      value={chipParams.pvalueThreshold}
                      onChange={(e) => setChipParams({ ...chipParams, pvalueThreshold: parseFloat(e.target.value) })}
                    />
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="fold-enrichment">Fold Enrichment</Label>
                    <Input
                      id="fold-enrichment"
                      type="number"
                      step="0.5"
                      value={chipParams.foldEnrichment}
                      onChange={(e) => setChipParams({ ...chipParams, foldEnrichment: parseFloat(e.target.value) })}
                    />
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="genome-size">Genome Size</Label>
                    <Select
                      value={chipParams.genomeSize}
                      onValueChange={(value) => setChipParams({ ...chipParams, genomeSize: value })}
                    >
                      <SelectTrigger id="genome-size">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="hs">Human (2.7e9)</SelectItem>
                        <SelectItem value="mm">Mouse (1.87e9)</SelectItem>
                        <SelectItem value="ce">C. elegans (9e7)</SelectItem>
                        <SelectItem value="dm">Drosophila (1.2e8)</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="fragment-length">Fragment Length (bp)</Label>
                    <Input
                      id="fragment-length"
                      type="number"
                      value={chipParams.fragmentLength}
                      onChange={(e) => setChipParams({ ...chipParams, fragmentLength: parseInt(e.target.value) })}
                    />
                  </div>
                </div>

                <Button
                  onClick={handleChIPSeq}
                  disabled={loading || !chipParams.sampleId || !chipParams.controlId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Running ChIP-seq Analysis...
                    </>
                  ) : (
                    'Run ChIP-seq Analysis'
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {/* ATAC-seq Tab */}
        <TabsContent value="atacseq" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>ATAC-seq Analysis</CardTitle>
              <CardDescription>
                Identify open chromatin regions and analyze chromatin accessibility
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="atac-sample">Sample ID</Label>
                  <Input
                    id="atac-sample"
                    placeholder="Enter ATAC-seq sample ID"
                    value={atacParams.sampleId}
                    onChange={(e) => setAtacParams({ ...atacParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="space-y-2">
                  <Label htmlFor="atac-peak-caller">Peak Caller</Label>
                  <Select
                    value={atacParams.peakCaller}
                    onValueChange={(value: any) => setAtacParams({ ...atacParams, peakCaller: value })}
                  >
                    <SelectTrigger id="atac-peak-caller">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="macs2">MACS2</SelectItem>
                      <SelectItem value="genrich">Genrich (ATAC-optimized)</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="min-frag">Min Fragment Length (bp)</Label>
                    <Input
                      id="min-frag"
                      type="number"
                      value={atacParams.minFragmentLength}
                      onChange={(e) => setAtacParams({ ...atacParams, minFragmentLength: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Nucleosome-free: &lt;100 bp</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="max-frag">Max Fragment Length (bp)</Label>
                    <Input
                      id="max-frag"
                      type="number"
                      value={atacParams.maxFragmentLength}
                      onChange={(e) => setAtacParams({ ...atacParams, maxFragmentLength: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Typical range: 38-2000</p>
                  </div>
                </div>

                <div className="space-y-4">
                  <div className="flex items-center space-x-2">
                    <Checkbox
                      id="tss-enrichment"
                      checked={atacParams.tssEnrichment}
                      onCheckedChange={(checked) => setAtacParams({ ...atacParams, tssEnrichment: checked as boolean })}
                    />
                    <label
                      htmlFor="tss-enrichment"
                      className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70"
                    >
                      Calculate TSS Enrichment Score
                    </label>
                  </div>
                  <p className="text-xs text-muted-foreground ml-6">
                    Quality metric: signal enrichment at transcription start sites
                  </p>

                  <div className="flex items-center space-x-2">
                    <Checkbox
                      id="nucleosome-free"
                      checked={atacParams.nucleosomeFree}
                      onCheckedChange={(checked) => setAtacParams({ ...atacParams, nucleosomeFree: checked as boolean })}
                    />
                    <label
                      htmlFor="nucleosome-free"
                      className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70"
                    >
                      Extract Nucleosome-Free Regions
                    </label>
                  </div>
                  <p className="text-xs text-muted-foreground ml-6">
                    Focus on fragments &lt;100 bp for transcription factor footprinting
                  </p>
                </div>

                <Button
                  onClick={handleATACSeq}
                  disabled={loading || !atacParams.sampleId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Running ATAC-seq Analysis...
                    </>
                  ) : (
                    'Run ATAC-seq Analysis'
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {/* DNA Methylation Tab */}
        <TabsContent value="methylation" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>DNA Methylation Analysis</CardTitle>
              <CardDescription>
                Analyze DNA methylation patterns and identify differentially methylated regions (DMRs)
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="methyl-sample">Sample ID</Label>
                  <Input
                    id="methyl-sample"
                    placeholder="Enter sample ID or project ID"
                    value={methylParams.sampleId}
                    onChange={(e) => setMethylParams({ ...methylParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="space-y-2">
                  <Label htmlFor="analysis-type">Analysis Type</Label>
                  <Select
                    value={methylParams.analysisType}
                    onValueChange={(value: any) => setMethylParams({ ...methylParams, analysisType: value })}
                  >
                    <SelectTrigger id="analysis-type">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="wgbs">WGBS (Whole Genome Bisulfite Seq)</SelectItem>
                      <SelectItem value="rrbs">RRBS (Reduced Representation)</SelectItem>
                      <SelectItem value="array">Methylation Array (450K/EPIC)</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="min-coverage">Min Coverage</Label>
                    <Input
                      id="min-coverage"
                      type="number"
                      value={methylParams.minCoverage}
                      onChange={(e) => setMethylParams({ ...methylParams, minCoverage: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Min reads per CpG site</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="min-samples">Min Samples per Group</Label>
                    <Input
                      id="min-samples"
                      type="number"
                      value={methylParams.minSamples}
                      onChange={(e) => setMethylParams({ ...methylParams, minSamples: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Statistical power threshold</p>
                  </div>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="diff-method">Differential Methylation Method</Label>
                  <Select
                    value={methylParams.diffMethod}
                    onValueChange={(value: any) => setMethylParams({ ...methylParams, diffMethod: value })}
                  >
                    <SelectTrigger id="diff-method">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="dss">DSS (Dispersion Shrinkage)</SelectItem>
                      <SelectItem value="methylkit">MethylKit</SelectItem>
                      <SelectItem value="dmrfinder">DMRfinder</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="qvalue">Q-value Threshold</Label>
                    <Input
                      id="qvalue"
                      type="number"
                      step="0.001"
                      value={methylParams.qvalueThreshold}
                      onChange={(e) => setMethylParams({ ...methylParams, qvalueThreshold: parseFloat(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">FDR-corrected significance</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="methyl-diff">Methylation Difference</Label>
                    <Input
                      id="methyl-diff"
                      type="number"
                      step="0.05"
                      value={methylParams.methylDiffThreshold}
                      onChange={(e) => setMethylParams({ ...methylParams, methylDiffThreshold: parseFloat(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Min % difference (0-1)</p>
                  </div>
                </div>

                <Button
                  onClick={handleMethylation}
                  disabled={loading || !methylParams.sampleId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Running Methylation Analysis...
                    </>
                  ) : (
                    'Run Methylation Analysis'
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

export default EpigenomicsAnalysis;
