import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Loader2, BarChart4, TrendingUp, Network } from 'lucide-react';

interface QuantParams {
  sampleId: string;
  method: 'salmon' | 'kallisto' | 'featurecounts';
  reference: string;
  strandedness: 'unstranded' | 'forward' | 'reverse';
  bootstraps: number;
}

interface DEParams {
  controlGroup: string;
  treatmentGroup: string;
  method: 'deseq2' | 'edger' | 'limma';
  pvalThreshold: number;
  logfcThreshold: number;
  adjustMethod: string;
}

interface PathwayParams {
  deResultsId: string;
  database: 'kegg' | 'go' | 'reactome';
  organism: string;
  pvalCutoff: number;
  qvalCutoff: number;
}

const TranscriptomicsAnalysis: React.FC = () => {
  const [activeTab, setActiveTab] = useState<string>('quantification');
  const [loading, setLoading] = useState<boolean>(false);
  const [result, setResult] = useState<any>(null);
  const [error, setError] = useState<string>('');

  // Quantification State
  const [quantParams, setQuantParams] = useState<QuantParams>({
    sampleId: '',
    method: 'salmon',
    reference: 'gencode_v44',
    strandedness: 'unstranded',
    bootstraps: 100
  });

  // Differential Expression State
  const [deParams, setDeParams] = useState<DEParams>({
    controlGroup: '',
    treatmentGroup: '',
    method: 'deseq2',
    pvalThreshold: 0.05,
    logfcThreshold: 1.0,
    adjustMethod: 'BH'
  });

  // Pathway Enrichment State
  const [pathwayParams, setPathwayParams] = useState<PathwayParams>({
    deResultsId: '',
    database: 'kegg',
    organism: 'hsa',
    pvalCutoff: 0.05,
    qvalCutoff: 0.2
  });

  const handleQuantification = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/transcriptomics/quantification', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(quantParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Quantification failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during quantification');
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
      const response = await fetch('http://localhost:8000/api/v1/transcriptomics/differential', {
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

  const handlePathwayEnrichment = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/transcriptomics/enrichment', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(pathwayParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Pathway enrichment failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during pathway enrichment');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container mx-auto py-8 px-4">
      <div className="mb-8">
        <h1 className="text-4xl font-bold mb-2">Transcriptomics Analysis</h1>
        <p className="text-muted-foreground">
          Perform transcript quantification, differential expression, and pathway enrichment analysis on RNA-seq data
        </p>
      </div>

      <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="quantification" className="flex items-center gap-2">
            <BarChart4 className="h-4 w-4" />
            Quantification
          </TabsTrigger>
          <TabsTrigger value="differential" className="flex items-center gap-2">
            <TrendingUp className="h-4 w-4" />
            Differential Expression
          </TabsTrigger>
          <TabsTrigger value="enrichment" className="flex items-center gap-2">
            <Network className="h-4 w-4" />
            Pathway Enrichment
          </TabsTrigger>
        </TabsList>

        {/* Quantification Tab */}
        <TabsContent value="quantification" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Transcript Quantification</CardTitle>
              <CardDescription>
                Estimate gene/transcript expression levels using Salmon, Kallisto, or featureCounts
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="quant-sample">Sample ID</Label>
                  <Input
                    id="quant-sample"
                    placeholder="Enter sample ID or batch ID"
                    value={quantParams.sampleId}
                    onChange={(e) => setQuantParams({ ...quantParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="quant-method">Quantification Method</Label>
                    <Select
                      value={quantParams.method}
                      onValueChange={(value: any) => setQuantParams({ ...quantParams, method: value })}
                    >
                      <SelectTrigger id="quant-method">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="salmon">Salmon (fast, accurate)</SelectItem>
                        <SelectItem value="kallisto">Kallisto (ultra-fast)</SelectItem>
                        <SelectItem value="featurecounts">featureCounts (alignment-based)</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="reference">Reference Transcriptome</Label>
                    <Select
                      value={quantParams.reference}
                      onValueChange={(value) => setQuantParams({ ...quantParams, reference: value })}
                    >
                      <SelectTrigger id="reference">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="gencode_v44">GENCODE v44 (Human)</SelectItem>
                        <SelectItem value="gencode_m32">GENCODE M32 (Mouse)</SelectItem>
                        <SelectItem value="ensembl_110">Ensembl 110</SelectItem>
                        <SelectItem value="refseq">RefSeq</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="strandedness">Library Strandedness</Label>
                    <Select
                      value={quantParams.strandedness}
                      onValueChange={(value: any) => setQuantParams({ ...quantParams, strandedness: value })}
                    >
                      <SelectTrigger id="strandedness">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="unstranded">Unstranded</SelectItem>
                        <SelectItem value="forward">Forward (RF/fr-firststrand)</SelectItem>
                        <SelectItem value="reverse">Reverse (FR/fr-secondstrand)</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="bootstraps">Bootstrap Samples</Label>
                    <Input
                      id="bootstraps"
                      type="number"
                      value={quantParams.bootstraps}
                      onChange={(e) => setQuantParams({ ...quantParams, bootstraps: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">For uncertainty estimation (Salmon/Kallisto)</p>
                  </div>
                </div>

                <Button
                  onClick={handleQuantification}
                  disabled={loading || !quantParams.sampleId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Quantifying...
                    </>
                  ) : (
                    'Run Quantification'
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Differential Expression Tab */}
        <TabsContent value="differential" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Differential Expression Analysis</CardTitle>
              <CardDescription>
                Identify differentially expressed genes between conditions using DESeq2, edgeR, or limma-voom
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="control">Control Group</Label>
                    <Input
                      id="control"
                      placeholder="Control sample IDs (comma-separated)"
                      value={deParams.controlGroup}
                      onChange={(e) => setDeParams({ ...deParams, controlGroup: e.target.value })}
                    />
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="treatment">Treatment Group</Label>
                    <Input
                      id="treatment"
                      placeholder="Treatment sample IDs (comma-separated)"
                      value={deParams.treatmentGroup}
                      onChange={(e) => setDeParams({ ...deParams, treatmentGroup: e.target.value })}
                    />
                  </div>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="de-method">Statistical Method</Label>
                  <Select
                    value={deParams.method}
                    onValueChange={(value: any) => setDeParams({ ...deParams, method: value })}
                  >
                    <SelectTrigger id="de-method">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="deseq2">DESeq2 (negative binomial)</SelectItem>
                      <SelectItem value="edger">edgeR (exact test)</SelectItem>
                      <SelectItem value="limma">limma-voom (microarray/RNA-seq)</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="de-pval">P-value Threshold</Label>
                    <Input
                      id="de-pval"
                      type="number"
                      step="0.001"
                      value={deParams.pvalThreshold}
                      onChange={(e) => setDeParams({ ...deParams, pvalThreshold: parseFloat(e.target.value) })}
                    />
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="de-logfc">Log2 Fold Change Threshold</Label>
                    <Input
                      id="de-logfc"
                      type="number"
                      step="0.1"
                      value={deParams.logfcThreshold}
                      onChange={(e) => setDeParams({ ...deParams, logfcThreshold: parseFloat(e.target.value) })}
                    />
                  </div>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="adjust-method">P-value Adjustment</Label>
                  <Select
                    value={deParams.adjustMethod}
                    onValueChange={(value) => setDeParams({ ...deParams, adjustMethod: value })}
                  >
                    <SelectTrigger id="adjust-method">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="BH">Benjamini-Hochberg (FDR)</SelectItem>
                      <SelectItem value="bonferroni">Bonferroni</SelectItem>
                      <SelectItem value="holm">Holm</SelectItem>
                      <SelectItem value="BY">Benjamini-Yekutieli</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <Button
                  onClick={handleDifferentialExpression}
                  disabled={loading || !deParams.controlGroup || !deParams.treatmentGroup}
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

        {/* Pathway Enrichment Tab */}
        <TabsContent value="enrichment" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Pathway Enrichment Analysis</CardTitle>
              <CardDescription>
                Perform gene set enrichment analysis (GSEA) using KEGG, GO, or Reactome databases
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="de-results">DE Results ID</Label>
                  <Input
                    id="de-results"
                    placeholder="Differential expression results ID"
                    value={pathwayParams.deResultsId}
                    onChange={(e) => setPathwayParams({ ...pathwayParams, deResultsId: e.target.value })}
                  />
                  <p className="text-xs text-muted-foreground">Use results from differential expression analysis</p>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="database">Pathway Database</Label>
                    <Select
                      value={pathwayParams.database}
                      onValueChange={(value: any) => setPathwayParams({ ...pathwayParams, database: value })}
                    >
                      <SelectTrigger id="database">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="kegg">KEGG Pathway</SelectItem>
                        <SelectItem value="go">Gene Ontology (GO)</SelectItem>
                        <SelectItem value="reactome">Reactome</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="organism">Organism</Label>
                    <Select
                      value={pathwayParams.organism}
                      onValueChange={(value) => setPathwayParams({ ...pathwayParams, organism: value })}
                    >
                      <SelectTrigger id="organism">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="hsa">Human (hsa)</SelectItem>
                        <SelectItem value="mmu">Mouse (mmu)</SelectItem>
                        <SelectItem value="rno">Rat (rno)</SelectItem>
                        <SelectItem value="dme">Drosophila (dme)</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="pval-cutoff">P-value Cutoff</Label>
                    <Input
                      id="pval-cutoff"
                      type="number"
                      step="0.001"
                      value={pathwayParams.pvalCutoff}
                      onChange={(e) => setPathwayParams({ ...pathwayParams, pvalCutoff: parseFloat(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Significance threshold</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="qval-cutoff">Q-value Cutoff</Label>
                    <Input
                      id="qval-cutoff"
                      type="number"
                      step="0.01"
                      value={pathwayParams.qvalCutoff}
                      onChange={(e) => setPathwayParams({ ...pathwayParams, qvalCutoff: parseFloat(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">FDR threshold</p>
                  </div>
                </div>

                <Button
                  onClick={handlePathwayEnrichment}
                  disabled={loading || !pathwayParams.deResultsId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Running Enrichment...
                    </>
                  ) : (
                    'Run Pathway Enrichment'
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

export default TranscriptomicsAnalysis;
