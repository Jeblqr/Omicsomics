import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Checkbox } from '@/components/ui/checkbox';
import { Loader2, FlaskConical, Activity, TrendingUp } from 'lucide-react';

interface SearchParams {
  sampleId: string;
  searchEngine: 'mascot' | 'sequest' | 'maxquant';
  database: string;
  precursorTolerance: number;
  fragmentTolerance: number;
  enzyme: string;
  missedCleavages: number;
  fixedMods: string[];
  variableMods: string[];
}

interface QuantParams {
  sampleId: string;
  method: 'label-free' | 'itraq' | 'tmt' | 'silac';
  normalization: string;
  imputationMethod: string;
}

interface DiffParams {
  control: string;
  treatment: string;
  method: 'ttest' | 'limma' | 'deseq2';
  pvalThreshold: number;
  foldChangeThreshold: number;
  adjustMethod: string;
}

const ProteomicsAnalysis: React.FC = () => {
  const [activeTab, setActiveTab] = useState<string>('search');
  const [loading, setLoading] = useState<boolean>(false);
  const [result, setResult] = useState<any>(null);
  const [error, setError] = useState<string>('');

  // Peptide Search State
  const [searchParams, setSearchParams] = useState<SearchParams>({
    sampleId: '',
    searchEngine: 'mascot',
    database: 'uniprot_human',
    precursorTolerance: 10,
    fragmentTolerance: 0.5,
    enzyme: 'trypsin',
    missedCleavages: 2,
    fixedMods: ['carbamidomethyl_c'],
    variableMods: ['oxidation_m']
  });

  // Quantification State
  const [quantParams, setQuantParams] = useState<QuantParams>({
    sampleId: '',
    method: 'label-free',
    normalization: 'median',
    imputationMethod: 'knn'
  });

  // Differential Expression State
  const [diffParams, setDiffParams] = useState<DiffParams>({
    control: '',
    treatment: '',
    method: 'limma',
    pvalThreshold: 0.05,
    foldChangeThreshold: 1.5,
    adjustMethod: 'BH'
  });

  // Available modifications
  const fixedModOptions = [
    { value: 'carbamidomethyl_c', label: 'Carbamidomethyl (C)' },
    { value: 'propionamide_c', label: 'Propionamide (C)' },
    { value: 'pyro-glu_n-term_q', label: 'Pyro-glu (N-term Q)' },
    { value: 'acetyl_n-term', label: 'Acetyl (N-term)' }
  ];

  const variableModOptions = [
    { value: 'oxidation_m', label: 'Oxidation (M)' },
    { value: 'phosphorylation_sty', label: 'Phosphorylation (STY)' },
    { value: 'acetylation_k', label: 'Acetylation (K)' },
    { value: 'deamidation_nq', label: 'Deamidation (NQ)' },
    { value: 'methylation_k', label: 'Methylation (K)' }
  ];

  const handlePeptideSearch = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/proteomics/peptide-search', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(searchParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Peptide search failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during peptide search');
    } finally {
      setLoading(false);
    }
  };

  const handleQuantification = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/proteomics/quantification', {
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

  const handleDifferentialAnalysis = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/proteomics/differential', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(diffParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Differential analysis failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during differential analysis');
    } finally {
      setLoading(false);
    }
  };

  const toggleFixedMod = (mod: string) => {
    setSearchParams(prev => ({
      ...prev,
      fixedMods: prev.fixedMods.includes(mod)
        ? prev.fixedMods.filter(m => m !== mod)
        : [...prev.fixedMods, mod]
    }));
  };

  const toggleVariableMod = (mod: string) => {
    setSearchParams(prev => ({
      ...prev,
      variableMods: prev.variableMods.includes(mod)
        ? prev.variableMods.filter(m => m !== mod)
        : [...prev.variableMods, mod]
    }));
  };

  return (
    <div className="container mx-auto py-8 px-4">
      <div className="mb-8">
        <h1 className="text-4xl font-bold mb-2">Proteomics Analysis</h1>
        <p className="text-muted-foreground">
          Perform peptide identification, protein quantification, and differential expression analysis
        </p>
      </div>

      <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="search" className="flex items-center gap-2">
            <FlaskConical className="h-4 w-4" />
            Peptide Search
          </TabsTrigger>
          <TabsTrigger value="quantification" className="flex items-center gap-2">
            <Activity className="h-4 w-4" />
            Quantification
          </TabsTrigger>
          <TabsTrigger value="differential" className="flex items-center gap-2">
            <TrendingUp className="h-4 w-4" />
            Differential Expression
          </TabsTrigger>
        </TabsList>

        {/* Peptide Search Tab */}
        <TabsContent value="search" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Peptide Identification</CardTitle>
              <CardDescription>
                Identify peptides from MS/MS spectra using search engines like Mascot, SEQUEST, or MaxQuant
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="sample-id-search">Sample ID</Label>
                  <Input
                    id="sample-id-search"
                    placeholder="Enter sample ID"
                    value={searchParams.sampleId}
                    onChange={(e) => setSearchParams({ ...searchParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="search-engine">Search Engine</Label>
                    <Select
                      value={searchParams.searchEngine}
                      onValueChange={(value: any) => setSearchParams({ ...searchParams, searchEngine: value })}
                    >
                      <SelectTrigger id="search-engine">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="mascot">Mascot</SelectItem>
                        <SelectItem value="sequest">SEQUEST</SelectItem>
                        <SelectItem value="maxquant">MaxQuant</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="database">Database</Label>
                    <Select
                      value={searchParams.database}
                      onValueChange={(value) => setSearchParams({ ...searchParams, database: value })}
                    >
                      <SelectTrigger id="database">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="uniprot_human">UniProt Human</SelectItem>
                        <SelectItem value="uniprot_mouse">UniProt Mouse</SelectItem>
                        <SelectItem value="uniprot_rat">UniProt Rat</SelectItem>
                        <SelectItem value="refseq">RefSeq</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="precursor-tol">Precursor Tolerance (ppm)</Label>
                    <Input
                      id="precursor-tol"
                      type="number"
                      value={searchParams.precursorTolerance}
                      onChange={(e) => setSearchParams({ ...searchParams, precursorTolerance: parseFloat(e.target.value) })}
                    />
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="fragment-tol">Fragment Tolerance (Da)</Label>
                    <Input
                      id="fragment-tol"
                      type="number"
                      step="0.1"
                      value={searchParams.fragmentTolerance}
                      onChange={(e) => setSearchParams({ ...searchParams, fragmentTolerance: parseFloat(e.target.value) })}
                    />
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="enzyme">Enzyme</Label>
                    <Select
                      value={searchParams.enzyme}
                      onValueChange={(value) => setSearchParams({ ...searchParams, enzyme: value })}
                    >
                      <SelectTrigger id="enzyme">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="trypsin">Trypsin</SelectItem>
                        <SelectItem value="chymotrypsin">Chymotrypsin</SelectItem>
                        <SelectItem value="lysc">LysC</SelectItem>
                        <SelectItem value="gluc">GluC</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="missed-cleavages">Missed Cleavages</Label>
                    <Input
                      id="missed-cleavages"
                      type="number"
                      value={searchParams.missedCleavages}
                      onChange={(e) => setSearchParams({ ...searchParams, missedCleavages: parseInt(e.target.value) })}
                    />
                  </div>
                </div>

                <div className="space-y-4">
                  <div className="space-y-2">
                    <Label>Fixed Modifications</Label>
                    <div className="grid grid-cols-2 gap-2">
                      {fixedModOptions.map((mod) => (
                        <div key={mod.value} className="flex items-center space-x-2">
                          <Checkbox
                            id={`fixed-${mod.value}`}
                            checked={searchParams.fixedMods.includes(mod.value)}
                            onCheckedChange={() => toggleFixedMod(mod.value)}
                          />
                          <label
                            htmlFor={`fixed-${mod.value}`}
                            className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70"
                          >
                            {mod.label}
                          </label>
                        </div>
                      ))}
                    </div>
                  </div>

                  <div className="space-y-2">
                    <Label>Variable Modifications</Label>
                    <div className="grid grid-cols-2 gap-2">
                      {variableModOptions.map((mod) => (
                        <div key={mod.value} className="flex items-center space-x-2">
                          <Checkbox
                            id={`var-${mod.value}`}
                            checked={searchParams.variableMods.includes(mod.value)}
                            onCheckedChange={() => toggleVariableMod(mod.value)}
                          />
                          <label
                            htmlFor={`var-${mod.value}`}
                            className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70"
                          >
                            {mod.label}
                          </label>
                        </div>
                      ))}
                    </div>
                  </div>
                </div>

                <Button
                  onClick={handlePeptideSearch}
                  disabled={loading || !searchParams.sampleId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Searching...
                    </>
                  ) : (
                    'Run Peptide Search'
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Quantification Tab */}
        <TabsContent value="quantification" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Protein Quantification</CardTitle>
              <CardDescription>
                Quantify protein abundance using label-free or labeling methods (iTRAQ, TMT, SILAC)
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="sample-id-quant">Sample ID</Label>
                  <Input
                    id="sample-id-quant"
                    placeholder="Enter sample ID"
                    value={quantParams.sampleId}
                    onChange={(e) => setQuantParams({ ...quantParams, sampleId: e.target.value })}
                  />
                </div>

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
                      <SelectItem value="label-free">Label-Free (LFQ)</SelectItem>
                      <SelectItem value="itraq">iTRAQ</SelectItem>
                      <SelectItem value="tmt">TMT (Tandem Mass Tag)</SelectItem>
                      <SelectItem value="silac">SILAC</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="normalization">Normalization Method</Label>
                  <Select
                    value={quantParams.normalization}
                    onValueChange={(value) => setQuantParams({ ...quantParams, normalization: value })}
                  >
                    <SelectTrigger id="normalization">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="median">Median</SelectItem>
                      <SelectItem value="mean">Mean</SelectItem>
                      <SelectItem value="quantile">Quantile</SelectItem>
                      <SelectItem value="vsn">VSN (Variance Stabilization)</SelectItem>
                      <SelectItem value="loess">LOESS</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="imputation">Imputation Method</Label>
                  <Select
                    value={quantParams.imputationMethod}
                    onValueChange={(value) => setQuantParams({ ...quantParams, imputationMethod: value })}
                  >
                    <SelectTrigger id="imputation">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="knn">KNN</SelectItem>
                      <SelectItem value="min">Min Value</SelectItem>
                      <SelectItem value="zero">Zero</SelectItem>
                      <SelectItem value="mle">Maximum Likelihood</SelectItem>
                    </SelectContent>
                  </Select>
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
                Identify differentially expressed proteins between conditions using statistical methods
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="control-group">Control Group</Label>
                    <Input
                      id="control-group"
                      placeholder="Control sample IDs"
                      value={diffParams.control}
                      onChange={(e) => setDiffParams({ ...diffParams, control: e.target.value })}
                    />
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="treatment-group">Treatment Group</Label>
                    <Input
                      id="treatment-group"
                      placeholder="Treatment sample IDs"
                      value={diffParams.treatment}
                      onChange={(e) => setDiffParams({ ...diffParams, treatment: e.target.value })}
                    />
                  </div>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="diff-method">Statistical Method</Label>
                  <Select
                    value={diffParams.method}
                    onValueChange={(value: any) => setDiffParams({ ...diffParams, method: value })}
                  >
                    <SelectTrigger id="diff-method">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="ttest">T-test</SelectItem>
                      <SelectItem value="limma">Limma</SelectItem>
                      <SelectItem value="deseq2">DESeq2</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="pval-threshold">P-value Threshold</Label>
                    <Input
                      id="pval-threshold"
                      type="number"
                      step="0.001"
                      value={diffParams.pvalThreshold}
                      onChange={(e) => setDiffParams({ ...diffParams, pvalThreshold: parseFloat(e.target.value) })}
                    />
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="fc-threshold">Fold Change Threshold</Label>
                    <Input
                      id="fc-threshold"
                      type="number"
                      step="0.1"
                      value={diffParams.foldChangeThreshold}
                      onChange={(e) => setDiffParams({ ...diffParams, foldChangeThreshold: parseFloat(e.target.value) })}
                    />
                  </div>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="adjust-method">P-value Adjustment</Label>
                  <Select
                    value={diffParams.adjustMethod}
                    onValueChange={(value) => setDiffParams({ ...diffParams, adjustMethod: value })}
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
                  onClick={handleDifferentialAnalysis}
                  disabled={loading || !diffParams.control || !diffParams.treatment}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Analyzing...
                    </>
                  ) : (
                    'Run Differential Analysis'
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

export default ProteomicsAnalysis;
