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

export default function GWASAnalysis() {
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<WorkflowResult | null>(null);
  const [error, setError] = useState<string | null>(null);

  // QC form state
  const [qcForm, setQcForm] = useState({
    sample_id: 1,
    bed_file: '',
    bim_file: '',
    fam_file: '',
    output_prefix: '',
    geno: 0.02,
    mind: 0.02,
    maf: 0.01,
    hwe: 0.000001,
  });

  // Association test form state
  const [assocForm, setAssocForm] = useState({
    sample_id: 1,
    bed_file: '',
    phenotype_file: '',
    output_prefix: '',
    covariates_file: '',
    binary_trait: false,
  });

  // MTAG form state
  const [mtagForm, setMtagForm] = useState({
    sample_id: 1,
    trait1_name: '',
    trait1_file: '',
    trait2_name: '',
    trait2_file: '',
    output_dir: '',
  });

  const runQC = async () => {
    setLoading(true);
    setError(null);
    try {
      const response = await fetch('/api/gwas/qc', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify(qcForm),
      });

      if (!response.ok) {
        throw new Error('QC request failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error');
    } finally {
      setLoading(false);
    }
  };

  const runAssociation = async () => {
    setLoading(true);
    setError(null);
    try {
      const response = await fetch('/api/gwas/association', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify(assocForm),
      });

      if (!response.ok) {
        throw new Error('Association test request failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error');
    } finally {
      setLoading(false);
    }
  };

  const runMTAG = async () => {
    setLoading(true);
    setError(null);
    try {
      const summary_stats_files: Record<string, string> = {};
      if (mtagForm.trait1_name && mtagForm.trait1_file) {
        summary_stats_files[mtagForm.trait1_name] = mtagForm.trait1_file;
      }
      if (mtagForm.trait2_name && mtagForm.trait2_file) {
        summary_stats_files[mtagForm.trait2_name] = mtagForm.trait2_file;
      }

      const response = await fetch('/api/gwas/mtag', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify({
          sample_id: mtagForm.sample_id,
          summary_stats_files,
          output_dir: mtagForm.output_dir,
        }),
      });

      if (!response.ok) {
        throw new Error('MTAG request failed');
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
          <h1 className="text-3xl font-bold">GWAS Analysis</h1>
          <p className="text-muted-foreground">
            Genome-Wide Association Studies and cross-trait analysis
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

      <Tabs defaultValue="qc" className="w-full">
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="qc">Quality Control</TabsTrigger>
          <TabsTrigger value="association">Association Test</TabsTrigger>
          <TabsTrigger value="mtag">MTAG (Cross-Trait)</TabsTrigger>
        </TabsList>

        <TabsContent value="qc">
          <Card>
            <CardHeader>
              <CardTitle>PLINK Quality Control</CardTitle>
              <CardDescription>
                Filter SNPs and samples based on quality metrics
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="bed_file">BED File</Label>
                  <Input
                    id="bed_file"
                    value={qcForm.bed_file}
                    onChange={(e) => setQcForm({ ...qcForm, bed_file: e.target.value })}
                    placeholder="/data/genotypes.bed"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="bim_file">BIM File</Label>
                  <Input
                    id="bim_file"
                    value={qcForm.bim_file}
                    onChange={(e) => setQcForm({ ...qcForm, bim_file: e.target.value })}
                    placeholder="/data/genotypes.bim"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="fam_file">FAM File</Label>
                  <Input
                    id="fam_file"
                    value={qcForm.fam_file}
                    onChange={(e) => setQcForm({ ...qcForm, fam_file: e.target.value })}
                    placeholder="/data/genotypes.fam"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="output_prefix">Output Prefix</Label>
                  <Input
                    id="output_prefix"
                    value={qcForm.output_prefix}
                    onChange={(e) => setQcForm({ ...qcForm, output_prefix: e.target.value })}
                    placeholder="/output/qc_filtered"
                  />
                </div>
              </div>

              <div className="grid grid-cols-4 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="geno">SNP Missing Rate (--geno)</Label>
                  <Input
                    id="geno"
                    type="number"
                    step="0.01"
                    value={qcForm.geno}
                    onChange={(e) => setQcForm({ ...qcForm, geno: parseFloat(e.target.value) })}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="mind">Individual Missing Rate (--mind)</Label>
                  <Input
                    id="mind"
                    type="number"
                    step="0.01"
                    value={qcForm.mind}
                    onChange={(e) => setQcForm({ ...qcForm, mind: parseFloat(e.target.value) })}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="maf">Minor Allele Frequency (--maf)</Label>
                  <Input
                    id="maf"
                    type="number"
                    step="0.01"
                    value={qcForm.maf}
                    onChange={(e) => setQcForm({ ...qcForm, maf: parseFloat(e.target.value) })}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="hwe">HWE p-value (--hwe)</Label>
                  <Input
                    id="hwe"
                    type="number"
                    step="0.000001"
                    value={qcForm.hwe}
                    onChange={(e) => setQcForm({ ...qcForm, hwe: parseFloat(e.target.value) })}
                  />
                </div>
              </div>

              <Button onClick={runQC} disabled={loading} className="w-full">
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <PlayCircle className="mr-2 h-4 w-4" />
                    Run Quality Control
                  </>
                )}
              </Button>
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="association">
          <Card>
            <CardHeader>
              <CardTitle>GWAS Association Test</CardTitle>
              <CardDescription>
                Test for genetic associations with phenotype
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="assoc_bed">BED File (QC-filtered)</Label>
                  <Input
                    id="assoc_bed"
                    value={assocForm.bed_file}
                    onChange={(e) => setAssocForm({ ...assocForm, bed_file: e.target.value })}
                    placeholder="/output/qc_filtered.bed"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="phenotype">Phenotype File</Label>
                  <Input
                    id="phenotype"
                    value={assocForm.phenotype_file}
                    onChange={(e) => setAssocForm({ ...assocForm, phenotype_file: e.target.value })}
                    placeholder="/data/phenotype.phen"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="assoc_output">Output Prefix</Label>
                  <Input
                    id="assoc_output"
                    value={assocForm.output_prefix}
                    onChange={(e) => setAssocForm({ ...assocForm, output_prefix: e.target.value })}
                    placeholder="/output/gwas_results"
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="covariates">Covariates File (optional)</Label>
                  <Input
                    id="covariates"
                    value={assocForm.covariates_file}
                    onChange={(e) => setAssocForm({ ...assocForm, covariates_file: e.target.value })}
                    placeholder="/data/covariates.txt"
                  />
                </div>
              </div>

              <div className="space-y-2">
                <Label htmlFor="trait_type">Trait Type</Label>
                <Select
                  value={assocForm.binary_trait ? 'binary' : 'quantitative'}
                  onValueChange={(value) =>
                    setAssocForm({ ...assocForm, binary_trait: value === 'binary' })
                  }
                >
                  <SelectTrigger>
                    <SelectValue />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="quantitative">Quantitative (Linear)</SelectItem>
                    <SelectItem value="binary">Binary (Logistic)</SelectItem>
                  </SelectContent>
                </Select>
              </div>

              <Button onClick={runAssociation} disabled={loading} className="w-full">
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <PlayCircle className="mr-2 h-4 w-4" />
                    Run Association Test
                  </>
                )}
              </Button>
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="mtag">
          <Card>
            <CardHeader>
              <CardTitle>MTAG (Multi-Trait Analysis of GWAS)</CardTitle>
              <CardDescription>
                Cross-trait meta-analysis to boost power and identify shared genetics
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-4">
                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="trait1_name">Trait 1 Name</Label>
                    <Input
                      id="trait1_name"
                      value={mtagForm.trait1_name}
                      onChange={(e) => setMtagForm({ ...mtagForm, trait1_name: e.target.value })}
                      placeholder="height"
                    />
                  </div>
                  <div className="space-y-2">
                    <Label htmlFor="trait1_file">Trait 1 Summary Stats</Label>
                    <Input
                      id="trait1_file"
                      value={mtagForm.trait1_file}
                      onChange={(e) => setMtagForm({ ...mtagForm, trait1_file: e.target.value })}
                      placeholder="/data/height_sumstats.txt"
                    />
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="trait2_name">Trait 2 Name</Label>
                    <Input
                      id="trait2_name"
                      value={mtagForm.trait2_name}
                      onChange={(e) => setMtagForm({ ...mtagForm, trait2_name: e.target.value })}
                      placeholder="bmi"
                    />
                  </div>
                  <div className="space-y-2">
                    <Label htmlFor="trait2_file">Trait 2 Summary Stats</Label>
                    <Input
                      id="trait2_file"
                      value={mtagForm.trait2_file}
                      onChange={(e) => setMtagForm({ ...mtagForm, trait2_file: e.target.value })}
                      placeholder="/data/bmi_sumstats.txt"
                    />
                  </div>
                </div>

                <div className="space-y-2">
                  <Label htmlFor="mtag_output">Output Directory</Label>
                  <Input
                    id="mtag_output"
                    value={mtagForm.output_dir}
                    onChange={(e) => setMtagForm({ ...mtagForm, output_dir: e.target.value })}
                    placeholder="/output/mtag_results"
                  />
                </div>
              </div>

              <Button onClick={runMTAG} disabled={loading} className="w-full">
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <PlayCircle className="mr-2 h-4 w-4" />
                    Run MTAG Analysis
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
