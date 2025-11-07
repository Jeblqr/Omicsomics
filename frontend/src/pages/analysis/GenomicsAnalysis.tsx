import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Checkbox } from '@/components/ui/checkbox';
import { Loader2, Dna, Search, BarChart3 } from 'lucide-react';

interface QCParams {
  sampleId: string;
  minQuality: number;
  maxNs: number;
  trimAdapters: boolean;
  removeDuplicates: boolean;
}

interface AlignmentParams {
  sampleId: string;
  reference: string;
  aligner: 'bwa' | 'bowtie2' | 'star';
  threads: number;
  markDuplicates: boolean;
}

interface VariantParams {
  sampleId: string;
  caller: 'gatk' | 'freebayes' | 'bcftools';
  minDepth: number;
  minQuality: number;
  filterSNPs: boolean;
  filterIndels: boolean;
}

const GenomicsAnalysis: React.FC = () => {
  const [activeTab, setActiveTab] = useState<string>('qc');
  const [loading, setLoading] = useState<boolean>(false);
  const [result, setResult] = useState<any>(null);
  const [error, setError] = useState<string>('');

  // Quality Control State
  const [qcParams, setQcParams] = useState<QCParams>({
    sampleId: '',
    minQuality: 30,
    maxNs: 5,
    trimAdapters: true,
    removeDuplicates: true
  });

  // Alignment State
  const [alignParams, setAlignParams] = useState<AlignmentParams>({
    sampleId: '',
    reference: 'hg38',
    aligner: 'bwa',
    threads: 8,
    markDuplicates: true
  });

  // Variant Calling State
  const [variantParams, setVariantParams] = useState<VariantParams>({
    sampleId: '',
    caller: 'gatk',
    minDepth: 10,
    minQuality: 30,
    filterSNPs: true,
    filterIndels: true
  });

  const handleQualityControl = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/genomics/qc', {
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

  const handleAlignment = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/genomics/alignment', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(alignParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Alignment failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during alignment');
    } finally {
      setLoading(false);
    }
  };

  const handleVariantCalling = async () => {
    setLoading(true);
    setError('');
    setResult(null);

    try {
      const token = localStorage.getItem('token');
      const response = await fetch('http://localhost:8000/api/v1/genomics/variant-calling', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(variantParams)
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Variant calling failed');
      }

      const data = await response.json();
      setResult(data);
    } catch (err: any) {
      setError(err.message || 'An error occurred during variant calling');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container mx-auto py-8 px-4">
      <div className="mb-8">
        <h1 className="text-4xl font-bold mb-2">Genomics Analysis</h1>
        <p className="text-muted-foreground">
          Perform quality control, read alignment, and variant calling on whole genome/exome sequencing data
        </p>
      </div>

      <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="qc" className="flex items-center gap-2">
            <BarChart3 className="h-4 w-4" />
            Quality Control
          </TabsTrigger>
          <TabsTrigger value="alignment" className="flex items-center gap-2">
            <Dna className="h-4 w-4" />
            Alignment
          </TabsTrigger>
          <TabsTrigger value="variants" className="flex items-center gap-2">
            <Search className="h-4 w-4" />
            Variant Calling
          </TabsTrigger>
        </TabsList>

        {/* Quality Control Tab */}
        <TabsContent value="qc" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Sequencing Quality Control</CardTitle>
              <CardDescription>
                Assess read quality, trim adapters, and remove low-quality sequences using FastQC and Trimmomatic
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
                    <Label htmlFor="min-quality">Min Quality Score (Phred)</Label>
                    <Input
                      id="min-quality"
                      type="number"
                      value={qcParams.minQuality}
                      onChange={(e) => setQcParams({ ...qcParams, minQuality: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Typical: 20-30</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="max-ns">Max N Bases</Label>
                    <Input
                      id="max-ns"
                      type="number"
                      value={qcParams.maxNs}
                      onChange={(e) => setQcParams({ ...qcParams, maxNs: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Ambiguous base threshold</p>
                  </div>
                </div>

                <div className="space-y-4">
                  <div className="flex items-center space-x-2">
                    <Checkbox
                      id="trim-adapters"
                      checked={qcParams.trimAdapters}
                      onCheckedChange={(checked) => setQcParams({ ...qcParams, trimAdapters: checked as boolean })}
                    />
                    <label
                      htmlFor="trim-adapters"
                      className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70"
                    >
                      Trim Illumina Adapters
                    </label>
                  </div>

                  <div className="flex items-center space-x-2">
                    <Checkbox
                      id="remove-duplicates"
                      checked={qcParams.removeDuplicates}
                      onCheckedChange={(checked) => setQcParams({ ...qcParams, removeDuplicates: checked as boolean })}
                    />
                    <label
                      htmlFor="remove-duplicates"
                      className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70"
                    >
                      Remove PCR Duplicates
                    </label>
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

        {/* Alignment Tab */}
        <TabsContent value="alignment" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Read Alignment</CardTitle>
              <CardDescription>
                Align sequencing reads to reference genome using BWA, Bowtie2, or STAR
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="align-sample">Sample ID</Label>
                  <Input
                    id="align-sample"
                    placeholder="Enter sample ID"
                    value={alignParams.sampleId}
                    onChange={(e) => setAlignParams({ ...alignParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="reference">Reference Genome</Label>
                    <Select
                      value={alignParams.reference}
                      onValueChange={(value) => setAlignParams({ ...alignParams, reference: value })}
                    >
                      <SelectTrigger id="reference">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="hg38">Human GRCh38/hg38</SelectItem>
                        <SelectItem value="hg19">Human GRCh37/hg19</SelectItem>
                        <SelectItem value="mm10">Mouse GRCm38/mm10</SelectItem>
                        <SelectItem value="mm39">Mouse GRCm39/mm39</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="aligner">Alignment Tool</Label>
                    <Select
                      value={alignParams.aligner}
                      onValueChange={(value: any) => setAlignParams({ ...alignParams, aligner: value })}
                    >
                      <SelectTrigger id="aligner">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="bwa">BWA-MEM (DNA, fast)</SelectItem>
                        <SelectItem value="bowtie2">Bowtie2 (sensitive)</SelectItem>
                        <SelectItem value="star">STAR (RNA-seq)</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="threads">CPU Threads</Label>
                    <Input
                      id="threads"
                      type="number"
                      value={alignParams.threads}
                      onChange={(e) => setAlignParams({ ...alignParams, threads: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Parallel processing cores</p>
                  </div>

                  <div className="space-y-2 flex items-center pt-6">
                    <Checkbox
                      id="mark-duplicates"
                      checked={alignParams.markDuplicates}
                      onCheckedChange={(checked) => setAlignParams({ ...alignParams, markDuplicates: checked as boolean })}
                    />
                    <label
                      htmlFor="mark-duplicates"
                      className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70 ml-2"
                    >
                      Mark Duplicates (Picard)
                    </label>
                  </div>
                </div>

                <Button
                  onClick={handleAlignment}
                  disabled={loading || !alignParams.sampleId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Aligning Reads...
                    </>
                  ) : (
                    'Run Alignment'
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Variant Calling Tab */}
        <TabsContent value="variants" className="space-y-4 mt-4">
          <Card>
            <CardHeader>
              <CardTitle>Variant Calling</CardTitle>
              <CardDescription>
                Identify SNPs, indels, and structural variants using GATK, FreeBayes, or BCFtools
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="variant-sample">Sample ID</Label>
                  <Input
                    id="variant-sample"
                    placeholder="Enter sample ID"
                    value={variantParams.sampleId}
                    onChange={(e) => setVariantParams({ ...variantParams, sampleId: e.target.value })}
                  />
                </div>

                <div className="space-y-2">
                  <Label htmlFor="caller">Variant Caller</Label>
                  <Select
                    value={variantParams.caller}
                    onValueChange={(value: any) => setVariantParams({ ...variantParams, caller: value })}
                  >
                    <SelectTrigger id="caller">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="gatk">GATK HaplotypeCaller (gold standard)</SelectItem>
                      <SelectItem value="freebayes">FreeBayes (haplotype-based)</SelectItem>
                      <SelectItem value="bcftools">BCFtools mpileup (fast)</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <Label htmlFor="min-depth">Min Read Depth</Label>
                    <Input
                      id="min-depth"
                      type="number"
                      value={variantParams.minDepth}
                      onChange={(e) => setVariantParams({ ...variantParams, minDepth: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Min coverage for variant call</p>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="min-var-quality">Min Quality Score</Label>
                    <Input
                      id="min-var-quality"
                      type="number"
                      value={variantParams.minQuality}
                      onChange={(e) => setVariantParams({ ...variantParams, minQuality: parseInt(e.target.value) })}
                    />
                    <p className="text-xs text-muted-foreground">Phred-scaled quality</p>
                  </div>
                </div>

                <div className="space-y-4">
                  <div className="flex items-center space-x-2">
                    <Checkbox
                      id="filter-snps"
                      checked={variantParams.filterSNPs}
                      onCheckedChange={(checked) => setVariantParams({ ...variantParams, filterSNPs: checked as boolean })}
                    />
                    <label
                      htmlFor="filter-snps"
                      className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70"
                    >
                      Apply SNP Filters (VQSR)
                    </label>
                  </div>
                  <p className="text-xs text-muted-foreground ml-6">
                    Variant Quality Score Recalibration for SNPs
                  </p>

                  <div className="flex items-center space-x-2">
                    <Checkbox
                      id="filter-indels"
                      checked={variantParams.filterIndels}
                      onCheckedChange={(checked) => setVariantParams({ ...variantParams, filterIndels: checked as boolean })}
                    />
                    <label
                      htmlFor="filter-indels"
                      className="text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70"
                    >
                      Apply Indel Filters (VQSR)
                    </label>
                  </div>
                  <p className="text-xs text-muted-foreground ml-6">
                    Variant Quality Score Recalibration for insertions/deletions
                  </p>
                </div>

                <Button
                  onClick={handleVariantCalling}
                  disabled={loading || !variantParams.sampleId}
                  className="w-full"
                >
                  {loading ? (
                    <>
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      Calling Variants...
                    </>
                  ) : (
                    'Run Variant Calling'
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

export default GenomicsAnalysis;
