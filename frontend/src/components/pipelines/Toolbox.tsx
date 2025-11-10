import { useState } from 'react';

export interface ToolDefinition {
  id: string;
  name: string;
  type: 'input' | 'process' | 'filter' | 'transform' | 'analysis' | 'output';
  tool: string;
  version: string;
  category: string;
  description: string;
  parameters?: Record<string, any>;
}

interface ToolboxProps {
  onDragStart: (tool: ToolDefinition) => void;
}

const availableTools: ToolDefinition[] = [
  // Input tools
  {
    id: 'fastq-input',
    name: 'FASTQ Input',
    type: 'input',
    tool: 'file-input',
    version: '1.0',
    category: 'Input',
    description: 'Read FASTQ sequencing files',
    parameters: { file_type: 'fastq' },
  },
  {
    id: 'vcf-input',
    name: 'VCF Input',
    type: 'input',
    tool: 'file-input',
    version: '1.0',
    category: 'Input',
    description: 'Read VCF variant files',
    parameters: { file_type: 'vcf' },
  },
  {
    id: 'csv-input',
    name: 'CSV Input',
    type: 'input',
    tool: 'file-input',
    version: '1.0',
    category: 'Input',
    description: 'Read CSV data files',
    parameters: { file_type: 'csv' },
  },

  // Quality Control
  {
    id: 'fastqc',
    name: 'Quality Control',
    type: 'process',
    tool: 'fastqc',
    version: '0.11.9',
    category: 'Quality Control',
    description: 'Assess sequence quality',
    parameters: { threads: 4 },
  },
  {
    id: 'trim-galore',
    name: 'Adapter Trimming',
    type: 'filter',
    tool: 'trim_galore',
    version: '0.6.7',
    category: 'Quality Control',
    description: 'Remove adapters and low-quality bases',
    parameters: { quality: 20, length: 20 },
  },

  // Alignment
  {
    id: 'bwa-mem',
    name: 'BWA Alignment',
    type: 'process',
    tool: 'bwa',
    version: '0.7.17',
    category: 'Alignment',
    description: 'Align reads to reference genome',
    parameters: { threads: 8, algorithm: 'mem' },
  },
  {
    id: 'star',
    name: 'STAR Alignment',
    type: 'process',
    tool: 'star',
    version: '2.7.10',
    category: 'Alignment',
    description: 'RNA-seq read alignment',
    parameters: { threads: 8, outSAMtype: 'BAM SortedByCoordinate' },
  },

  // Analysis
  {
    id: 'variant-calling',
    name: 'Variant Calling',
    type: 'analysis',
    tool: 'haplotypecaller',
    version: '4.2.0',
    category: 'Variant Analysis',
    description: 'Call genomic variants',
    parameters: { ploidy: 2 },
  },
  {
    id: 'deseq2',
    name: 'Differential Expression',
    type: 'analysis',
    tool: 'deseq2',
    version: '1.34.0',
    category: 'RNA-seq',
    description: 'Identify differentially expressed genes',
    parameters: { alpha: 0.05, lfcThreshold: 0 },
  },
  {
    id: 'feature-counts',
    name: 'Feature Counting',
    type: 'process',
    tool: 'featureCounts',
    version: '2.0.1',
    category: 'RNA-seq',
    description: 'Count reads per gene',
    parameters: { threads: 4, isPairedEnd: true },
  },

  // Transformation
  {
    id: 'samtools-sort',
    name: 'Sort BAM',
    type: 'transform',
    tool: 'samtools',
    version: '1.15',
    category: 'Data Processing',
    description: 'Sort BAM alignment files',
    parameters: { threads: 4 },
  },
  {
    id: 'normalize',
    name: 'Normalize Data',
    type: 'transform',
    tool: 'custom',
    version: '1.0',
    category: 'Data Processing',
    description: 'Normalize expression data',
    parameters: { method: 'TMM' },
  },

  // Output
  {
    id: 'vcf-output',
    name: 'VCF Output',
    type: 'output',
    tool: 'file-output',
    version: '1.0',
    category: 'Output',
    description: 'Save variants in VCF format',
    parameters: { format: 'vcf' },
  },
  {
    id: 'csv-output',
    name: 'CSV Output',
    type: 'output',
    tool: 'file-output',
    version: '1.0',
    category: 'Output',
    description: 'Save results in CSV format',
    parameters: { format: 'csv' },
  },
  {
    id: 'report-output',
    name: 'HTML Report',
    type: 'output',
    tool: 'report-generator',
    version: '1.0',
    category: 'Output',
    description: 'Generate analysis report',
    parameters: { format: 'html' },
  },
];

const Toolbox = ({ onDragStart }: ToolboxProps) => {
  const [selectedCategory, setSelectedCategory] = useState<string>('All');
  const [searchQuery, setSearchQuery] = useState('');

  const categories = ['All', ...Array.from(new Set(availableTools.map((t) => t.category)))];

  const filteredTools = availableTools.filter((tool) => {
    const matchesCategory = selectedCategory === 'All' || tool.category === selectedCategory;
    const matchesSearch =
      searchQuery === '' ||
      tool.name.toLowerCase().includes(searchQuery.toLowerCase()) ||
      tool.description.toLowerCase().includes(searchQuery.toLowerCase());
    return matchesCategory && matchesSearch;
  });

  const handleDragStart = (event: React.DragEvent, tool: ToolDefinition) => {
    event.dataTransfer.effectAllowed = 'move';
    event.dataTransfer.setData('application/reactflow', JSON.stringify(tool));
    onDragStart(tool);
  };

  return (
    <div
      style={{
        width: '280px',
        height: '100%',
        backgroundColor: '#1f2937',
        borderRight: '1px solid #374151',
        display: 'flex',
        flexDirection: 'column',
        overflow: 'hidden',
      }}
    >
      {/* Header */}
      <div style={{ padding: '1rem', borderBottom: '1px solid #374151' }}>
        <h3 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.1rem', marginBottom: '0.75rem' }}>
          🧰 Tool Library
        </h3>

        {/* Search */}
        <input
          type="text"
          placeholder="Search tools..."
          value={searchQuery}
          onChange={(e) => setSearchQuery(e.target.value)}
          style={{
            width: '100%',
            padding: '0.5rem',
            backgroundColor: '#374151',
            border: '1px solid #4b5563',
            borderRadius: '4px',
            color: '#e5e7eb',
            fontSize: '0.875rem',
            marginBottom: '0.75rem',
          }}
        />

        {/* Category Filter */}
        <select
          value={selectedCategory}
          onChange={(e) => setSelectedCategory(e.target.value)}
          style={{
            width: '100%',
            padding: '0.5rem',
            backgroundColor: '#374151',
            border: '1px solid #4b5563',
            borderRadius: '4px',
            color: '#e5e7eb',
            fontSize: '0.875rem',
          }}
        >
          {categories.map((cat) => (
            <option key={cat} value={cat}>
              {cat}
            </option>
          ))}
        </select>
      </div>

      {/* Tools List */}
      <div style={{ flex: 1, overflowY: 'auto', padding: '0.5rem' }}>
        {filteredTools.length === 0 ? (
          <div style={{ padding: '2rem 1rem', textAlign: 'center', color: '#9ca3af', fontSize: '0.875rem' }}>
            No tools found
          </div>
        ) : (
          filteredTools.map((tool) => (
            <div
              key={tool.id}
              draggable
              onDragStart={(e) => handleDragStart(e, tool)}
              style={{
                padding: '0.75rem',
                marginBottom: '0.5rem',
                backgroundColor: '#374151',
                borderRadius: '6px',
                cursor: 'grab',
                border: '1px solid #4b5563',
                transition: 'all 0.2s',
              }}
              onMouseEnter={(e) => {
                e.currentTarget.style.backgroundColor = '#4b5563';
                e.currentTarget.style.transform = 'translateX(4px)';
              }}
              onMouseLeave={(e) => {
                e.currentTarget.style.backgroundColor = '#374151';
                e.currentTarget.style.transform = 'translateX(0)';
              }}
            >
              <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginBottom: '0.25rem' }}>
                <span style={{ fontSize: '1rem' }}>
                  {tool.type === 'input' && '📥'}
                  {tool.type === 'process' && '⚙️'}
                  {tool.type === 'filter' && '🔍'}
                  {tool.type === 'transform' && '🔄'}
                  {tool.type === 'analysis' && '📊'}
                  {tool.type === 'output' && '📤'}
                </span>
                <div style={{ flex: 1 }}>
                  <div style={{ fontWeight: 500, color: '#f3f4f6', fontSize: '0.875rem' }}>
                    {tool.name}
                  </div>
                  <div style={{ fontSize: '0.7rem', color: '#9ca3af' }}>
                    {tool.tool} {tool.version}
                  </div>
                </div>
              </div>
              <div style={{ fontSize: '0.75rem', color: '#9ca3af', lineHeight: '1.3' }}>
                {tool.description}
              </div>
            </div>
          ))
        )}
      </div>

      {/* Footer */}
      <div
        style={{
          padding: '0.75rem',
          borderTop: '1px solid #374151',
          backgroundColor: '#111827',
          fontSize: '0.75rem',
          color: '#6b7280',
          textAlign: 'center',
        }}
      >
        Drag tools to canvas to build pipeline
      </div>
    </div>
  );
};

export default Toolbox;
export { availableTools };
