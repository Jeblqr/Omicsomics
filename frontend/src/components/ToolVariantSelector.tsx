import { useState, useEffect } from 'react';

interface ToolVariant {
  tool_id: string;
  tool_name: string;
  runtime: string;
  language: string;
  method: string;
  strengths: string;
  use_case: string;
  popularity_score: number;
  data_types: string[];
}

interface AnalysisFunction {
  function_id: string;
  display_name: string;
  category: string;
  description: string;
  data_types: string[];
  variants: ToolVariant[];
}

interface ToolVariantSelectorProps {
  functionId: string;
  onSelect: (toolId: string) => void;
  onClose: () => void;
}

const ToolVariantSelector = ({ functionId, onSelect, onClose }: ToolVariantSelectorProps) => {
  const [analysisFunction, setAnalysisFunction] = useState<AnalysisFunction | null>(null);
  const [selectedTool, setSelectedTool] = useState<string>('');
  const [compareMode, setCompareMode] = useState(false);
  const [compareTools, setCompareTools] = useState<Set<string>>(new Set());
  const [filterRuntime, setFilterRuntime] = useState<string>('all');

  useEffect(() => {
    // TODO: Load from API
    // Mock data for demonstration
    const mockFunction: AnalysisFunction = {
      function_id: 'quality_control',
      display_name: 'Quality Control',
      category: 'preprocessing',
      description: 'Assess data quality and filter low-quality samples/features',
      data_types: ['rna_seq', 'single_cell', 'proteomics'],
      variants: [
        {
          tool_id: 'seurat_qc',
          tool_name: 'Seurat QC',
          runtime: 'r',
          language: 'R',
          method: 'Feature/Cell filtering',
          strengths: 'Integrated workflow, rich visualization, widely used',
          use_case: 'Single-cell RNA-seq standard analysis',
          popularity_score: 95,
          data_types: ['single_cell'],
        },
        {
          tool_id: 'scanpy_qc',
          tool_name: 'scanpy QC',
          runtime: 'python',
          language: 'Python',
          method: 'Statistical filtering',
          strengths: 'Fast, memory-efficient, scalable to large datasets',
          use_case: 'Large-scale single-cell analysis (>100k cells)',
          popularity_score: 88,
          data_types: ['single_cell'],
        },
        {
          tool_id: 'fastqc',
          tool_name: 'FastQC',
          runtime: 'binary',
          language: 'Java',
          method: 'Raw sequencing QC',
          strengths: 'Industry standard, detailed HTML reports, handles all NGS data',
          use_case: 'Raw FASTQ quality assessment',
          popularity_score: 98,
          data_types: ['rna_seq', 'genomics'],
        },
        {
          tool_id: 'multiqc',
          tool_name: 'MultiQC',
          runtime: 'python',
          language: 'Python',
          method: 'Aggregate multiple QC reports',
          strengths: 'Combines multiple QC tools, beautiful reports',
          use_case: 'Batch sample QC summarization',
          popularity_score: 85,
          data_types: ['rna_seq', 'proteomics'],
        },
      ],
    };
    setAnalysisFunction(mockFunction);
  }, [functionId]);

  if (!analysisFunction) return null;

  const filteredVariants = analysisFunction.variants.filter((v) => {
    if (filterRuntime !== 'all' && v.runtime !== filterRuntime) return false;
    return true;
  });

  const sortedVariants = [...filteredVariants].sort((a, b) => b.popularity_score - a.popularity_score);

  const getRuntimeColor = (runtime: string) => {
    switch (runtime) {
      case 'r':
        return '#3b82f6';
      case 'python':
        return '#10b981';
      case 'binary':
        return '#f59e0b';
      default:
        return '#6b7280';
    }
  };

  const getRuntimeIcon = (runtime: string) => {
    switch (runtime) {
      case 'r':
        return '📊';
      case 'python':
        return '🐍';
      case 'binary':
        return '⚙️';
      default:
        return '🔧';
    }
  };

  const getPopularityStars = (score: number) => {
    const stars = Math.round(score / 20);
    return '★'.repeat(stars) + '☆'.repeat(5 - stars);
  };

  const handleCompareToggle = (toolId: string) => {
    const newCompare = new Set(compareTools);
    if (newCompare.has(toolId)) {
      newCompare.delete(toolId);
    } else {
      newCompare.add(toolId);
    }
    setCompareTools(newCompare);
  };

  return (
    <div
      style={{
        position: 'fixed',
        top: 0,
        left: 0,
        right: 0,
        bottom: 0,
        backgroundColor: 'rgba(0,0,0,0.85)',
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        zIndex: 4000,
        padding: '2rem',
      }}
      onClick={onClose}
    >
      <div
        style={{
          backgroundColor: '#1f2937',
          borderRadius: '16px',
          maxWidth: '1200px',
          width: '100%',
          maxHeight: '90vh',
          overflow: 'hidden',
          border: '2px solid #374151',
          display: 'flex',
          flexDirection: 'column',
        }}
        onClick={(e) => e.stopPropagation()}
      >
        {/* Header */}
        <div style={{ padding: '2rem', borderBottom: '2px solid #374151' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
            <div>
              <h2 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.75rem', fontWeight: 700 }}>
                Select Tool for {analysisFunction.display_name}
              </h2>
              <p style={{ color: '#9ca3af', fontSize: '1rem', marginTop: '0.5rem', marginBottom: 0 }}>
                {analysisFunction.description}
              </p>
            </div>
            <button
              onClick={onClose}
              style={{
                padding: '0.5rem 1rem',
                backgroundColor: '#374151',
                color: '#f3f4f6',
                border: 'none',
                borderRadius: '6px',
                cursor: 'pointer',
                fontWeight: 600,
              }}
            >
              ✕ Close
            </button>
          </div>

          {/* Filters */}
          <div style={{ marginTop: '1.5rem', display: 'flex', gap: '1rem', alignItems: 'center' }}>
            <div style={{ display: 'flex', gap: '0.5rem' }}>
              <button
                onClick={() => setFilterRuntime('all')}
                style={{
                  padding: '0.5rem 1rem',
                  backgroundColor: filterRuntime === 'all' ? '#3b82f6' : '#374151',
                  color: '#f3f4f6',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '0.875rem',
                  fontWeight: 600,
                }}
              >
                All Runtimes
              </button>
              <button
                onClick={() => setFilterRuntime('r')}
                style={{
                  padding: '0.5rem 1rem',
                  backgroundColor: filterRuntime === 'r' ? '#3b82f6' : '#374151',
                  color: '#f3f4f6',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '0.875rem',
                  fontWeight: 600,
                }}
              >
                📊 R
              </button>
              <button
                onClick={() => setFilterRuntime('python')}
                style={{
                  padding: '0.5rem 1rem',
                  backgroundColor: filterRuntime === 'python' ? '#3b82f6' : '#374151',
                  color: '#f3f4f6',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '0.875rem',
                  fontWeight: 600,
                }}
              >
                🐍 Python
              </button>
              <button
                onClick={() => setFilterRuntime('binary')}
                style={{
                  padding: '0.5rem 1rem',
                  backgroundColor: filterRuntime === 'binary' ? '#3b82f6' : '#374151',
                  color: '#f3f4f6',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '0.875rem',
                  fontWeight: 600,
                }}
              >
                ⚙️ Binary
              </button>
            </div>

            <div style={{ flex: 1 }} />

            {compareTools.size > 1 && (
              <button
                onClick={() => setCompareMode(true)}
                style={{
                  padding: '0.5rem 1rem',
                  backgroundColor: '#8b5cf6',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '0.875rem',
                  fontWeight: 600,
                }}
              >
                🔍 Compare ({compareTools.size})
              </button>
            )}
          </div>
        </div>

        {/* Tool List */}
        <div style={{ flex: 1, overflowY: 'auto', padding: '2rem' }}>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
            {sortedVariants.map((variant, index) => (
              <div
                key={variant.tool_id}
                style={{
                  backgroundColor: selectedTool === variant.tool_id ? '#1e3a8a' : '#111827',
                  border: `2px solid ${selectedTool === variant.tool_id ? '#3b82f6' : '#374151'}`,
                  borderRadius: '12px',
                  padding: '1.5rem',
                  cursor: 'pointer',
                  transition: 'all 0.2s',
                  position: 'relative',
                }}
                onClick={() => setSelectedTool(variant.tool_id)}
              >
                {/* Most Popular Badge */}
                {index === 0 && (
                  <div
                    style={{
                      position: 'absolute',
                      top: '-10px',
                      right: '1rem',
                      backgroundColor: '#f59e0b',
                      color: '#000',
                      padding: '0.25rem 0.75rem',
                      borderRadius: '12px',
                      fontSize: '0.75rem',
                      fontWeight: 700,
                    }}
                  >
                    🏆 Most Popular
                  </div>
                )}

                <div style={{ display: 'flex', gap: '1.5rem' }}>
                  {/* Radio Button */}
                  <div style={{ paddingTop: '0.25rem' }}>
                    <div
                      style={{
                        width: '24px',
                        height: '24px',
                        borderRadius: '50%',
                        border: `3px solid ${selectedTool === variant.tool_id ? '#3b82f6' : '#6b7280'}`,
                        backgroundColor: selectedTool === variant.tool_id ? '#3b82f6' : 'transparent',
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                      }}
                    >
                      {selectedTool === variant.tool_id && (
                        <div style={{ width: '10px', height: '10px', borderRadius: '50%', backgroundColor: 'white' }} />
                      )}
                    </div>
                  </div>

                  {/* Content */}
                  <div style={{ flex: 1 }}>
                    {/* Header */}
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '0.75rem' }}>
                      <div>
                        <h3 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.25rem', fontWeight: 700 }}>
                          {variant.tool_name}
                        </h3>
                        <div style={{ marginTop: '0.5rem', display: 'flex', gap: '0.75rem', alignItems: 'center' }}>
                          <span
                            style={{
                              padding: '0.25rem 0.75rem',
                              backgroundColor: getRuntimeColor(variant.runtime) + '33',
                              color: getRuntimeColor(variant.runtime),
                              borderRadius: '6px',
                              fontSize: '0.875rem',
                              fontWeight: 600,
                            }}
                          >
                            {getRuntimeIcon(variant.runtime)} {variant.language}
                          </span>
                          <span style={{ color: '#fbbf24', fontSize: '0.875rem' }}>
                            {getPopularityStars(variant.popularity_score)}
                          </span>
                          <span style={{ color: '#9ca3af', fontSize: '0.875rem' }}>
                            ({variant.popularity_score}/100)
                          </span>
                        </div>
                      </div>

                      <input
                        type="checkbox"
                        checked={compareTools.has(variant.tool_id)}
                        onChange={(e) => {
                          e.stopPropagation();
                          handleCompareToggle(variant.tool_id);
                        }}
                        style={{ width: '20px', height: '20px', cursor: 'pointer' }}
                        title="Add to comparison"
                      />
                    </div>

                    {/* Details */}
                    <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem', marginBottom: '1rem' }}>
                      <div>
                        <div style={{ color: '#9ca3af', fontSize: '0.75rem', marginBottom: '0.25rem' }}>Method</div>
                        <div style={{ color: '#d1d5db', fontSize: '0.875rem' }}>{variant.method}</div>
                      </div>
                      <div>
                        <div style={{ color: '#9ca3af', fontSize: '0.75rem', marginBottom: '0.25rem' }}>Use Case</div>
                        <div style={{ color: '#d1d5db', fontSize: '0.875rem' }}>{variant.use_case}</div>
                      </div>
                    </div>

                    {/* Strengths */}
                    <div>
                      <div style={{ color: '#9ca3af', fontSize: '0.75rem', marginBottom: '0.25rem' }}>✨ Strengths</div>
                      <div style={{ color: '#d1d5db', fontSize: '0.875rem', lineHeight: 1.5 }}>{variant.strengths}</div>
                    </div>

                    {/* Data Types */}
                    <div style={{ marginTop: '0.75rem', display: 'flex', gap: '0.5rem', flexWrap: 'wrap' }}>
                      {variant.data_types.map((type) => (
                        <span
                          key={type}
                          style={{
                            padding: '0.25rem 0.5rem',
                            backgroundColor: '#374151',
                            color: '#9ca3af',
                            borderRadius: '4px',
                            fontSize: '0.75rem',
                          }}
                        >
                          {type}
                        </span>
                      ))}
                    </div>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>

        {/* Footer */}
        <div style={{ padding: '1.5rem', borderTop: '2px solid #374151', backgroundColor: '#111827' }}>
          <div style={{ display: 'flex', gap: '1rem', justifyContent: 'flex-end' }}>
            <button
              onClick={onClose}
              style={{
                padding: '0.75rem 1.5rem',
                backgroundColor: '#374151',
                color: '#f3f4f6',
                border: 'none',
                borderRadius: '6px',
                cursor: 'pointer',
                fontSize: '1rem',
                fontWeight: 600,
              }}
            >
              Cancel
            </button>
            <button
              onClick={() => {
                if (selectedTool) {
                  onSelect(selectedTool);
                  onClose();
                }
              }}
              disabled={!selectedTool}
              style={{
                padding: '0.75rem 2rem',
                backgroundColor: selectedTool ? '#10b981' : '#374151',
                color: 'white',
                border: 'none',
                borderRadius: '6px',
                cursor: selectedTool ? 'pointer' : 'not-allowed',
                fontSize: '1rem',
                fontWeight: 600,
                opacity: selectedTool ? 1 : 0.5,
              }}
            >
              Add to Pipeline
            </button>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ToolVariantSelector;
