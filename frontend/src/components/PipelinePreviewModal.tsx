import { useState } from 'react';
import { Node, Edge } from 'reactflow';
import { PipelineNodeData } from './pipelines/PipelineNode';

interface PreviewResult {
  nodeId: string;
  nodeLabel: string;
  inputSample: any[];
  outputSample: any[];
  impact: {
    inputRows: number;
    outputRows: number;
    rowsFiltered: number;
    filterRate: number;
    columnsAdded?: number;
    columnsRemoved?: number;
    transformations?: string[];
  };
}

interface PipelinePreviewModalProps {
  isOpen: boolean;
  onClose: () => void;
  nodes: Node<PipelineNodeData>[];
  edges: Edge[];
  onExecute?: () => void;
}

const PipelinePreviewModal = ({ isOpen, onClose, nodes, edges, onExecute }: PipelinePreviewModalProps) => {
  const [previewResults, setPreviewResults] = useState<PreviewResult[]>([]);
  const [isGenerating, setIsGenerating] = useState(false);

  if (!isOpen) return null;

  // 生成模拟预览结果
  const generatePreview = () => {
    setIsGenerating(true);
    
    // 模拟异步处理
    setTimeout(() => {
      const results: PreviewResult[] = [];
      
      // 为每个节点生成预览
      nodes.forEach((node) => {
        const inputSample = Array.from({ length: 5 }, (_, i) => ({
          id: i + 1,
          sample_name: `Sample_${i + 1}`,
          value: Math.random() * 100,
          quality_score: Math.random() * 10,
        }));

        // 根据节点类型模拟不同的输出
        let outputRows = inputSample.length;
        let rowsFiltered = 0;
        let transformations: string[] = [];

        if (node.data.nodeType === 'filter' || node.data.tool?.includes('QC')) {
          rowsFiltered = Math.floor(Math.random() * 2);
          outputRows = inputSample.length - rowsFiltered;
          transformations.push('Quality filtering applied');
        } else if (node.data.nodeType === 'transform') {
          transformations.push('Data transformation applied');
        } else if (node.data.nodeType === 'analysis') {
          transformations.push('Statistical analysis performed');
        }

        const outputSample = inputSample.slice(0, outputRows).map((row) => ({
          ...row,
          processed: true,
          timestamp: new Date().toISOString(),
        }));

        results.push({
          nodeId: node.id,
          nodeLabel: node.data.label,
          inputSample,
          outputSample,
          impact: {
            inputRows: inputSample.length,
            outputRows,
            rowsFiltered,
            filterRate: rowsFiltered / inputSample.length,
            transformations,
          },
        });
      });

      setPreviewResults(results);
      setIsGenerating(false);
    }, 1000);
  };

  const handleExecute = () => {
    if (onExecute) {
      onExecute();
    }
    onClose();
  };

  return (
    <div
      style={{
        position: 'fixed',
        top: 0,
        left: 0,
        right: 0,
        bottom: 0,
        backgroundColor: 'rgba(0,0,0,0.8)',
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        zIndex: 3000,
      }}
      onClick={onClose}
    >
      <div
        style={{
          backgroundColor: '#1f2937',
          borderRadius: '12px',
          padding: '2rem',
          maxWidth: '1200px',
          width: '90%',
          maxHeight: '90vh',
          overflow: 'auto',
          border: '2px solid #374151',
        }}
        onClick={(e) => e.stopPropagation()}
      >
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1.5rem' }}>
          <h2 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.5rem', fontWeight: 700 }}>
            🔍 Pipeline Preview (Dry Run)
          </h2>
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

        <div style={{ marginBottom: '1.5rem', padding: '1rem', backgroundColor: '#111827', borderRadius: '8px', border: '1px solid #374151' }}>
          <p style={{ margin: 0, color: '#9ca3af', fontSize: '0.875rem', lineHeight: 1.6 }}>
            This is a <strong style={{ color: '#f3f4f6' }}>dry run preview</strong> that shows what your pipeline will do without actually executing it.
            Review the sample data and impact summary for each node before running the full pipeline.
          </p>
        </div>

        {previewResults.length === 0 ? (
          <div style={{ textAlign: 'center', padding: '3rem' }}>
            <button
              onClick={generatePreview}
              disabled={isGenerating}
              style={{
                padding: '1rem 2rem',
                backgroundColor: isGenerating ? '#6b7280' : '#3b82f6',
                color: 'white',
                border: 'none',
                borderRadius: '8px',
                cursor: isGenerating ? 'not-allowed' : 'pointer',
                fontSize: '1rem',
                fontWeight: 600,
              }}
            >
              {isGenerating ? '⏳ Generating Preview...' : '🚀 Generate Preview'}
            </button>
          </div>
        ) : (
          <div>
            {/* Overall Summary */}
            <div style={{ marginBottom: '2rem', padding: '1.5rem', backgroundColor: '#111827', borderRadius: '8px', border: '1px solid #374151' }}>
              <h3 style={{ margin: '0 0 1rem 0', color: '#f3f4f6', fontSize: '1.125rem', fontWeight: 600 }}>
                📊 Overall Impact
              </h3>
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '1rem' }}>
                <div style={{ padding: '1rem', backgroundColor: '#1f2937', borderRadius: '6px', border: '1px solid #374151' }}>
                  <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem', textTransform: 'uppercase', fontWeight: 600 }}>
                    Total Nodes
                  </div>
                  <div style={{ fontSize: '1.5rem', color: '#3b82f6', fontWeight: 700 }}>
                    {previewResults.length}
                  </div>
                </div>
                <div style={{ padding: '1rem', backgroundColor: '#1f2937', borderRadius: '6px', border: '1px solid #374151' }}>
                  <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem', textTransform: 'uppercase', fontWeight: 600 }}>
                    Total Rows Filtered
                  </div>
                  <div style={{ fontSize: '1.5rem', color: '#ef4444', fontWeight: 700 }}>
                    {previewResults.reduce((sum, r) => sum + r.impact.rowsFiltered, 0)}
                  </div>
                </div>
                <div style={{ padding: '1rem', backgroundColor: '#1f2937', borderRadius: '6px', border: '1px solid #374151' }}>
                  <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem', textTransform: 'uppercase', fontWeight: 600 }}>
                    Avg Filter Rate
                  </div>
                  <div style={{ fontSize: '1.5rem', color: '#f59e0b', fontWeight: 700 }}>
                    {(previewResults.reduce((sum, r) => sum + r.impact.filterRate, 0) / previewResults.length * 100).toFixed(1)}%
                  </div>
                </div>
              </div>
            </div>

            {/* Node Results */}
            <div style={{ display: 'flex', flexDirection: 'column', gap: '1.5rem' }}>
              {previewResults.map((result) => (
                <div
                  key={result.nodeId}
                  style={{
                    padding: '1.5rem',
                    backgroundColor: '#111827',
                    borderRadius: '8px',
                    border: '1px solid #374151',
                  }}
                >
                  <h4 style={{ margin: '0 0 1rem 0', color: '#f3f4f6', fontSize: '1rem', fontWeight: 600 }}>
                    {result.nodeLabel}
                  </h4>

                  {/* Impact Summary */}
                  <div style={{ marginBottom: '1rem', display: 'flex', gap: '1rem', flexWrap: 'wrap' }}>
                    <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>
                      Input: <span style={{ color: '#3b82f6', fontWeight: 600 }}>{result.impact.inputRows} rows</span>
                    </div>
                    <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>
                      Output: <span style={{ color: '#10b981', fontWeight: 600 }}>{result.impact.outputRows} rows</span>
                    </div>
                    {result.impact.rowsFiltered > 0 && (
                      <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>
                        Filtered: <span style={{ color: '#ef4444', fontWeight: 600 }}>{result.impact.rowsFiltered} rows ({(result.impact.filterRate * 100).toFixed(1)}%)</span>
                      </div>
                    )}
                  </div>

                  {/* Transformations */}
                  {result.impact.transformations && result.impact.transformations.length > 0 && (
                    <div style={{ marginBottom: '1rem' }}>
                      <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.5rem', textTransform: 'uppercase', fontWeight: 600 }}>
                        Transformations
                      </div>
                      <div style={{ display: 'flex', gap: '0.5rem', flexWrap: 'wrap' }}>
                        {result.impact.transformations.map((t, idx) => (
                          <span
                            key={idx}
                            style={{
                              padding: '0.25rem 0.75rem',
                              backgroundColor: '#1f2937',
                              color: '#8b5cf6',
                              fontSize: '0.75rem',
                              borderRadius: '4px',
                              border: '1px solid #374151',
                            }}
                          >
                            {t}
                          </span>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* Sample Data Preview */}
                  <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem' }}>
                    {/* Input Sample */}
                    <div>
                      <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.5rem', textTransform: 'uppercase', fontWeight: 600 }}>
                        Input Sample (First 5 rows)
                      </div>
                      <div style={{ backgroundColor: '#1f2937', borderRadius: '6px', overflow: 'hidden', border: '1px solid #374151' }}>
                        <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.75rem' }}>
                          <thead>
                            <tr style={{ backgroundColor: '#111827' }}>
                              {Object.keys(result.inputSample[0] || {}).map((key) => (
                                <th key={key} style={{ padding: '0.5rem', color: '#f3f4f6', fontWeight: 600, textAlign: 'left', borderBottom: '1px solid #374151' }}>
                                  {key}
                                </th>
                              ))}
                            </tr>
                          </thead>
                          <tbody>
                            {result.inputSample.map((row, idx) => (
                              <tr key={idx} style={{ backgroundColor: idx % 2 === 0 ? '#1f2937' : '#111827' }}>
                                {Object.values(row).map((val: any, i) => (
                                  <td key={i} style={{ padding: '0.5rem', color: '#e5e7eb', borderBottom: '1px solid #374151' }}>
                                    {typeof val === 'number' ? val.toFixed(2) : String(val)}
                                  </td>
                                ))}
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                    </div>

                    {/* Output Sample */}
                    <div>
                      <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.5rem', textTransform: 'uppercase', fontWeight: 600 }}>
                        Output Sample (First 5 rows)
                      </div>
                      <div style={{ backgroundColor: '#1f2937', borderRadius: '6px', overflow: 'hidden', border: '1px solid #374151' }}>
                        <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.75rem' }}>
                          <thead>
                            <tr style={{ backgroundColor: '#111827' }}>
                              {Object.keys(result.outputSample[0] || {}).map((key) => (
                                <th key={key} style={{ padding: '0.5rem', color: '#f3f4f6', fontWeight: 600, textAlign: 'left', borderBottom: '1px solid #374151' }}>
                                  {key}
                                </th>
                              ))}
                            </tr>
                          </thead>
                          <tbody>
                            {result.outputSample.map((row, idx) => (
                              <tr key={idx} style={{ backgroundColor: idx % 2 === 0 ? '#1f2937' : '#111827' }}>
                                {Object.values(row).map((val: any, i) => (
                                  <td key={i} style={{ padding: '0.5rem', color: '#e5e7eb', borderBottom: '1px solid #374151' }}>
                                    {typeof val === 'number' ? val.toFixed(2) : String(val).length > 20 ? String(val).substring(0, 20) + '...' : String(val)}
                                  </td>
                                ))}
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                    </div>
                  </div>
                </div>
              ))}
            </div>

            {/* Action Buttons */}
            <div style={{ marginTop: '2rem', display: 'flex', gap: '1rem', justifyContent: 'flex-end' }}>
              <button
                onClick={onClose}
                style={{
                  padding: '0.75rem 1.5rem',
                  backgroundColor: '#374151',
                  color: '#f3f4f6',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontWeight: 600,
                }}
              >
                Cancel
              </button>
              <button
                onClick={handleExecute}
                style={{
                  padding: '0.75rem 1.5rem',
                  backgroundColor: '#10b981',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontWeight: 600,
                }}
              >
                ✓ Execute Pipeline
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default PipelinePreviewModal;
