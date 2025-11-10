import { memo } from 'react';
import { Handle, Position, NodeProps } from 'reactflow';

export interface PipelineNodeData {
  label: string;
  tool?: string;
  version?: string;
  parameters?: Record<string, any>;
  nodeType: 'input' | 'process' | 'filter' | 'transform' | 'analysis' | 'output';
  onConfigure?: (nodeId: string) => void;
  onDelete?: (nodeId: string) => void;
}

const PipelineNode = memo(({ data, id, selected }: NodeProps<PipelineNodeData>) => {
  const getNodeColor = (type: string): string => {
    const colors: Record<string, string> = {
      input: '#10b981',
      process: '#3b82f6',
      filter: '#f59e0b',
      transform: '#8b5cf6',
      analysis: '#ef4444',
      output: '#6b7280',
    };
    return colors[type] || '#3b82f6';
  };

  const getNodeIcon = (type: string): string => {
    const icons: Record<string, string> = {
      input: '📥',
      process: '⚙️',
      filter: '🔍',
      transform: '🔄',
      analysis: '📊',
      output: '📤',
    };
    return icons[type] || '⚙️';
  };

  const backgroundColor = getNodeColor(data.nodeType);
  const borderColor = selected ? '#60a5fa' : backgroundColor;

  return (
    <div
      style={{
        background: '#1f2937',
        border: `2px solid ${borderColor}`,
        borderRadius: '8px',
        minWidth: '180px',
        boxShadow: selected
          ? '0 0 0 2px #3b82f6, 0 10px 25px -5px rgba(0, 0, 0, 0.5)'
          : '0 4px 6px -1px rgba(0, 0, 0, 0.3)',
        transition: 'all 0.2s',
      }}
    >
      {/* Input Handle */}
      {data.nodeType !== 'input' && (
        <Handle
          type="target"
          position={Position.Top}
          style={{
            background: backgroundColor,
            width: '12px',
            height: '12px',
            border: '2px solid #1f2937',
          }}
        />
      )}

      {/* Node Header */}
      <div
        style={{
          background: backgroundColor,
          padding: '0.5rem 0.75rem',
          borderTopLeftRadius: '6px',
          borderTopRightRadius: '6px',
          display: 'flex',
          alignItems: 'center',
          gap: '0.5rem',
        }}
      >
        <span style={{ fontSize: '1.2rem' }}>{getNodeIcon(data.nodeType)}</span>
        <span style={{ fontSize: '0.75rem', fontWeight: 600, color: 'white', textTransform: 'uppercase' }}>
          {data.nodeType}
        </span>
      </div>

      {/* Node Body */}
      <div style={{ padding: '0.75rem', color: '#e5e7eb' }}>
        <div style={{ fontWeight: 500, marginBottom: '0.25rem', fontSize: '0.9rem' }}>
          {data.label}
        </div>
        {data.tool && (
          <div style={{ fontSize: '0.75rem', color: '#9ca3af' }}>
            {data.tool} {data.version && `v${data.version}`}
          </div>
        )}
        {data.parameters && Object.keys(data.parameters).length > 0 && (
          <div style={{ fontSize: '0.7rem', color: '#6b7280', marginTop: '0.25rem' }}>
            {Object.keys(data.parameters).length} parameters
          </div>
        )}
      </div>

      {/* Node Actions */}
      <div
        style={{
          display: 'flex',
          gap: '0.25rem',
          padding: '0.5rem',
          borderTop: '1px solid #374151',
        }}
      >
        <button
          onClick={(e) => {
            e.stopPropagation();
            data.onConfigure?.(id);
          }}
          style={{
            flex: 1,
            padding: '0.25rem 0.5rem',
            fontSize: '0.7rem',
            backgroundColor: '#374151',
            color: '#e5e7eb',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontWeight: 500,
          }}
          title="Configure parameters"
        >
          ⚙️ Config
        </button>
        <button
          onClick={(e) => {
            e.stopPropagation();
            data.onDelete?.(id);
          }}
          style={{
            padding: '0.25rem 0.5rem',
            fontSize: '0.7rem',
            backgroundColor: '#7f1d1d',
            color: '#fecaca',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontWeight: 500,
          }}
          title="Delete node"
        >
          🗑️
        </button>
      </div>

      {/* Output Handle */}
      {data.nodeType !== 'output' && (
        <Handle
          type="source"
          position={Position.Bottom}
          style={{
            background: backgroundColor,
            width: '12px',
            height: '12px',
            border: '2px solid #1f2937',
          }}
        />
      )}
    </div>
  );
});

PipelineNode.displayName = 'PipelineNode';

export default PipelineNode;
