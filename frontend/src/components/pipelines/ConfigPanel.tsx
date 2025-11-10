import { useState, useEffect } from 'react';
import { Node } from 'reactflow';
import { PipelineNodeData } from './PipelineNode';

interface ConfigPanelProps {
  selectedNode: Node<PipelineNodeData> | null;
  onUpdate: (nodeId: string, updates: Partial<PipelineNodeData>) => void;
  onClose: () => void;
}

const ConfigPanel = ({ selectedNode, onUpdate, onClose }: ConfigPanelProps) => {
  const [label, setLabel] = useState('');
  const [tool, setTool] = useState('');
  const [version, setVersion] = useState('');
  const [parameters, setParameters] = useState<Record<string, any>>({});
  const [newParamKey, setNewParamKey] = useState('');
  const [newParamValue, setNewParamValue] = useState('');

  useEffect(() => {
    if (selectedNode) {
      setLabel(selectedNode.data.label || '');
      setTool(selectedNode.data.tool || '');
      setVersion(selectedNode.data.version || '');
      setParameters(selectedNode.data.parameters || {});
    }
  }, [selectedNode]);

  if (!selectedNode) {
    return (
      <div
        style={{
          width: '320px',
          height: '100%',
          backgroundColor: '#1f2937',
          borderLeft: '1px solid #374151',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          color: '#9ca3af',
          fontSize: '0.875rem',
          textAlign: 'center',
          padding: '2rem',
        }}
      >
        Select a node to configure its parameters
      </div>
    );
  }

  const handleSave = () => {
    onUpdate(selectedNode.id, {
      label,
      tool,
      version,
      parameters,
    });
    onClose();
  };

  const handleAddParameter = () => {
    if (newParamKey.trim()) {
      setParameters({
        ...parameters,
        [newParamKey]: newParamValue || '',
      });
      setNewParamKey('');
      setNewParamValue('');
    }
  };

  const handleRemoveParameter = (key: string) => {
    const newParams = { ...parameters };
    delete newParams[key];
    setParameters(newParams);
  };

  const handleParamValueChange = (key: string, value: string) => {
    setParameters({
      ...parameters,
      [key]: value,
    });
  };

  return (
    <div
      style={{
        width: '320px',
        height: '100%',
        backgroundColor: '#1f2937',
        borderLeft: '1px solid #374151',
        display: 'flex',
        flexDirection: 'column',
        overflow: 'hidden',
      }}
    >
      {/* Header */}
      <div style={{ padding: '1rem', borderBottom: '1px solid #374151' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '0.5rem' }}>
          <h3 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.1rem' }}>
            ⚙️ Configure Node
          </h3>
          <button
            onClick={onClose}
            style={{
              background: 'none',
              border: 'none',
              color: '#9ca3af',
              cursor: 'pointer',
              fontSize: '1.5rem',
              padding: 0,
              lineHeight: 1,
            }}
            title="Close panel"
          >
            ×
          </button>
        </div>
        <div
          style={{
            fontSize: '0.75rem',
            color: '#6b7280',
            textTransform: 'uppercase',
            fontWeight: 600,
          }}
        >
          {selectedNode.data.nodeType} NODE
        </div>
      </div>

      {/* Form */}
      <div style={{ flex: 1, overflowY: 'auto', padding: '1rem' }}>
        {/* Label */}
        <div style={{ marginBottom: '1rem' }}>
          <label
            style={{
              display: 'block',
              marginBottom: '0.5rem',
              fontSize: '0.875rem',
              fontWeight: 500,
              color: '#f3f4f6',
            }}
          >
            Step Name
          </label>
          <input
            type="text"
            value={label}
            onChange={(e) => setLabel(e.target.value)}
            style={{
              width: '100%',
              padding: '0.5rem',
              backgroundColor: '#374151',
              border: '1px solid #4b5563',
              borderRadius: '4px',
              color: '#e5e7eb',
              fontSize: '0.875rem',
            }}
            placeholder="e.g., Quality Control"
          />
        </div>

        {/* Tool */}
        <div style={{ marginBottom: '1rem' }}>
          <label
            style={{
              display: 'block',
              marginBottom: '0.5rem',
              fontSize: '0.875rem',
              fontWeight: 500,
              color: '#f3f4f6',
            }}
          >
            Tool
          </label>
          <input
            type="text"
            value={tool}
            onChange={(e) => setTool(e.target.value)}
            style={{
              width: '100%',
              padding: '0.5rem',
              backgroundColor: '#374151',
              border: '1px solid #4b5563',
              borderRadius: '4px',
              color: '#e5e7eb',
              fontSize: '0.875rem',
            }}
            placeholder="e.g., fastqc"
          />
        </div>

        {/* Version */}
        <div style={{ marginBottom: '1rem' }}>
          <label
            style={{
              display: 'block',
              marginBottom: '0.5rem',
              fontSize: '0.875rem',
              fontWeight: 500,
              color: '#f3f4f6',
            }}
          >
            Version
          </label>
          <input
            type="text"
            value={version}
            onChange={(e) => setVersion(e.target.value)}
            style={{
              width: '100%',
              padding: '0.5rem',
              backgroundColor: '#374151',
              border: '1px solid #4b5563',
              borderRadius: '4px',
              color: '#e5e7eb',
              fontSize: '0.875rem',
            }}
            placeholder="e.g., 0.11.9"
          />
        </div>

        {/* Parameters */}
        <div style={{ marginBottom: '1rem' }}>
          <div style={{ marginBottom: '0.5rem', fontSize: '0.875rem', fontWeight: 500, color: '#f3f4f6' }}>
            Parameters
          </div>

          {/* Existing Parameters */}
          {Object.entries(parameters).map(([key, value]) => (
            <div
              key={key}
              style={{
                display: 'flex',
                gap: '0.5rem',
                marginBottom: '0.5rem',
                alignItems: 'center',
              }}
            >
              <input
                type="text"
                value={key}
                disabled
                style={{
                  flex: 1,
                  padding: '0.4rem',
                  backgroundColor: '#2d3748',
                  border: '1px solid #4b5563',
                  borderRadius: '4px',
                  color: '#9ca3af',
                  fontSize: '0.8rem',
                }}
              />
              <input
                type="text"
                value={String(value)}
                onChange={(e) => handleParamValueChange(key, e.target.value)}
                style={{
                  flex: 1,
                  padding: '0.4rem',
                  backgroundColor: '#374151',
                  border: '1px solid #4b5563',
                  borderRadius: '4px',
                  color: '#e5e7eb',
                  fontSize: '0.8rem',
                }}
              />
              <button
                onClick={() => handleRemoveParameter(key)}
                style={{
                  padding: '0.4rem 0.6rem',
                  backgroundColor: '#7f1d1d',
                  color: '#fecaca',
                  border: 'none',
                  borderRadius: '4px',
                  cursor: 'pointer',
                  fontSize: '0.75rem',
                }}
                title="Remove parameter"
              >
                🗑️
              </button>
            </div>
          ))}

          {/* Add New Parameter */}
          <div style={{ marginTop: '0.75rem', paddingTop: '0.75rem', borderTop: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.5rem' }}>
              Add New Parameter
            </div>
            <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
              <input
                type="text"
                value={newParamKey}
                onChange={(e) => setNewParamKey(e.target.value)}
                style={{
                  padding: '0.4rem',
                  backgroundColor: '#374151',
                  border: '1px solid #4b5563',
                  borderRadius: '4px',
                  color: '#e5e7eb',
                  fontSize: '0.8rem',
                }}
                placeholder="Parameter name"
              />
              <input
                type="text"
                value={newParamValue}
                onChange={(e) => setNewParamValue(e.target.value)}
                onKeyPress={(e) => {
                  if (e.key === 'Enter') {
                    handleAddParameter();
                  }
                }}
                style={{
                  padding: '0.4rem',
                  backgroundColor: '#374151',
                  border: '1px solid #4b5563',
                  borderRadius: '4px',
                  color: '#e5e7eb',
                  fontSize: '0.8rem',
                }}
                placeholder="Parameter value"
              />
              <button
                onClick={handleAddParameter}
                disabled={!newParamKey.trim()}
                style={{
                  padding: '0.5rem',
                  backgroundColor: newParamKey.trim() ? '#059669' : '#374151',
                  color: newParamKey.trim() ? 'white' : '#6b7280',
                  border: 'none',
                  borderRadius: '4px',
                  cursor: newParamKey.trim() ? 'pointer' : 'not-allowed',
                  fontSize: '0.8rem',
                  fontWeight: 500,
                }}
              >
                + Add Parameter
              </button>
            </div>
          </div>
        </div>
      </div>

      {/* Footer */}
      <div
        style={{
          padding: '1rem',
          borderTop: '1px solid #374151',
          backgroundColor: '#111827',
          display: 'flex',
          gap: '0.5rem',
        }}
      >
        <button
          onClick={onClose}
          style={{
            flex: 1,
            padding: '0.75rem',
            backgroundColor: '#374151',
            color: '#e5e7eb',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            fontWeight: 500,
          }}
        >
          Cancel
        </button>
        <button
          onClick={handleSave}
          style={{
            flex: 1,
            padding: '0.75rem',
            backgroundColor: '#3b82f6',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            fontWeight: 500,
          }}
        >
          Save
        </button>
      </div>
    </div>
  );
};

export default ConfigPanel;
