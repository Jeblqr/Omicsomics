import { useState, useEffect } from 'react';
import { Node } from 'reactflow';
import { PipelineNodeData } from './PipelineNode';
import { getToolById, ToolParameter } from './ToolDefinitions';

interface ConfigPanelV2Props {
  selectedNode: Node<PipelineNodeData> | null;
  onUpdate: (nodeId: string, updates: Partial<PipelineNodeData>) => void;
  onClose: () => void;
}

const ConfigPanelV2 = ({ selectedNode, onUpdate, onClose }: ConfigPanelV2Props) => {
  const [label, setLabel] = useState('');
  const [tool, setTool] = useState('');
  const [version, setVersion] = useState('');
  const [parameters, setParameters] = useState<Record<string, any>>({});
  const [toolDefinition, setToolDefinition] = useState<any>(null);

  useEffect(() => {
    if (selectedNode) {
      setLabel(selectedNode.data.label || '');
      setTool(selectedNode.data.tool || '');
      setVersion(selectedNode.data.version || '');
      setParameters(selectedNode.data.parameters || {});

      // Try to find tool definition
      const toolDef = getToolById(selectedNode.id.replace(/^node_\d+_/, ''));
      if (toolDef) {
        setToolDefinition(toolDef);
        // Initialize missing parameters with defaults
        const newParams = { ...selectedNode.data.parameters };
        toolDef.parameterTemplate.forEach((param: ToolParameter) => {
          if (!(param.name in newParams) && param.default !== undefined) {
            newParams[param.name] = param.default;
          }
        });
        setParameters(newParams);
      } else {
        setToolDefinition(null);
      }
    }
  }, [selectedNode]);

  if (!selectedNode) {
    return (
      <div
        style={{
          width: '360px',
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

  const handleParamChange = (paramName: string, value: any) => {
    setParameters({
      ...parameters,
      [paramName]: value,
    });
  };

  const renderParameterInput = (param: ToolParameter) => {
    const value = parameters[param.name] ?? param.default ?? '';

    switch (param.type) {
      case 'boolean':
        return (
          <label
            style={{
              display: 'flex',
              alignItems: 'center',
              gap: '0.5rem',
              cursor: 'pointer',
            }}
          >
            <input
              type="checkbox"
              checked={value === true || value === 'true'}
              onChange={(e) => handleParamChange(param.name, e.target.checked)}
              style={{ width: '18px', height: '18px' }}
            />
            <span style={{ fontSize: '0.875rem', color: '#e5e7eb' }}>
              {param.label}
              {param.required && <span style={{ color: '#ef4444' }}> *</span>}
            </span>
          </label>
        );

      case 'select':
        return (
          <div>
            <label style={{ display: 'block', marginBottom: '0.25rem', fontSize: '0.875rem', color: '#e5e7eb' }}>
              {param.label}
              {param.required && <span style={{ color: '#ef4444' }}> *</span>}
            </label>
            <select
              value={value}
              onChange={(e) => handleParamChange(param.name, e.target.value)}
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
              {param.options?.map((option) => (
                <option key={option} value={option}>
                  {option}
                </option>
              ))}
            </select>
          </div>
        );

      case 'number':
        return (
          <div>
            <label style={{ display: 'block', marginBottom: '0.25rem', fontSize: '0.875rem', color: '#e5e7eb' }}>
              {param.label}
              {param.required && <span style={{ color: '#ef4444' }}> *</span>}
              {param.unit && <span style={{ color: '#9ca3af' }}> ({param.unit})</span>}
            </label>
            <input
              type="number"
              value={value}
              onChange={(e) => handleParamChange(param.name, parseFloat(e.target.value) || 0)}
              min={param.min}
              max={param.max}
              step={param.type === 'number' && param.min !== undefined && param.min < 1 ? 0.1 : 1}
              style={{
                width: '100%',
                padding: '0.5rem',
                backgroundColor: '#374151',
                border: '1px solid #4b5563',
                borderRadius: '4px',
                color: '#e5e7eb',
                fontSize: '0.875rem',
              }}
            />
            {(param.min !== undefined || param.max !== undefined) && (
              <div style={{ fontSize: '0.7rem', color: '#6b7280', marginTop: '0.25rem' }}>
                Range: {param.min ?? '-∞'} to {param.max ?? '∞'}
              </div>
            )}
          </div>
        );

      case 'file':
        return (
          <div>
            <label style={{ display: 'block', marginBottom: '0.25rem', fontSize: '0.875rem', color: '#e5e7eb' }}>
              {param.label}
              {param.required && <span style={{ color: '#ef4444' }}> *</span>}
            </label>
            <input
              type="text"
              value={value}
              onChange={(e) => handleParamChange(param.name, e.target.value)}
              placeholder="Enter file path..."
              style={{
                width: '100%',
                padding: '0.5rem',
                backgroundColor: '#374151',
                border: '1px solid #4b5563',
                borderRadius: '4px',
                color: '#e5e7eb',
                fontSize: '0.875rem',
              }}
            />
          </div>
        );

      case 'string':
      default:
        return (
          <div>
            <label style={{ display: 'block', marginBottom: '0.25rem', fontSize: '0.875rem', color: '#e5e7eb' }}>
              {param.label}
              {param.required && <span style={{ color: '#ef4444' }}> *</span>}
            </label>
            <input
              type="text"
              value={value}
              onChange={(e) => handleParamChange(param.name, e.target.value)}
              style={{
                width: '100%',
                padding: '0.5rem',
                backgroundColor: '#374151',
                border: '1px solid #4b5563',
                borderRadius: '4px',
                color: '#e5e7eb',
                fontSize: '0.875rem',
              }}
            />
          </div>
        );
    }
  };

  return (
    <div
      style={{
        width: '360px',
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
              padding: '0.25rem 0.5rem',
              backgroundColor: '#374151',
              color: '#9ca3af',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontSize: '0.875rem',
            }}
          >
            ✕
          </button>
        </div>
        <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>
          Node ID: {selectedNode.id}
        </div>
      </div>

      {/* Content */}
      <div style={{ flex: 1, overflowY: 'auto', padding: '1rem' }}>
        {/* Basic Info */}
        <div style={{ marginBottom: '1.5rem' }}>
          <h4 style={{ margin: '0 0 0.75rem 0', color: '#f3f4f6', fontSize: '0.9rem', fontWeight: 600 }}>
            Basic Information
          </h4>

          <div style={{ marginBottom: '1rem' }}>
            <label style={{ display: 'block', marginBottom: '0.25rem', fontSize: '0.875rem', color: '#e5e7eb' }}>
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
            />
          </div>

          <div style={{ marginBottom: '1rem' }}>
            <label style={{ display: 'block', marginBottom: '0.25rem', fontSize: '0.875rem', color: '#e5e7eb' }}>
              Tool
            </label>
            <input
              type="text"
              value={tool}
              disabled
              style={{
                width: '100%',
                padding: '0.5rem',
                backgroundColor: '#1f2937',
                border: '1px solid #4b5563',
                borderRadius: '4px',
                color: '#9ca3af',
                fontSize: '0.875rem',
              }}
            />
          </div>

          <div>
            <label style={{ display: 'block', marginBottom: '0.25rem', fontSize: '0.875rem', color: '#e5e7eb' }}>
              Version
            </label>
            <input
              type="text"
              value={version}
              disabled
              style={{
                width: '100%',
                padding: '0.5rem',
                backgroundColor: '#1f2937',
                border: '1px solid #4b5563',
                borderRadius: '4px',
                color: '#9ca3af',
                fontSize: '0.875rem',
              }}
            />
          </div>
        </div>

        {/* Parameters */}
        <div>
          <h4 style={{ margin: '0 0 0.75rem 0', color: '#f3f4f6', fontSize: '0.9rem', fontWeight: 600 }}>
            Parameters
          </h4>

          {toolDefinition && toolDefinition.parameterTemplate.length > 0 ? (
            <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
              {toolDefinition.parameterTemplate.map((param: ToolParameter) => (
                <div
                  key={param.name}
                  style={{
                    padding: '0.75rem',
                    backgroundColor: '#111827',
                    borderRadius: '6px',
                    border: '1px solid #374151',
                  }}
                >
                  {renderParameterInput(param)}
                  {param.description && (
                    <div style={{ marginTop: '0.5rem', fontSize: '0.7rem', color: '#6b7280', lineHeight: '1.4' }}>
                      ℹ️ {param.description}
                    </div>
                  )}
                </div>
              ))}
            </div>
          ) : (
            <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '6px', textAlign: 'center' }}>
              <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem' }}>
                No parameters defined for this tool
              </div>
              <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>
                This tool may not require configuration
              </div>
            </div>
          )}
        </div>

        {/* Documentation */}
        {toolDefinition?.documentation && (
          <div style={{ marginTop: '1.5rem' }}>
            <h4 style={{ margin: '0 0 0.5rem 0', color: '#f3f4f6', fontSize: '0.9rem', fontWeight: 600 }}>
              📖 Documentation
            </h4>
            <div
              style={{
                padding: '0.75rem',
                backgroundColor: '#111827',
                borderRadius: '6px',
                fontSize: '0.75rem',
                color: '#9ca3af',
                lineHeight: '1.5',
              }}
            >
              {toolDefinition.documentation}
            </div>
          </div>
        )}
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
          onClick={handleSave}
          style={{
            flex: 1,
            padding: '0.75rem',
            backgroundColor: '#059669',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontWeight: 600,
            fontSize: '0.875rem',
          }}
        >
          💾 Save Changes
        </button>
        <button
          onClick={onClose}
          style={{
            padding: '0.75rem 1rem',
            backgroundColor: '#374151',
            color: '#e5e7eb',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontWeight: 500,
            fontSize: '0.875rem',
          }}
        >
          Cancel
        </button>
      </div>
    </div>
  );
};

export default ConfigPanelV2;
