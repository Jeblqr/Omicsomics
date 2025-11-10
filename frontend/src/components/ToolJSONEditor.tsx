/**
 * Tool JSON Editor
 * 在线编辑工具JSON定义
 * 支持语法高亮、验证、模板
 */

import React, { useState } from 'react';
import { ToolDefinitionSchema, validateToolSchema, TOOL_SCHEMA_EXAMPLE, TOOL_SCHEMA_VERSION } from '../schemas/ToolSchema';
import { toolRegistry } from '../utils/toolRegistry';

interface ToolJSONEditorProps {
  toolId?: string;  // 如果提供，则编辑现有工具；否则创建新工具
  onSave?: (tool: ToolDefinitionSchema) => void;
  onCancel?: () => void;
}

const ToolJSONEditor: React.FC<ToolJSONEditorProps> = ({ toolId, onSave, onCancel }) => {
  const [jsonText, setJsonText] = useState(() => {
    if (toolId) {
      const tool = toolRegistry.getTool(toolId);
      return JSON.stringify(tool, null, 2);
    }
    return JSON.stringify(TOOL_SCHEMA_EXAMPLE, null, 2);
  });

  const [validationErrors, setValidationErrors] = useState<string[]>([]);
  const [isValid, setIsValid] = useState(true);

  const handleValidate = () => {
    try {
      const tool = JSON.parse(jsonText);
      const validation = validateToolSchema(tool);
      setValidationErrors(validation.errors);
      setIsValid(validation.valid);
      return validation.valid;
    } catch (error) {
      setValidationErrors([`JSON Parse Error: ${error}`]);
      setIsValid(false);
      return false;
    }
  };

  const handleSave = () => {
    if (!handleValidate()) {
      alert('Please fix validation errors before saving');
      return;
    }

    try {
      const tool = JSON.parse(jsonText);
      
      if (toolId) {
        // 更新现有工具
        const result = toolRegistry.updateCustomTool(toolId, tool);
        if (!result.success) {
          alert(`Failed to update tool: ${result.error}`);
          return;
        }
      } else {
        // 添加新工具
        const result = toolRegistry.addCustomTool(tool);
        if (!result.success) {
          alert(`Failed to add tool: ${result.error}`);
          return;
        }
      }

      if (onSave) {
        onSave(tool);
      } else {
        alert('✅ Tool saved successfully!');
      }
    } catch (error) {
      alert(`Error: ${error}`);
    }
  };

  const handleFormat = () => {
    try {
      const tool = JSON.parse(jsonText);
      setJsonText(JSON.stringify(tool, null, 2));
    } catch (error) {
      alert('Cannot format invalid JSON');
    }
  };

  const handleLoadTemplate = () => {
    setJsonText(JSON.stringify(TOOL_SCHEMA_EXAMPLE, null, 2));
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100%', backgroundColor: '#0f172a' }}>
      {/* Toolbar */}
      <div
        style={{
          padding: '1rem',
          backgroundColor: '#1f2937',
          borderBottom: '1px solid #374151',
          display: 'flex',
          gap: '0.5rem',
          flexWrap: 'wrap',
        }}
      >
        <button
          onClick={handleValidate}
          style={{
            padding: '0.5rem 1rem',
            backgroundColor: isValid ? '#10b981' : '#f59e0b',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            fontWeight: 600,
          }}
        >
          {isValid ? '✓ Valid' : '⚠️ Validate'}
        </button>
        <button
          onClick={handleFormat}
          style={{
            padding: '0.5rem 1rem',
            backgroundColor: '#3b82f6',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            fontWeight: 600,
          }}
        >
          🎨 Format
        </button>
        <button
          onClick={handleLoadTemplate}
          style={{
            padding: '0.5rem 1rem',
            backgroundColor: '#8b5cf6',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            fontWeight: 600,
          }}
        >
          📄 Load Template
        </button>
        <div style={{ flex: 1 }} />
        <button
          onClick={onCancel}
          style={{
            padding: '0.5rem 1rem',
            backgroundColor: '#6b7280',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            fontWeight: 600,
          }}
        >
          Cancel
        </button>
        <button
          onClick={handleSave}
          style={{
            padding: '0.5rem 1rem',
            backgroundColor: '#059669',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            fontWeight: 600,
          }}
        >
          💾 Save Tool
        </button>
      </div>

      {/* Validation Errors */}
      {validationErrors.length > 0 && (
        <div
          style={{
            padding: '1rem',
            backgroundColor: '#7f1d1d',
            borderBottom: '1px solid #991b1b',
            color: '#fecaca',
            fontSize: '0.875rem',
          }}
        >
          <div style={{ fontWeight: 600, marginBottom: '0.5rem' }}>❌ Validation Errors:</div>
          <ul style={{ margin: 0, paddingLeft: '1.5rem' }}>
            {validationErrors.map((error, idx) => (
              <li key={idx}>{error}</li>
            ))}
          </ul>
        </div>
      )}

      {/* JSON Editor */}
      <div style={{ flex: 1, overflow: 'auto', padding: '1rem' }}>
        <textarea
          value={jsonText}
          onChange={(e) => setJsonText(e.target.value)}
          spellCheck={false}
          style={{
            width: '100%',
            height: '100%',
            minHeight: '500px',
            padding: '1rem',
            backgroundColor: '#111827',
            border: '1px solid #374151',
            borderRadius: '4px',
            color: '#f3f4f6',
            fontSize: '0.875rem',
            fontFamily: 'Monaco, Menlo, "Ubuntu Mono", monospace',
            lineHeight: 1.6,
            resize: 'vertical',
          }}
          placeholder="Enter tool JSON definition here..."
        />
      </div>

      {/* Help Panel */}
      <div
        style={{
          padding: '1rem',
          backgroundColor: '#1f2937',
          borderTop: '1px solid #374151',
          color: '#9ca3af',
          fontSize: '0.75rem',
        }}
      >
        <details>
          <summary style={{ cursor: 'pointer', fontWeight: 600, marginBottom: '0.5rem' }}>
            📖 JSON Schema Help
          </summary>
          <div style={{ lineHeight: 1.6, paddingLeft: '1rem' }}>
            <p><strong>Required Fields:</strong></p>
            <ul style={{ margin: '0.5rem 0', paddingLeft: '1.5rem' }}>
              <li>$schema: "{TOOL_SCHEMA_VERSION}"</li>
              <li>id: Unique identifier (e.g., "author.toolname")</li>
              <li>name: Tool display name</li>
              <li>version: Tool version (e.g., "1.0.0")</li>
              <li>omicsType: genomics, transcriptomics, proteomics, metabolomics, multiomics, general</li>
              <li>category: qc, alignment, peak_calling, etc.</li>
              <li>description: Short description</li>
              <li>parameters: Array of parameter definitions</li>
            </ul>
            <p><strong>Parameter Types:</strong></p>
            <ul style={{ margin: '0.5rem 0', paddingLeft: '1.5rem' }}>
              <li>string, number, integer, boolean</li>
              <li>select, multiselect</li>
              <li>file, directory</li>
              <li>color, date, range</li>
              <li>textarea, json, array</li>
            </ul>
            <p>Click "Load Template" to see a complete example.</p>
          </div>
        </details>
      </div>
    </div>
  );
};

export default ToolJSONEditor;
