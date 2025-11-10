/**
 * Config Panel V3 (通用版)
 * 基于JSON Schema的动态参数配置面板
 * 支持内置工具和用户自定义工具
 */

import { useState, useEffect } from 'react';
import { Node } from 'reactflow';
import { PipelineNodeData } from './PipelineNode';
import { toolRegistry } from '../../utils/toolRegistry';
import { ToolDefinitionSchema } from '../../schemas/ToolSchema';
import DynamicParameterRenderer from './DynamicParameterRenderer';
import {
  getPresetsForTool,
  savePreset,
  deletePreset,
  initializeDefaultPresets,
  type ParameterPreset,
} from '../../utils/parameterPresets';

interface ConfigPanelV3Props {
  selectedNode: Node<PipelineNodeData> | null;
  onUpdate: (nodeId: string, updates: Partial<PipelineNodeData>) => void;
  onClose: () => void;
}

const ConfigPanelV3 = ({ selectedNode, onUpdate, onClose }: ConfigPanelV3Props) => {
  const [label, setLabel] = useState('');
  const [tool, setTool] = useState('');
  const [version, setVersion] = useState('');
  const [parameters, setParameters] = useState<Record<string, any>>({});
  const [toolDefinition, setToolDefinition] = useState<ToolDefinitionSchema | null>(null);
  const [availablePresets, setAvailablePresets] = useState<ParameterPreset[]>([]);
  const [showPresetDialog, setShowPresetDialog] = useState(false);
  const [newPresetName, setNewPresetName] = useState('');
  const [newPresetDescription, setNewPresetDescription] = useState('');
  const [collapsedGroups, setCollapsedGroups] = useState<Set<string>>(new Set());

  useEffect(() => {
    initializeDefaultPresets();
  }, []);

  useEffect(() => {
    if (selectedNode) {
      setLabel(selectedNode.data.label || '');
      setTool(selectedNode.data.tool || '');
      setVersion(selectedNode.data.version || '');
      setParameters(selectedNode.data.parameters || {});

      // 从toolRegistry获取工具定义
      const toolId = selectedNode.id.replace(/^node_\d+_/, '');
      const toolDef = toolRegistry.getTool(toolId);
      
      if (toolDef) {
        setToolDefinition(toolDef);
        
        // 初始化默认参数
        const newParams = { ...selectedNode.data.parameters };
        toolDef.parameters.forEach((param) => {
          if (!(param.name in newParams) && param.default !== undefined) {
            newParams[param.name] = param.default;
          }
        });
        setParameters(newParams);
        
        // 加载预设
        setAvailablePresets(getPresetsForTool(toolDef.id));
      } else {
        setToolDefinition(null);
        setAvailablePresets([]);
      }
    }
  }, [selectedNode]);

  const handleSave = () => {
    if (!selectedNode) return;

    onUpdate(selectedNode.id, {
      label,
      tool,
      version,
      parameters,
    });

    onClose();
  };

  const handleParameterChange = (paramName: string, value: any) => {
    setParameters((prev) => ({
      ...prev,
      [paramName]: value,
    }));
  };

  const handleLoadPreset = (preset: ParameterPreset) => {
    setParameters(preset.parameters);
  };

  const handleSavePreset = () => {
    if (!toolDefinition || !newPresetName.trim()) return;

    savePreset({
      toolId: toolDefinition.id,
      name: newPresetName,
      description: newPresetDescription,
      parameters,
    });

    setAvailablePresets(getPresetsForTool(toolDefinition.id));
    setShowPresetDialog(false);
    setNewPresetName('');
    setNewPresetDescription('');
  };

  const handleDeletePreset = (presetId: string) => {
    if (!confirm('Are you sure you want to delete this preset?')) return;
    
    deletePreset(presetId);
    if (toolDefinition) {
      setAvailablePresets(getPresetsForTool(toolDefinition.id));
    }
  };

  const handleResetDefaults = () => {
    if (!toolDefinition) return;
    
    const defaults: Record<string, any> = {};
    toolDefinition.parameters.forEach((param) => {
      if (param.default !== undefined) {
        defaults[param.name] = param.default;
      }
    });
    setParameters(defaults);
  };

  const toggleGroup = (groupName: string) => {
    setCollapsedGroups((prev) => {
      const newSet = new Set(prev);
      if (newSet.has(groupName)) {
        newSet.delete(groupName);
      } else {
        newSet.add(groupName);
      }
      return newSet;
    });
  };

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

  // 按组分类参数
  const parametersByGroup: Record<string, any[]> = {};
  const ungroupedParameters: any[] = [];

  if (toolDefinition) {
    toolDefinition.parameters.forEach((param) => {
      if (param.group) {
        if (!parametersByGroup[param.group]) {
          parametersByGroup[param.group] = [];
        }
        parametersByGroup[param.group].push(param);
      } else {
        ungroupedParameters.push(param);
      }
    });
  }

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
      <div
        style={{
          padding: '1.5rem',
          borderBottom: '1px solid #374151',
          backgroundColor: '#111827',
        }}
      >
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
          <h3 style={{ fontSize: '1.125rem', fontWeight: 600, color: '#f3f4f6', margin: 0 }}>
            Configure Node
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
          >
            ✕
          </button>
        </div>
      </div>

      {/* Content */}
      <div style={{ flex: 1, overflowY: 'auto', padding: '1.5rem' }}>
        {/* Basic Info */}
        <div style={{ marginBottom: '1.5rem' }}>
          <label style={{ display: 'block', fontSize: '0.875rem', fontWeight: 600, color: '#e5e7eb', marginBottom: '0.5rem' }}>
            Node Label
          </label>
          <input
            type="text"
            value={label}
            onChange={(e) => setLabel(e.target.value)}
            style={{
              width: '100%',
              padding: '0.5rem',
              backgroundColor: '#111827',
              border: '1px solid #374151',
              borderRadius: '4px',
              color: '#f3f4f6',
              fontSize: '0.875rem',
            }}
          />
        </div>

        {/* Tool Info */}
        {toolDefinition && (
          <div
            style={{
              padding: '1rem',
              backgroundColor: '#111827',
              borderRadius: '6px',
              marginBottom: '1.5rem',
              border: '1px solid #374151',
            }}
          >
            <div style={{ fontSize: '0.875rem', fontWeight: 600, color: '#e5e7eb', marginBottom: '0.5rem' }}>
              {toolDefinition.name}
            </div>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.5rem' }}>
              v{toolDefinition.version} • {toolDefinition.omicsType} • {toolDefinition.category}
            </div>
            {toolDefinition.description && (
              <div style={{ fontSize: '0.75rem', color: '#6b7280', lineHeight: 1.5 }}>
                {toolDefinition.description}
              </div>
            )}
            {toolDefinition.documentation && (
              <a
                href={toolDefinition.documentation}
                target="_blank"
                rel="noopener noreferrer"
                style={{
                  display: 'inline-block',
                  marginTop: '0.5rem',
                  fontSize: '0.75rem',
                  color: '#3b82f6',
                  textDecoration: 'none',
                }}
              >
                📖 Documentation →
              </a>
            )}
          </div>
        )}

        {/* Parameter Presets */}
        {availablePresets.length > 0 && (
          <div style={{ marginBottom: '1.5rem' }}>
            <div style={{ fontSize: '0.875rem', fontWeight: 600, color: '#e5e7eb', marginBottom: '0.75rem' }}>
              📋 Parameter Presets
            </div>
            <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
              {availablePresets.map((preset) => (
                <div
                  key={preset.id}
                  style={{
                    padding: '0.75rem',
                    backgroundColor: '#111827',
                    borderRadius: '4px',
                    border: '1px solid #374151',
                  }}
                >
                  <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'start', marginBottom: '0.25rem' }}>
                    <div style={{ fontSize: '0.75rem', fontWeight: 600, color: '#e5e7eb' }}>
                      {preset.name}
                    </div>
                    <div style={{ display: 'flex', gap: '0.25rem' }}>
                      <button
                        onClick={() => handleLoadPreset(preset)}
                        style={{
                          padding: '0.25rem 0.5rem',
                          backgroundColor: '#3b82f6',
                          color: 'white',
                          border: 'none',
                          borderRadius: '3px',
                          cursor: 'pointer',
                          fontSize: '0.7rem',
                        }}
                      >
                        Load
                      </button>
                      <button
                        onClick={() => handleDeletePreset(preset.id!)}
                        style={{
                          padding: '0.25rem 0.5rem',
                          backgroundColor: '#dc2626',
                          color: 'white',
                          border: 'none',
                          borderRadius: '3px',
                          cursor: 'pointer',
                          fontSize: '0.7rem',
                        }}
                      >
                        ✕
                      </button>
                    </div>
                  </div>
                  {preset.description && (
                    <div style={{ fontSize: '0.7rem', color: '#6b7280' }}>
                      {preset.description}
                    </div>
                  )}
                </div>
              ))}
            </div>
          </div>
        )}

        {/* Preset Actions */}
        <div style={{ display: 'flex', gap: '0.5rem', marginBottom: '1.5rem' }}>
          <button
            onClick={() => setShowPresetDialog(true)}
            style={{
              flex: 1,
              padding: '0.5rem',
              backgroundColor: '#8b5cf6',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontSize: '0.75rem',
              fontWeight: 600,
            }}
          >
            💾 Save as Preset
          </button>
          <button
            onClick={handleResetDefaults}
            style={{
              flex: 1,
              padding: '0.5rem',
              backgroundColor: '#374151',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontSize: '0.75rem',
              fontWeight: 600,
            }}
          >
            🔄 Reset Defaults
          </button>
        </div>

        {/* Parameters */}
        {toolDefinition && (
          <div>
            <div style={{ fontSize: '0.875rem', fontWeight: 600, color: '#e5e7eb', marginBottom: '0.75rem' }}>
              ⚙️ Parameters
            </div>

            {/* Ungrouped Parameters */}
            {ungroupedParameters.length > 0 && (
              <div style={{ marginBottom: '1rem' }}>
                {ungroupedParameters.map((param) => (
                  <DynamicParameterRenderer
                    key={param.name}
                    parameter={param}
                    value={parameters[param.name]}
                    onChange={(value) => handleParameterChange(param.name, value)}
                  />
                ))}
              </div>
            )}

            {/* Grouped Parameters */}
            {Object.keys(parametersByGroup).map((groupName) => {
              const isCollapsed = collapsedGroups.has(groupName);
              return (
                <div key={groupName} style={{ marginBottom: '1rem' }}>
                  <button
                    onClick={() => toggleGroup(groupName)}
                    style={{
                      width: '100%',
                      padding: '0.75rem',
                      backgroundColor: '#111827',
                      border: '1px solid #374151',
                      borderRadius: '4px',
                      color: '#e5e7eb',
                      fontSize: '0.875rem',
                      fontWeight: 600,
                      cursor: 'pointer',
                      display: 'flex',
                      justifyContent: 'space-between',
                      alignItems: 'center',
                      marginBottom: '0.5rem',
                    }}
                  >
                    <span>{groupName}</span>
                    <span>{isCollapsed ? '▼' : '▲'}</span>
                  </button>
                  {!isCollapsed && (
                    <div style={{ paddingLeft: '0.5rem' }}>
                      {parametersByGroup[groupName].map((param) => (
                        <DynamicParameterRenderer
                          key={param.name}
                          parameter={param}
                          value={parameters[param.name]}
                          onChange={(value) => handleParameterChange(param.name, value)}
                        />
                      ))}
                    </div>
                  )}
                </div>
              );
            })}
          </div>
        )}
      </div>

      {/* Footer */}
      <div
        style={{
          padding: '1.5rem',
          borderTop: '1px solid #374151',
          backgroundColor: '#111827',
          display: 'flex',
          gap: '0.75rem',
        }}
      >
        <button
          onClick={onClose}
          style={{
            flex: 1,
            padding: '0.75rem',
            backgroundColor: '#374151',
            color: 'white',
            border: 'none',
            borderRadius: '6px',
            cursor: 'pointer',
            fontWeight: 600,
            fontSize: '0.875rem',
          }}
        >
          Cancel
        </button>
        <button
          onClick={handleSave}
          style={{
            flex: 1,
            padding: '0.75rem',
            backgroundColor: '#059669',
            color: 'white',
            border: 'none',
            borderRadius: '6px',
            cursor: 'pointer',
            fontWeight: 600,
            fontSize: '0.875rem',
          }}
        >
          Save Changes
        </button>
      </div>

      {/* Save Preset Dialog */}
      {showPresetDialog && (
        <div
          style={{
            position: 'fixed',
            top: 0,
            left: 0,
            right: 0,
            bottom: 0,
            backgroundColor: 'rgba(0,0,0,0.7)',
            display: 'flex',
            justifyContent: 'center',
            alignItems: 'center',
            zIndex: 3000,
          }}
          onClick={() => setShowPresetDialog(false)}
        >
          <div
            style={{
              backgroundColor: '#1f2937',
              borderRadius: '8px',
              padding: '2rem',
              maxWidth: '400px',
              width: '90%',
              border: '2px solid #374151',
            }}
            onClick={(e) => e.stopPropagation()}
          >
            <h3 style={{ margin: '0 0 1rem 0', color: '#f3f4f6' }}>💾 Save Parameter Preset</h3>
            
            <div style={{ marginBottom: '1rem' }}>
              <label style={{ display: 'block', fontSize: '0.875rem', color: '#e5e7eb', marginBottom: '0.5rem' }}>
                Preset Name *
              </label>
              <input
                type="text"
                value={newPresetName}
                onChange={(e) => setNewPresetName(e.target.value)}
                placeholder="e.g., High Quality QC"
                style={{
                  width: '100%',
                  padding: '0.5rem',
                  backgroundColor: '#111827',
                  border: '1px solid #374151',
                  borderRadius: '4px',
                  color: '#f3f4f6',
                  fontSize: '0.875rem',
                }}
              />
            </div>

            <div style={{ marginBottom: '1.5rem' }}>
              <label style={{ display: 'block', fontSize: '0.875rem', color: '#e5e7eb', marginBottom: '0.5rem' }}>
                Description (optional)
              </label>
              <textarea
                value={newPresetDescription}
                onChange={(e) => setNewPresetDescription(e.target.value)}
                placeholder="Describe this preset..."
                rows={3}
                style={{
                  width: '100%',
                  padding: '0.5rem',
                  backgroundColor: '#111827',
                  border: '1px solid #374151',
                  borderRadius: '4px',
                  color: '#f3f4f6',
                  fontSize: '0.875rem',
                  resize: 'vertical',
                }}
              />
            </div>

            <div style={{ display: 'flex', gap: '0.5rem' }}>
              <button
                onClick={() => setShowPresetDialog(false)}
                style={{
                  flex: 1,
                  padding: '0.75rem',
                  backgroundColor: '#374151',
                  color: 'white',
                  border: 'none',
                  borderRadius: '4px',
                  cursor: 'pointer',
                  fontWeight: 600,
                }}
              >
                Cancel
              </button>
              <button
                onClick={handleSavePreset}
                disabled={!newPresetName.trim()}
                style={{
                  flex: 1,
                  padding: '0.75rem',
                  backgroundColor: newPresetName.trim() ? '#3b82f6' : '#374151',
                  color: 'white',
                  border: 'none',
                  borderRadius: '4px',
                  cursor: newPresetName.trim() ? 'pointer' : 'not-allowed',
                  fontWeight: 600,
                }}
              >
                Save
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default ConfigPanelV3;
