/**
 * Tool Manager Page
 * 工具管理页面 - 浏览、上传、编辑、删除自定义工具
 * 
 * 功能：
 * 1. 显示所有工具（内置+自定义）
 * 2. 上传自定义工具JSON
 * 3. 在线编辑工具定义
 * 4. 导出工具为JSON
 * 5. 启用/禁用工具
 * 6. 搜索和过滤
 */

import React, { useState, useEffect } from 'react';
import { toolRegistry } from '../utils/toolRegistry';
import { ToolDefinitionSchema, TOOL_SCHEMA_EXAMPLE } from '../schemas/ToolSchema';

const ToolManagerPage: React.FC = () => {
  const [tools, setTools] = useState<ToolDefinitionSchema[]>([]);
  const [searchQuery, setSearchQuery] = useState('');
  const [filterOmics, setFilterOmics] = useState<string>('all');
  const [filterCategory, setFilterCategory] = useState<string>('all');
  const [showOnlyCustom, setShowOnlyCustom] = useState(false);
  const [selectedTool, setSelectedTool] = useState<ToolDefinitionSchema | null>(null);
  const [showEditor, setShowEditor] = useState(false);
  const [showUpload, setShowUpload] = useState(false);
  const [stats, setStats] = useState(toolRegistry.getStats());

  useEffect(() => {
    loadTools();
  }, [filterOmics, filterCategory, showOnlyCustom]);

  const loadTools = () => {
    let filteredTools = toolRegistry.getAllTools();

    if (showOnlyCustom) {
      filteredTools = filteredTools.filter((t) => t.id.startsWith('custom.') || t.id.includes('.'));
    }

    if (filterOmics !== 'all') {
      filteredTools = filteredTools.filter((t) => t.omicsType === filterOmics);
    }

    if (filterCategory !== 'all') {
      filteredTools = filteredTools.filter((t) => t.category === filterCategory);
    }

    if (searchQuery) {
      filteredTools = toolRegistry.searchTools(searchQuery);
    }

    setTools(filteredTools);
    setStats(toolRegistry.getStats());
  };

  const handleUploadJSON = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    const result = await toolRegistry.importToolFromJSON(file);
    if (result.success) {
      alert(`✅ Successfully imported tool: ${result.tool?.name || 'Unknown'}`);
      loadTools();
      setShowUpload(false);
    } else {
      alert(`❌ Failed to import tool: ${result.error}`);
    }

    // Reset input
    e.target.value = '';
  };

  const handleDeleteTool = (id: string) => {
    if (!confirm(`Are you sure you want to delete tool "${id}"?`)) return;

    const result = toolRegistry.deleteCustomTool(id);
    if (result.success) {
      alert('✅ Tool deleted successfully');
      loadTools();
      setSelectedTool(null);
    } else {
      alert(`❌ Failed to delete tool: ${result.error}`);
    }
  };

  const handleToggleEnabled = (id: string, enabled: boolean) => {
    toolRegistry.setToolEnabled(id, enabled);
    loadTools();
  };

  const handleExportTool = (id: string) => {
    toolRegistry.downloadToolJSON(id);
  };

  const handleExportAll = () => {
    toolRegistry.downloadAllCustomTools();
  };

  const omicsTypes = ['all', 'genomics', 'transcriptomics', 'proteomics', 'metabolomics', 'multiomics', 'general'];
  const categories = ['all', 'input', 'qc', 'alignment', 'peak_calling', 'variant_calling', 'feature_extraction', 'differential_expression', 'normalization', 'data_processing', 'statistical_analysis', 'integration', 'visualization', 'output'];

  return (
    <div style={{ padding: '2rem', backgroundColor: '#0f172a', minHeight: '100vh', color: '#f3f4f6' }}>
      {/* Header */}
      <div style={{ marginBottom: '2rem' }}>
        <h1 style={{ fontSize: '2rem', fontWeight: 'bold', marginBottom: '0.5rem' }}>🔧 Tool Manager</h1>
        <p style={{ color: '#9ca3af' }}>Manage built-in and custom bioinformatics tools</p>
      </div>

      {/* Stats */}
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '1rem', marginBottom: '2rem' }}>
        <div style={{ padding: '1rem', backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151' }}>
          <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>Total Tools</div>
          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#3b82f6' }}>{stats.total}</div>
        </div>
        <div style={{ padding: '1rem', backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151' }}>
          <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>Built-in</div>
          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#10b981' }}>{stats.builtIn}</div>
        </div>
        <div style={{ padding: '1rem', backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151' }}>
          <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>Custom</div>
          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#8b5cf6' }}>{stats.custom}</div>
        </div>
        <div style={{ padding: '1rem', backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151' }}>
          <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>Enabled</div>
          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#059669' }}>{stats.enabled}</div>
        </div>
      </div>

      {/* Actions */}
      <div style={{ display: 'flex', gap: '1rem', marginBottom: '2rem', flexWrap: 'wrap' }}>
        <button
          onClick={() => setShowUpload(true)}
          style={{
            padding: '0.75rem 1.5rem',
            backgroundColor: '#3b82f6',
            color: 'white',
            border: 'none',
            borderRadius: '6px',
            cursor: 'pointer',
            fontWeight: 600,
            fontSize: '0.875rem',
          }}
        >
          📤 Upload Tool JSON
        </button>
        <button
          onClick={() => {
            setSelectedTool(TOOL_SCHEMA_EXAMPLE);
            setShowEditor(true);
          }}
          style={{
            padding: '0.75rem 1.5rem',
            backgroundColor: '#8b5cf6',
            color: 'white',
            border: 'none',
            borderRadius: '6px',
            cursor: 'pointer',
            fontWeight: 600,
            fontSize: '0.875rem',
          }}
        >
          ➕ Create New Tool
        </button>
        <button
          onClick={handleExportAll}
          disabled={stats.custom === 0}
          style={{
            padding: '0.75rem 1.5rem',
            backgroundColor: stats.custom > 0 ? '#059669' : '#374151',
            color: 'white',
            border: 'none',
            borderRadius: '6px',
            cursor: stats.custom > 0 ? 'pointer' : 'not-allowed',
            fontWeight: 600,
            fontSize: '0.875rem',
          }}
        >
          💾 Export All Custom Tools
        </button>
      </div>

      {/* Filters */}
      <div style={{ marginBottom: '2rem', display: 'flex', gap: '1rem', flexWrap: 'wrap', alignItems: 'center' }}>
        <input
          type="text"
          placeholder="🔍 Search tools..."
          value={searchQuery}
          onChange={(e) => {
            setSearchQuery(e.target.value);
            loadTools();
          }}
          style={{
            flex: 1,
            minWidth: '200px',
            padding: '0.75rem',
            backgroundColor: '#1f2937',
            border: '1px solid #374151',
            borderRadius: '6px',
            color: '#f3f4f6',
            fontSize: '0.875rem',
          }}
        />
        <select
          value={filterOmics}
          onChange={(e) => setFilterOmics(e.target.value)}
          style={{
            padding: '0.75rem',
            backgroundColor: '#1f2937',
            border: '1px solid #374151',
            borderRadius: '6px',
            color: '#f3f4f6',
            fontSize: '0.875rem',
          }}
        >
          {omicsTypes.map((type) => (
            <option key={type} value={type}>
              {type === 'all' ? 'All Omics' : type}
            </option>
          ))}
        </select>
        <select
          value={filterCategory}
          onChange={(e) => setFilterCategory(e.target.value)}
          style={{
            padding: '0.75rem',
            backgroundColor: '#1f2937',
            border: '1px solid #374151',
            borderRadius: '6px',
            color: '#f3f4f6',
            fontSize: '0.875rem',
          }}
        >
          {categories.map((cat) => (
            <option key={cat} value={cat}>
              {cat === 'all' ? 'All Categories' : cat}
            </option>
          ))}
        </select>
        <label style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', cursor: 'pointer' }}>
          <input
            type="checkbox"
            checked={showOnlyCustom}
            onChange={(e) => setShowOnlyCustom(e.target.checked)}
            style={{ width: '1rem', height: '1rem', cursor: 'pointer' }}
          />
          <span style={{ fontSize: '0.875rem' }}>Custom Only</span>
        </label>
      </div>

      {/* Tools Grid */}
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(300px, 1fr))', gap: '1rem' }}>
        {tools.map((tool) => {
          const isCustom = !tool.id.startsWith('legacy.');
          const isEnabled = toolRegistry.isToolEnabled(tool.id);

          return (
            <div
              key={tool.id}
              style={{
                padding: '1.5rem',
                backgroundColor: '#1f2937',
                border: `2px solid ${isEnabled ? '#374151' : '#7f1d1d'}`,
                borderRadius: '8px',
                opacity: isEnabled ? 1 : 0.6,
              }}
            >
              {/* Header */}
              <div style={{ marginBottom: '1rem' }}>
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'start', marginBottom: '0.5rem' }}>
                  <h3 style={{ fontSize: '1.125rem', fontWeight: 600, color: '#f3f4f6', margin: 0 }}>
                    {tool.name}
                  </h3>
                  {isCustom && (
                    <span
                      style={{
                        padding: '0.25rem 0.5rem',
                        backgroundColor: '#8b5cf6',
                        color: 'white',
                        borderRadius: '4px',
                        fontSize: '0.75rem',
                        fontWeight: 600,
                      }}
                    >
                      CUSTOM
                    </span>
                  )}
                </div>
                <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>
                  {tool.omicsType} • {tool.category} • v{tool.version}
                </div>
              </div>

              {/* Description */}
              <p style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '1rem', lineHeight: 1.5 }}>
                {tool.description}
              </p>

              {/* Parameters Count */}
              <div style={{ fontSize: '0.875rem', color: '#6b7280', marginBottom: '1rem' }}>
                📋 {tool.parameters.length} parameters
              </div>

              {/* Actions */}
              <div style={{ display: 'flex', gap: '0.5rem', flexWrap: 'wrap' }}>
                <button
                  onClick={() => {
                    setSelectedTool(tool);
                    setShowEditor(true);
                  }}
                  style={{
                    padding: '0.5rem 1rem',
                    backgroundColor: '#3b82f6',
                    color: 'white',
                    border: 'none',
                    borderRadius: '4px',
                    cursor: 'pointer',
                    fontSize: '0.75rem',
                    fontWeight: 600,
                  }}
                >
                  View
                </button>
                <button
                  onClick={() => handleExportTool(tool.id)}
                  style={{
                    padding: '0.5rem 1rem',
                    backgroundColor: '#059669',
                    color: 'white',
                    border: 'none',
                    borderRadius: '4px',
                    cursor: 'pointer',
                    fontSize: '0.75rem',
                    fontWeight: 600,
                  }}
                >
                  Export
                </button>
                <button
                  onClick={() => handleToggleEnabled(tool.id, !isEnabled)}
                  style={{
                    padding: '0.5rem 1rem',
                    backgroundColor: isEnabled ? '#f59e0b' : '#10b981',
                    color: 'white',
                    border: 'none',
                    borderRadius: '4px',
                    cursor: 'pointer',
                    fontSize: '0.75rem',
                    fontWeight: 600,
                  }}
                >
                  {isEnabled ? 'Disable' : 'Enable'}
                </button>
                {isCustom && (
                  <button
                    onClick={() => handleDeleteTool(tool.id)}
                    style={{
                      padding: '0.5rem 1rem',
                      backgroundColor: '#dc2626',
                      color: 'white',
                      border: 'none',
                      borderRadius: '4px',
                      cursor: 'pointer',
                      fontSize: '0.75rem',
                      fontWeight: 600,
                    }}
                  >
                    Delete
                  </button>
                )}
              </div>
            </div>
          );
        })}
      </div>

      {/* No Results */}
      {tools.length === 0 && (
        <div style={{ textAlign: 'center', padding: '4rem', color: '#6b7280' }}>
          <div style={{ fontSize: '3rem', marginBottom: '1rem' }}>🔍</div>
          <div style={{ fontSize: '1.25rem', fontWeight: 600, marginBottom: '0.5rem' }}>No tools found</div>
          <div style={{ fontSize: '0.875rem' }}>Try adjusting your filters or search query</div>
        </div>
      )}

      {/* Upload Modal */}
      {showUpload && (
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
            zIndex: 2000,
          }}
          onClick={() => setShowUpload(false)}
        >
          <div
            style={{
              backgroundColor: '#1f2937',
              borderRadius: '8px',
              padding: '2rem',
              maxWidth: '500px',
              width: '90%',
              border: '2px solid #374151',
            }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 style={{ margin: '0 0 1rem 0', color: '#f3f4f6' }}>📤 Upload Tool JSON</h2>
            <p style={{ color: '#9ca3af', fontSize: '0.875rem', marginBottom: '1.5rem' }}>
              Select a JSON file containing a tool definition or tool collection.
            </p>
            <input
              type="file"
              accept=".json"
              onChange={handleUploadJSON}
              style={{
                width: '100%',
                padding: '0.75rem',
                backgroundColor: '#111827',
                border: '1px solid #374151',
                borderRadius: '4px',
                color: '#f3f4f6',
                fontSize: '0.875rem',
                marginBottom: '1rem',
              }}
            />
            <div style={{ display: 'flex', gap: '0.5rem' }}>
              <button
                onClick={() => setShowUpload(false)}
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
            </div>
          </div>
        </div>
      )}

      {/* Editor Modal (View/Edit) */}
      {showEditor && selectedTool && (
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
            zIndex: 2000,
            padding: '2rem',
          }}
          onClick={() => setShowEditor(false)}
        >
          <div
            style={{
              backgroundColor: '#1f2937',
              borderRadius: '8px',
              padding: '2rem',
              maxWidth: '800px',
              width: '100%',
              maxHeight: '90vh',
              overflow: 'auto',
              border: '2px solid #374151',
            }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 style={{ margin: '0 0 1rem 0', color: '#f3f4f6' }}>🔧 {selectedTool.name}</h2>
            
            <div style={{ marginBottom: '1rem' }}>
              <pre
                style={{
                  backgroundColor: '#111827',
                  padding: '1rem',
                  borderRadius: '4px',
                  overflow: 'auto',
                  fontSize: '0.75rem',
                  color: '#e5e7eb',
                  fontFamily: 'monospace',
                  maxHeight: '60vh',
                }}
              >
                {JSON.stringify(selectedTool, null, 2)}
              </pre>
            </div>

            <button
              onClick={() => setShowEditor(false)}
              style={{
                padding: '0.75rem 1.5rem',
                backgroundColor: '#3b82f6',
                color: 'white',
                border: 'none',
                borderRadius: '4px',
                cursor: 'pointer',
                fontWeight: 600,
                width: '100%',
              }}
            >
              Close
            </button>
          </div>
        </div>
      )}
    </div>
  );
};

export default ToolManagerPage;
