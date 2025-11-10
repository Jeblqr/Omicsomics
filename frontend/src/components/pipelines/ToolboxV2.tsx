import { useState } from 'react';
import {
  ToolDefinition,
  allTools,
  getOmicsTypes,
  getCategoriesByOmicsType,
  getToolsByOmicsType,
} from './ToolDefinitions';

interface ToolboxV2Props {
  onDragStart: (tool: ToolDefinition) => void;
}

const ToolboxV2 = ({ onDragStart }: ToolboxV2Props) => {
  const [selectedOmicsType, setSelectedOmicsType] = useState<string>('All');
  const [selectedCategory, setSelectedCategory] = useState<string>('All');
  const [searchQuery, setSearchQuery] = useState('');
  const [expandedTool, setExpandedTool] = useState<string | null>(null);

  const omicsTypes = ['All', ...getOmicsTypes()];
  const categories =
    selectedOmicsType === 'All'
      ? ['All', ...Array.from(new Set(allTools.map((t) => t.operationCategory)))]
      : ['All', ...getCategoriesByOmicsType(selectedOmicsType)];

  // Filter tools based on selections
  const filteredTools = allTools.filter((tool) => {
    const matchesOmics =
      selectedOmicsType === 'All' ||
      tool.omicsType === selectedOmicsType ||
      tool.omicsType === 'General';
    const matchesCategory =
      selectedCategory === 'All' || tool.operationCategory === selectedCategory;
    const matchesSearch =
      searchQuery === '' ||
      tool.name.toLowerCase().includes(searchQuery.toLowerCase()) ||
      tool.description.toLowerCase().includes(searchQuery.toLowerCase()) ||
      tool.tool.toLowerCase().includes(searchQuery.toLowerCase());
    return matchesOmics && matchesCategory && matchesSearch;
  });

  // Group tools by category for display
  const groupedTools = filteredTools.reduce((acc, tool) => {
    const category = tool.operationCategory;
    if (!acc[category]) {
      acc[category] = [];
    }
    acc[category].push(tool);
    return acc;
  }, {} as Record<string, ToolDefinition[]>);

  const handleDragStart = (event: React.DragEvent, tool: ToolDefinition) => {
    event.dataTransfer.effectAllowed = 'move';
    event.dataTransfer.setData('application/reactflow', JSON.stringify(tool));
    onDragStart(tool);
  };

  const getCategoryIcon = (category: string): string => {
    const icons: Record<string, string> = {
      Input: '📥',
      QC: '🔬',
      'Quality Control': '🔬',
      Alignment: '🧬',
      'Peak Calling': '⛰️',
      'Variant Calling': '🔍',
      'Feature Extraction': '📊',
      'Differential Expression': '📈',
      Normalization: '⚖️',
      'Data Processing': '⚙️',
      'Statistical Analysis': '📉',
      Integration: '🔗',
      Visualization: '📊',
      Output: '📤',
    };
    return icons[category] || '🔧';
  };

  const getOmicsIcon = (omicsType: string): string => {
    const icons: Record<string, string> = {
      Genomics: '🧬',
      Transcriptomics: '📜',
      Proteomics: '🧪',
      Metabolomics: '⚗️',
      'Multi-omics': '🔬',
      General: '⚙️',
    };
    return icons[omicsType] || '🔧';
  };

  return (
    <div
      style={{
        width: '320px',
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

        {/* Level 1: Omics Type Filter */}
        <div style={{ marginBottom: '0.75rem' }}>
          <label
            style={{
              display: 'block',
              fontSize: '0.75rem',
              color: '#9ca3af',
              marginBottom: '0.25rem',
              fontWeight: 600,
              textTransform: 'uppercase',
            }}
          >
            Omics Type
          </label>
          <select
            value={selectedOmicsType}
            onChange={(e) => {
              setSelectedOmicsType(e.target.value);
              setSelectedCategory('All');
            }}
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
            {omicsTypes.map((type) => (
              <option key={type} value={type}>
                {getOmicsIcon(type)} {type}
              </option>
            ))}
          </select>
        </div>

        {/* Level 2: Operation Category Filter */}
        <div>
          <label
            style={{
              display: 'block',
              fontSize: '0.75rem',
              color: '#9ca3af',
              marginBottom: '0.25rem',
              fontWeight: 600,
              textTransform: 'uppercase',
            }}
          >
            Operation
          </label>
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
                {cat !== 'All' && getCategoryIcon(cat)} {cat}
              </option>
            ))}
          </select>
        </div>
      </div>

      {/* Tools List - Level 3: Specific Tools */}
      <div style={{ flex: 1, overflowY: 'auto', padding: '0.5rem' }}>
        {filteredTools.length === 0 ? (
          <div
            style={{
              padding: '2rem 1rem',
              textAlign: 'center',
              color: '#9ca3af',
              fontSize: '0.875rem',
            }}
          >
            No tools found
          </div>
        ) : (
          <>
            {Object.entries(groupedTools).map(([category, tools]) => (
              <div key={category} style={{ marginBottom: '1rem' }}>
                {/* Category Header */}
                <div
                  style={{
                    padding: '0.5rem 0.75rem',
                    backgroundColor: '#111827',
                    borderRadius: '4px',
                    marginBottom: '0.5rem',
                    borderLeft: '3px solid #3b82f6',
                  }}
                >
                  <div style={{ fontSize: '0.8rem', fontWeight: 600, color: '#60a5fa' }}>
                    {getCategoryIcon(category)} {category}
                  </div>
                  <div style={{ fontSize: '0.7rem', color: '#6b7280' }}>
                    {tools.length} tool{tools.length > 1 ? 's' : ''}
                  </div>
                </div>

                {/* Tools in Category */}
                {tools.map((tool) => (
                  <div key={tool.id} style={{ marginBottom: '0.5rem' }}>
                    <div
                      draggable
                      onDragStart={(e) => handleDragStart(e, tool)}
                      onClick={() => setExpandedTool(expandedTool === tool.id ? null : tool.id)}
                      style={{
                        padding: '0.75rem',
                        backgroundColor: '#374151',
                        borderRadius: '6px',
                        cursor: 'grab',
                        border: expandedTool === tool.id ? '2px solid #3b82f6' : '1px solid #4b5563',
                        transition: 'all 0.2s',
                      }}
                      onMouseEnter={(e) => {
                        if (expandedTool !== tool.id) {
                          e.currentTarget.style.backgroundColor = '#4b5563';
                          e.currentTarget.style.transform = 'translateX(4px)';
                        }
                      }}
                      onMouseLeave={(e) => {
                        if (expandedTool !== tool.id) {
                          e.currentTarget.style.backgroundColor = '#374151';
                          e.currentTarget.style.transform = 'translateX(0)';
                        }
                      }}
                    >
                      <div style={{ display: 'flex', alignItems: 'flex-start', gap: '0.5rem' }}>
                        <span style={{ fontSize: '1rem', marginTop: '0.1rem' }}>
                          {tool.type === 'input' && '📥'}
                          {tool.type === 'process' && '⚙️'}
                          {tool.type === 'filter' && '🔍'}
                          {tool.type === 'transform' && '🔄'}
                          {tool.type === 'analysis' && '📊'}
                          {tool.type === 'output' && '📤'}
                          {tool.type === 'visualization' && '📈'}
                        </span>
                        <div style={{ flex: 1 }}>
                          <div
                            style={{
                              fontWeight: 500,
                              color: '#f3f4f6',
                              fontSize: '0.875rem',
                              marginBottom: '0.25rem',
                            }}
                          >
                            {tool.name}
                          </div>
                          <div style={{ fontSize: '0.7rem', color: '#9ca3af', marginBottom: '0.25rem' }}>
                            {tool.tool} v{tool.version}
                          </div>
                          {tool.omicsType !== 'General' && (
                            <div
                              style={{
                                display: 'inline-block',
                                fontSize: '0.65rem',
                                padding: '0.15rem 0.4rem',
                                backgroundColor: '#1f2937',
                                color: '#60a5fa',
                                borderRadius: '3px',
                                marginBottom: '0.25rem',
                              }}
                            >
                              {getOmicsIcon(tool.omicsType)} {tool.omicsType}
                            </div>
                          )}
                          <div style={{ fontSize: '0.75rem', color: '#9ca3af', lineHeight: '1.3' }}>
                            {tool.description}
                          </div>

                          {/* Expanded Details */}
                          {expandedTool === tool.id && (
                            <div
                              style={{
                                marginTop: '0.75rem',
                                paddingTop: '0.75rem',
                                borderTop: '1px solid #4b5563',
                              }}
                            >
                              <div
                                style={{
                                  fontSize: '0.7rem',
                                  fontWeight: 600,
                                  color: '#9ca3af',
                                  marginBottom: '0.5rem',
                                  textTransform: 'uppercase',
                                }}
                              >
                                Parameters ({tool.parameterTemplate.length})
                              </div>
                              <div style={{ fontSize: '0.7rem', color: '#9ca3af', lineHeight: '1.4' }}>
                                {tool.parameterTemplate.slice(0, 5).map((param) => (
                                  <div
                                    key={param.name}
                                    style={{
                                      marginBottom: '0.25rem',
                                      display: 'flex',
                                      gap: '0.25rem',
                                    }}
                                  >
                                    <span style={{ color: '#60a5fa' }}>•</span>
                                    <span>
                                      {param.label}
                                      {param.required && (
                                        <span style={{ color: '#ef4444' }}> *</span>
                                      )}
                                    </span>
                                  </div>
                                ))}
                                {tool.parameterTemplate.length > 5 && (
                                  <div style={{ color: '#6b7280', fontStyle: 'italic' }}>
                                    +{tool.parameterTemplate.length - 5} more...
                                  </div>
                                )}
                              </div>
                              <div
                                style={{
                                  marginTop: '0.5rem',
                                  fontSize: '0.65rem',
                                  color: '#6b7280',
                                  fontStyle: 'italic',
                                }}
                              >
                                Click to collapse • Drag to add to pipeline
                              </div>
                            </div>
                          )}
                        </div>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            ))}
          </>
        )}
      </div>

      {/* Footer Stats */}
      <div
        style={{
          padding: '0.75rem',
          borderTop: '1px solid #374151',
          backgroundColor: '#111827',
          fontSize: '0.7rem',
          color: '#6b7280',
        }}
      >
        <div style={{ marginBottom: '0.25rem' }}>
          <strong style={{ color: '#9ca3af' }}>{filteredTools.length}</strong> tools available
        </div>
        <div style={{ fontSize: '0.65rem', color: '#6b7280' }}>
          Drag tools to canvas to build pipeline
        </div>
      </div>
    </div>
  );
};

export default ToolboxV2;
