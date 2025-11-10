import { useState, useEffect } from 'react';

interface VisualizationTool {
  tool_id: string;
  name: string;
  category: string;
  runtime: string;
  language: string;
  description: string;
  version: string;
  author: string;
  tags: string[];
  parameters: Record<string, any>;
}

const VisualizationToolsPage = () => {
  const [tools, setTools] = useState<VisualizationTool[]>([]);
  const [selectedTool, setSelectedTool] = useState<VisualizationTool | null>(null);
  const [filterRuntime, setFilterRuntime] = useState<string>('all');
  const [searchQuery, setSearchQuery] = useState('');
  const [showConfigModal, setShowConfigModal] = useState(false);

  // 加载可视化工具
  useEffect(() => {
    // TODO: 从 API 加载
    // 这里使用硬编码的示例数据
    const mockTools: VisualizationTool[] = [
      {
        tool_id: 'ggplot2_scatter',
        name: 'ggplot2 Scatter Plot',
        category: 'visualization',
        runtime: 'r',
        language: 'R',
        description: 'Create publication-ready scatter plots using ggplot2',
        version: '1.0.0',
        author: 'Built-in',
        tags: ['visualization', 'scatter', 'ggplot2', 'R'],
        parameters: {},
      },
      {
        tool_id: 'enriched_heatmap',
        name: 'EnrichedHeatmap',
        category: 'visualization',
        runtime: 'r',
        language: 'R',
        description: 'Create enriched heatmaps for genomic regions using EnrichedHeatmap package',
        version: '1.0.0',
        author: 'Built-in',
        tags: ['visualization', 'heatmap', 'genomics', 'EnrichedHeatmap', 'R'],
        parameters: {},
      },
      {
        tool_id: 'seaborn_heatmap',
        name: 'seaborn Heatmap',
        category: 'visualization',
        runtime: 'python',
        language: 'Python',
        description: 'Create statistical heatmaps using seaborn',
        version: '1.0.0',
        author: 'Built-in',
        tags: ['visualization', 'heatmap', 'seaborn', 'Python'],
        parameters: {},
      },
      {
        tool_id: 'plotly_3d_scatter',
        name: 'Plotly 3D Scatter Plot',
        category: 'visualization',
        runtime: 'python',
        language: 'Python',
        description: 'Create interactive 3D scatter plots using plotly',
        version: '1.0.0',
        author: 'Built-in',
        tags: ['visualization', 'scatter', '3D', 'interactive', 'plotly', 'Python'],
        parameters: {},
      },
      {
        tool_id: 'ggplot2_boxplot',
        name: 'ggplot2 Box Plot',
        category: 'visualization',
        runtime: 'r',
        language: 'R',
        description: 'Create box plots using ggplot2',
        version: '1.0.0',
        author: 'Built-in',
        tags: ['visualization', 'boxplot', 'ggplot2', 'R'],
        parameters: {},
      },
      {
        tool_id: 'pheatmap',
        name: 'pheatmap',
        category: 'visualization',
        runtime: 'r',
        language: 'R',
        description: 'Create annotated heatmaps using pheatmap package',
        version: '1.0.0',
        author: 'Built-in',
        tags: ['visualization', 'heatmap', 'pheatmap', 'R'],
        parameters: {},
      },
      {
        tool_id: 'seaborn_violin',
        name: 'seaborn Violin Plot',
        category: 'visualization',
        runtime: 'python',
        language: 'Python',
        description: 'Create violin plots using seaborn',
        version: '1.0.0',
        author: 'Built-in',
        tags: ['visualization', 'violin', 'seaborn', 'Python'],
        parameters: {},
      },
      {
        tool_id: 'plotly_interactive_heatmap',
        name: 'Plotly Interactive Heatmap',
        category: 'visualization',
        runtime: 'python',
        language: 'Python',
        description: 'Create interactive heatmaps using plotly',
        version: '1.0.0',
        author: 'Built-in',
        tags: ['visualization', 'heatmap', 'interactive', 'plotly', 'Python'],
        parameters: {},
      },
    ];
    setTools(mockTools);
  }, []);

  const filteredTools = tools.filter((tool) => {
    if (filterRuntime !== 'all' && tool.runtime !== filterRuntime) return false;
    if (
      searchQuery &&
      !tool.name.toLowerCase().includes(searchQuery.toLowerCase()) &&
      !tool.description.toLowerCase().includes(searchQuery.toLowerCase()) &&
      !tool.tags.some((tag) => tag.toLowerCase().includes(searchQuery.toLowerCase()))
    ) {
      return false;
    }
    return true;
  });

  const getRuntimeIcon = (runtime: string) => {
    switch (runtime) {
      case 'r':
        return '📊 R';
      case 'python':
        return '🐍 Python';
      case 'binary':
        return '⚙️ Binary';
      default:
        return runtime;
    }
  };

  const getRuntimeColor = (runtime: string) => {
    switch (runtime) {
      case 'r':
        return '#3b82f6'; // blue
      case 'python':
        return '#10b981'; // green
      case 'binary':
        return '#f59e0b'; // amber
      default:
        return '#6b7280'; // gray
    }
  };

  return (
    <div style={{ padding: '2rem', backgroundColor: '#0f172a', minHeight: '100vh', color: '#f3f4f6' }}>
      {/* Header */}
      <div style={{ marginBottom: '2rem' }}>
        <h1 style={{ marginBottom: '0.5rem', color: '#f3f4f6', fontSize: '2rem', fontWeight: 700 }}>
          🎨 Visualization Tools
        </h1>
        <p style={{ color: '#9ca3af', fontSize: '1rem' }}>
          Professional visualization tools from R and Python
        </p>
      </div>

      {/* Filters */}
      <div style={{ marginBottom: '2rem', display: 'flex', gap: '1rem', flexWrap: 'wrap', alignItems: 'center' }}>
        {/* Runtime Filter */}
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
        </div>

        {/* Search */}
        <input
          type="text"
          placeholder="Search tools..."
          value={searchQuery}
          onChange={(e) => setSearchQuery(e.target.value)}
          style={{
            flex: 1,
            minWidth: '200px',
            padding: '0.5rem 1rem',
            backgroundColor: '#374151',
            border: '1px solid #4b5563',
            borderRadius: '6px',
            color: '#f3f4f6',
            fontSize: '0.875rem',
          }}
        />

        {/* Stats */}
        <div style={{ color: '#9ca3af', fontSize: '0.875rem', fontWeight: 600 }}>
          {filteredTools.length} tool{filteredTools.length !== 1 ? 's' : ''}
        </div>
      </div>

      {/* Tools Grid */}
      <div
        style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fill, minmax(350px, 1fr))',
          gap: '1.5rem',
        }}
      >
        {filteredTools.map((tool) => (
          <div
            key={tool.tool_id}
            style={{
              backgroundColor: '#1f2937',
              borderRadius: '12px',
              border: '2px solid #374151',
              padding: '1.5rem',
              cursor: 'pointer',
              transition: 'all 0.2s',
            }}
            onClick={() => setSelectedTool(tool)}
            onMouseEnter={(e) => {
              e.currentTarget.style.borderColor = getRuntimeColor(tool.runtime);
              e.currentTarget.style.transform = 'translateY(-2px)';
            }}
            onMouseLeave={(e) => {
              e.currentTarget.style.borderColor = '#374151';
              e.currentTarget.style.transform = 'translateY(0)';
            }}
          >
            {/* Tool Header */}
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '1rem' }}>
              <div>
                <h3 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.25rem', fontWeight: 700 }}>
                  {tool.name}
                </h3>
                <div style={{ marginTop: '0.5rem', display: 'flex', gap: '0.5rem', alignItems: 'center' }}>
                  <span
                    style={{
                      padding: '0.25rem 0.5rem',
                      backgroundColor: getRuntimeColor(tool.runtime) + '33',
                      color: getRuntimeColor(tool.runtime),
                      borderRadius: '4px',
                      fontSize: '0.75rem',
                      fontWeight: 600,
                    }}
                  >
                    {getRuntimeIcon(tool.runtime)}
                  </span>
                  <span style={{ color: '#9ca3af', fontSize: '0.75rem' }}>v{tool.version}</span>
                </div>
              </div>
            </div>

            {/* Description */}
            <p style={{ color: '#d1d5db', fontSize: '0.875rem', lineHeight: 1.6, marginBottom: '1rem' }}>
              {tool.description}
            </p>

            {/* Tags */}
            <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.5rem', marginBottom: '1rem' }}>
              {tool.tags.slice(0, 4).map((tag) => (
                <span
                  key={tag}
                  style={{
                    padding: '0.25rem 0.5rem',
                    backgroundColor: '#374151',
                    color: '#9ca3af',
                    borderRadius: '4px',
                    fontSize: '0.75rem',
                  }}
                >
                  {tag}
                </span>
              ))}
              {tool.tags.length > 4 && (
                <span
                  style={{
                    padding: '0.25rem 0.5rem',
                    backgroundColor: '#374151',
                    color: '#9ca3af',
                    borderRadius: '4px',
                    fontSize: '0.75rem',
                  }}
                >
                  +{tool.tags.length - 4}
                </span>
              )}
            </div>

            {/* Actions */}
            <div style={{ display: 'flex', gap: '0.5rem' }}>
              <button
                onClick={(e) => {
                  e.stopPropagation();
                  setSelectedTool(tool);
                  setShowConfigModal(true);
                }}
                style={{
                  flex: 1,
                  padding: '0.5rem 1rem',
                  backgroundColor: '#3b82f6',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '0.875rem',
                  fontWeight: 600,
                }}
              >
                Configure
              </button>
              <button
                onClick={(e) => {
                  e.stopPropagation();
                  console.log('Add to pipeline:', tool.tool_id);
                }}
                style={{
                  padding: '0.5rem 1rem',
                  backgroundColor: '#10b981',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '0.875rem',
                  fontWeight: 600,
                }}
                title="Add to Pipeline"
              >
                ➕
              </button>
            </div>
          </div>
        ))}
      </div>

      {/* Tool Detail Modal */}
      {selectedTool && !showConfigModal && (
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
          onClick={() => setSelectedTool(null)}
        >
          <div
            style={{
              backgroundColor: '#1f2937',
              borderRadius: '12px',
              padding: '2rem',
              maxWidth: '800px',
              width: '90%',
              maxHeight: '90vh',
              overflow: 'auto',
              border: '2px solid #374151',
            }}
            onClick={(e) => e.stopPropagation()}
          >
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '1.5rem' }}>
              <div>
                <h2 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.75rem', fontWeight: 700 }}>
                  {selectedTool.name}
                </h2>
                <div style={{ marginTop: '0.75rem', display: 'flex', gap: '0.75rem', alignItems: 'center' }}>
                  <span
                    style={{
                      padding: '0.5rem 1rem',
                      backgroundColor: getRuntimeColor(selectedTool.runtime) + '33',
                      color: getRuntimeColor(selectedTool.runtime),
                      borderRadius: '6px',
                      fontSize: '0.875rem',
                      fontWeight: 600,
                    }}
                  >
                    {getRuntimeIcon(selectedTool.runtime)}
                  </span>
                  <span style={{ color: '#9ca3af', fontSize: '0.875rem' }}>
                    Version {selectedTool.version} • {selectedTool.author}
                  </span>
                </div>
              </div>
              <button
                onClick={() => setSelectedTool(null)}
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

            <div style={{ marginBottom: '1.5rem' }}>
              <h3 style={{ color: '#f3f4f6', fontSize: '1.125rem', marginBottom: '0.5rem' }}>Description</h3>
              <p style={{ color: '#d1d5db', lineHeight: 1.6 }}>{selectedTool.description}</p>
            </div>

            <div style={{ marginBottom: '1.5rem' }}>
              <h3 style={{ color: '#f3f4f6', fontSize: '1.125rem', marginBottom: '0.75rem' }}>Tags</h3>
              <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.5rem' }}>
                {selectedTool.tags.map((tag) => (
                  <span
                    key={tag}
                    style={{
                      padding: '0.5rem 0.75rem',
                      backgroundColor: '#374151',
                      color: '#d1d5db',
                      borderRadius: '6px',
                      fontSize: '0.875rem',
                    }}
                  >
                    {tag}
                  </span>
                ))}
              </div>
            </div>

            <div style={{ display: 'flex', gap: '1rem' }}>
              <button
                onClick={() => setShowConfigModal(true)}
                style={{
                  flex: 1,
                  padding: '0.75rem 1.5rem',
                  backgroundColor: '#3b82f6',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '1rem',
                  fontWeight: 600,
                }}
              >
                Configure & Use
              </button>
              <button
                onClick={() => console.log('Add to pipeline:', selectedTool.tool_id)}
                style={{
                  padding: '0.75rem 1.5rem',
                  backgroundColor: '#10b981',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '1rem',
                  fontWeight: 600,
                }}
              >
                Add to Pipeline
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Configuration Modal */}
      {showConfigModal && selectedTool && (
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
            zIndex: 3001,
          }}
          onClick={() => setShowConfigModal(false)}
        >
          <div
            style={{
              backgroundColor: '#1f2937',
              borderRadius: '12px',
              padding: '2rem',
              maxWidth: '900px',
              width: '90%',
              maxHeight: '90vh',
              overflow: 'auto',
              border: '2px solid #374151',
            }}
            onClick={(e) => e.stopPropagation()}
          >
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1.5rem' }}>
              <h2 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.75rem', fontWeight: 700 }}>
                Configure: {selectedTool.name}
              </h2>
              <button
                onClick={() => setShowConfigModal(false)}
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

            <div style={{ padding: '2rem', backgroundColor: '#111827', borderRadius: '8px', textAlign: 'center' }}>
              <p style={{ color: '#9ca3af', fontSize: '1rem' }}>
                Parameter configuration interface coming soon...
              </p>
              <p style={{ color: '#6b7280', fontSize: '0.875rem', marginTop: '0.5rem' }}>
                This will allow you to configure all tool parameters with a visual interface.
              </p>
            </div>

            <div style={{ marginTop: '1.5rem', display: 'flex', gap: '1rem' }}>
              <button
                onClick={() => console.log('Save configuration')}
                style={{
                  flex: 1,
                  padding: '0.75rem 1.5rem',
                  backgroundColor: '#10b981',
                  color: 'white',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '1rem',
                  fontWeight: 600,
                }}
              >
                Save & Add to Pipeline
              </button>
              <button
                onClick={() => setShowConfigModal(false)}
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
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default VisualizationToolsPage;
