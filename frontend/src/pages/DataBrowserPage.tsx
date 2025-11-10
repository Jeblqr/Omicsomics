import { useState, useEffect } from 'react';
import { useProjects } from '../contexts/ProjectsContext';
import QuickVisualizerModal from '../components/QuickVisualizerModal';

interface DataFile {
  id: string;
  name: string;
  type: string;
  size: string;
  category: 'raw_data' | 'run_output' | 'dataset' | 'archived';
  created_at: string;
  source?: string; // Run ID or upload source
  tags?: string[];
}

type ViewCategory = 'all' | 'raw_data' | 'run_outputs' | 'datasets' | 'archived';

const DataBrowserPage = () => {
  const { selectedProject } = useProjects();
  const [files, setFiles] = useState<DataFile[]>([]);
  const [selectedCategory, setSelectedCategory] = useState<ViewCategory>('all');
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedFiles, setSelectedFiles] = useState<Set<string>>(new Set());
  const [showPreview, setShowPreview] = useState(false);
  const [previewFile, setPreviewFile] = useState<DataFile | null>(null);
  const [showVisualizer, setShowVisualizer] = useState(false);
  const [visualizerFile, setVisualizerFile] = useState<DataFile | null>(null);

  // 模拟数据
  useEffect(() => {
    const mockFiles: DataFile[] = [
      {
        id: '1',
        name: 'proteomics_raw_batch1.csv',
        type: 'CSV',
        size: '2.3 MB',
        category: 'raw_data',
        created_at: '2025-01-08 10:30',
        tags: ['proteomics', 'batch1'],
      },
      {
        id: '2',
        name: 'qc_report.html',
        type: 'HTML',
        size: '156 KB',
        category: 'run_output',
        created_at: '2025-01-09 14:20',
        source: 'Run #123',
      },
      {
        id: '3',
        name: 'gene_expression.tsv',
        type: 'TSV',
        size: '5.1 MB',
        category: 'run_output',
        created_at: '2025-01-09 15:45',
        source: 'Run #124',
        tags: ['RNA-seq', 'expression'],
      },
    ];
    setFiles(mockFiles);
  }, [selectedProject]);

  const filteredFiles = files.filter((file) => {
    if (selectedCategory !== 'all' && file.category !== selectedCategory) return false;
    if (searchQuery && !file.name.toLowerCase().includes(searchQuery.toLowerCase())) return false;
    return true;
  });

  const stats = {
    totalFiles: files.length,
    totalSize: '15.3 GB',
    rawData: files.filter((f) => f.category === 'raw_data').length,
    runOutputs: files.filter((f) => f.category === 'run_output').length,
    datasets: files.filter((f) => f.category === 'dataset').length,
    archived: files.filter((f) => f.category === 'archived').length,
  };

  const handleSelectFile = (fileId: string) => {
    const newSelected = new Set(selectedFiles);
    if (newSelected.has(fileId)) {
      newSelected.delete(fileId);
    } else {
      newSelected.add(fileId);
    }
    setSelectedFiles(newSelected);
  };

  const handleSelectAll = () => {
    if (selectedFiles.size === filteredFiles.length) {
      setSelectedFiles(new Set());
    } else {
      setSelectedFiles(new Set(filteredFiles.map((f) => f.id)));
    }
  };

  const handlePreview = (file: DataFile) => {
    setPreviewFile(file);
    setShowPreview(true);
  };

  const handleDownload = (fileId: string) => {
    console.log('Download file:', fileId);
    // TODO: 实现下载功能
  };

  const handleDelete = (fileId: string) => {
    if (window.confirm('Are you sure you want to move this file to archive?')) {
      console.log('Archive file:', fileId);
      // TODO: 实现归档功能
    }
  };

  const handleVisualize = (file: DataFile) => {
    setVisualizerFile(file);
    setShowVisualizer(true);
  };

  if (!selectedProject) {
    return (
      <div style={{ padding: '2rem', backgroundColor: '#0f172a', minHeight: '100vh' }}>
        <div
          style={{
            padding: '2rem',
            backgroundColor: '#422006',
            borderRadius: '8px',
            border: '1px solid #78350f',
            color: '#fcd34d',
          }}
        >
          <strong style={{ color: '#fbbf24' }}>⚠️ No project selected</strong>
          <p style={{ marginBottom: 0, marginTop: '0.5rem', color: '#fcd34d' }}>
            Please select a project from the Projects page to view its data.
          </p>
        </div>
      </div>
    );
  }

  return (
    <div style={{ padding: '2rem', backgroundColor: '#0f172a', minHeight: '100vh', color: '#f3f4f6' }}>
      {/* Header */}
      <div style={{ marginBottom: '2rem' }}>
        <h1 style={{ marginBottom: '0.5rem', color: '#f3f4f6', fontSize: '2rem', fontWeight: 700 }}>
          📁 Data Browser
        </h1>
        <p style={{ color: '#9ca3af', fontSize: '1rem' }}>
          Browse, manage, and organize all project data
        </p>
      </div>

      {/* Statistics Cards */}
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '1rem', marginBottom: '2rem' }}>
        <div style={{ padding: '1.5rem', backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151' }}>
          <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem', fontWeight: 600 }}>
            Total Files
          </div>
          <div style={{ fontSize: '2rem', color: '#3b82f6', fontWeight: 700 }}>{stats.totalFiles}</div>
          <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginTop: '0.25rem' }}>{stats.totalSize}</div>
        </div>

        <div style={{ padding: '1.5rem', backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151' }}>
          <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem', fontWeight: 600 }}>
            📥 Raw Data
          </div>
          <div style={{ fontSize: '2rem', color: '#10b981', fontWeight: 700 }}>{stats.rawData}</div>
        </div>

        <div style={{ padding: '1.5rem', backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151' }}>
          <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem', fontWeight: 600 }}>
            🔬 Run Outputs
          </div>
          <div style={{ fontSize: '2rem', color: '#8b5cf6', fontWeight: 700 }}>{stats.runOutputs}</div>
        </div>

        <div style={{ padding: '1.5rem', backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151' }}>
          <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem', fontWeight: 600 }}>
            📊 Datasets
          </div>
          <div style={{ fontSize: '2rem', color: '#f59e0b', fontWeight: 700 }}>{stats.datasets}</div>
        </div>

        <div style={{ padding: '1.5rem', backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151' }}>
          <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem', fontWeight: 600 }}>
            🗑️ Archived
          </div>
          <div style={{ fontSize: '2rem', color: '#6b7280', fontWeight: 700 }}>{stats.archived}</div>
        </div>
      </div>

      {/* Controls */}
      <div style={{ marginBottom: '1.5rem', display: 'flex', gap: '1rem', flexWrap: 'wrap' }}>
        {/* Category Filters */}
        <div style={{ display: 'flex', gap: '0.5rem' }}>
          {(['all', 'raw_data', 'run_outputs', 'datasets', 'archived'] as ViewCategory[]).map((category) => (
            <button
              key={category}
              onClick={() => setSelectedCategory(category)}
              style={{
                padding: '0.5rem 1rem',
                backgroundColor: selectedCategory === category ? '#3b82f6' : '#374151',
                color: '#f3f4f6',
                border: 'none',
                borderRadius: '6px',
                cursor: 'pointer',
                fontSize: '0.875rem',
                fontWeight: 600,
              }}
            >
              {category === 'all' && 'All Files'}
              {category === 'raw_data' && '📥 Raw Data'}
              {category === 'run_outputs' && '🔬 Run Outputs'}
              {category === 'datasets' && '📊 Datasets'}
              {category === 'archived' && '🗑️ Archived'}
            </button>
          ))}
        </div>

        {/* Search */}
        <input
          type="text"
          placeholder="Search files..."
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

        {/* Batch Actions */}
        {selectedFiles.size > 0 && (
          <div style={{ display: 'flex', gap: '0.5rem' }}>
            <button
              onClick={() => console.log('Batch download')}
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
            >
              💾 Download ({selectedFiles.size})
            </button>
            <button
              onClick={() => console.log('Create dataset')}
              style={{
                padding: '0.5rem 1rem',
                backgroundColor: '#f59e0b',
                color: 'white',
                border: 'none',
                borderRadius: '6px',
                cursor: 'pointer',
                fontSize: '0.875rem',
                fontWeight: 600,
              }}
            >
              📊 Create Dataset
            </button>
            <button
              onClick={() => console.log('Batch delete')}
              style={{
                padding: '0.5rem 1rem',
                backgroundColor: '#dc2626',
                color: 'white',
                border: 'none',
                borderRadius: '6px',
                cursor: 'pointer',
                fontSize: '0.875rem',
                fontWeight: 600,
              }}
            >
              🗑️ Archive
            </button>
          </div>
        )}
      </div>

      {/* File List */}
      <div style={{ backgroundColor: '#1f2937', borderRadius: '8px', border: '1px solid #374151', overflow: 'hidden' }}>
        {/* Table Header */}
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: '40px 1fr 120px 120px 180px 200px 120px',
            padding: '1rem',
            backgroundColor: '#111827',
            borderBottom: '1px solid #374151',
            fontWeight: 600,
            fontSize: '0.875rem',
            color: '#f3f4f6',
          }}
        >
          <div>
            <input
              type="checkbox"
              checked={selectedFiles.size === filteredFiles.length && filteredFiles.length > 0}
              onChange={handleSelectAll}
              style={{ cursor: 'pointer' }}
            />
          </div>
          <div>Name</div>
          <div>Type</div>
          <div>Size</div>
          <div>Created</div>
          <div>Source/Tags</div>
          <div>Actions</div>
        </div>

        {/* Table Body */}
        <div style={{ maxHeight: 'calc(100vh - 500px)', overflowY: 'auto' }}>
          {filteredFiles.length === 0 ? (
            <div style={{ padding: '3rem', textAlign: 'center', color: '#9ca3af' }}>
              <p>No files found</p>
            </div>
          ) : (
            filteredFiles.map((file, idx) => (
              <div
                key={file.id}
                style={{
                  display: 'grid',
                  gridTemplateColumns: '40px 1fr 120px 120px 180px 200px 120px',
                  padding: '1rem',
                  backgroundColor: idx % 2 === 0 ? '#1f2937' : '#111827',
                  borderBottom: '1px solid #374151',
                  fontSize: '0.875rem',
                  color: '#e5e7eb',
                  alignItems: 'center',
                }}
              >
                <div>
                  <input
                    type="checkbox"
                    checked={selectedFiles.has(file.id)}
                    onChange={() => handleSelectFile(file.id)}
                    style={{ cursor: 'pointer' }}
                  />
                </div>
                <div
                  style={{ fontWeight: 500, color: '#f3f4f6', cursor: 'pointer' }}
                  onClick={() => handlePreview(file)}
                >
                  {file.name}
                </div>
                <div>
                  <span
                    style={{
                      padding: '0.25rem 0.5rem',
                      backgroundColor: '#374151',
                      borderRadius: '4px',
                      fontSize: '0.75rem',
                    }}
                  >
                    {file.type}
                  </span>
                </div>
                <div style={{ color: '#9ca3af' }}>{file.size}</div>
                <div style={{ color: '#9ca3af', fontSize: '0.75rem' }}>{file.created_at}</div>
                <div>
                  {file.source && (
                    <span style={{ fontSize: '0.75rem', color: '#8b5cf6' }}>{file.source}</span>
                  )}
                  {file.tags && (
                    <div style={{ display: 'flex', gap: '0.25rem', flexWrap: 'wrap', marginTop: '0.25rem' }}>
                      {file.tags.map((tag) => (
                        <span
                          key={tag}
                          style={{
                            padding: '0.125rem 0.375rem',
                            backgroundColor: '#1e3a8a',
                            color: '#dbeafe',
                            borderRadius: '3px',
                            fontSize: '0.625rem',
                          }}
                        >
                          {tag}
                        </span>
                      ))}
                    </div>
                  )}
                </div>
                <div style={{ display: 'flex', gap: '0.5rem' }}>
                  <button
                    onClick={() => handlePreview(file)}
                    style={{
                      padding: '0.25rem 0.5rem',
                      backgroundColor: '#3b82f6',
                      color: 'white',
                      border: 'none',
                      borderRadius: '4px',
                      cursor: 'pointer',
                      fontSize: '0.75rem',
                    }}
                    title="Preview"
                  >
                    👁️
                  </button>
                  {(file.type === 'CSV' || file.type === 'TSV' || file.type === 'XLSX') && (
                    <button
                      onClick={() => handleVisualize(file)}
                      style={{
                        padding: '0.25rem 0.5rem',
                        backgroundColor: '#8b5cf6',
                        color: 'white',
                        border: 'none',
                        borderRadius: '4px',
                        cursor: 'pointer',
                        fontSize: '0.75rem',
                      }}
                      title="Quick Visualization"
                    >
                      📊
                    </button>
                  )}
                  <button
                    onClick={() => handleDownload(file.id)}
                    style={{
                      padding: '0.25rem 0.5rem',
                      backgroundColor: '#10b981',
                      color: 'white',
                      border: 'none',
                      borderRadius: '4px',
                      cursor: 'pointer',
                      fontSize: '0.75rem',
                    }}
                    title="Download"
                  >
                    💾
                  </button>
                  <button
                    onClick={() => handleDelete(file.id)}
                    style={{
                      padding: '0.25rem 0.5rem',
                      backgroundColor: '#dc2626',
                      color: 'white',
                      border: 'none',
                      borderRadius: '4px',
                      cursor: 'pointer',
                      fontSize: '0.75rem',
                    }}
                    title="Archive"
                  >
                    🗑️
                  </button>
                </div>
              </div>
            ))
          )}
        </div>
      </div>

      {/* File Preview Modal */}
      {showPreview && previewFile && (
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
          onClick={() => setShowPreview(false)}
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
                📄 {previewFile.name}
              </h2>
              <button
                onClick={() => setShowPreview(false)}
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

            <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '8px' }}>
              <p style={{ color: '#9ca3af', textAlign: 'center' }}>
                File preview feature coming soon...
              </p>
              <p style={{ color: '#9ca3af', textAlign: 'center', fontSize: '0.875rem' }}>
                Type: {previewFile.type} • Size: {previewFile.size}
              </p>
            </div>
          </div>
        </div>
      )}

      {/* Quick Visualizer Modal */}
      {visualizerFile && (
        <QuickVisualizerModal
          open={showVisualizer}
          onClose={() => {
            setShowVisualizer(false);
            setVisualizerFile(null);
          }}
          filePath={`/path/to/${visualizerFile.name}`} // TODO: Use actual file path
          fileName={visualizerFile.name}
        />
      )}
    </div>
  );
};

export default DataBrowserPage;
