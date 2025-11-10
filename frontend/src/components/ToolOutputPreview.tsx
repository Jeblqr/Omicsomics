import { useState } from 'react';

interface OutputFile {
  name: string;
  type: 'text' | 'image' | 'table' | 'json';
  content: any;
  size?: string;
  timestamp?: string;
}

interface ToolOutputPreviewProps {
  isOpen: boolean;
  onClose: () => void;
  toolName: string;
  outputFiles: OutputFile[];
}

const ToolOutputPreview = ({ isOpen, onClose, toolName, outputFiles }: ToolOutputPreviewProps) => {
  const [selectedFileIndex, setSelectedFileIndex] = useState(0);

  if (!isOpen || outputFiles.length === 0) return null;

  const selectedFile = outputFiles[selectedFileIndex];

  const renderFileContent = () => {
    if (!selectedFile) return null;

    switch (selectedFile.type) {
      case 'text':
        return (
          <pre
            style={{
              backgroundColor: '#111827',
              padding: '1rem',
              borderRadius: '6px',
              overflow: 'auto',
              maxHeight: '500px',
              color: '#e5e7eb',
              fontSize: '0.875rem',
              fontFamily: 'monospace',
              lineHeight: 1.5,
            }}
          >
            {selectedFile.content}
          </pre>
        );

      case 'image':
        return (
          <div style={{ textAlign: 'center', padding: '1rem', backgroundColor: '#111827', borderRadius: '6px' }}>
            <img
              src={selectedFile.content}
              alt={selectedFile.name}
              style={{ maxWidth: '100%', maxHeight: '500px', borderRadius: '4px' }}
            />
          </div>
        );

      case 'table':
        return (
          <div style={{ overflow: 'auto', maxHeight: '500px', backgroundColor: '#111827', borderRadius: '6px' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.875rem' }}>
              <thead style={{ position: 'sticky', top: 0, backgroundColor: '#1f2937', zIndex: 1 }}>
                <tr>
                  {selectedFile.content.headers.map((header: string, idx: number) => (
                    <th
                      key={idx}
                      style={{
                        padding: '0.75rem',
                        textAlign: 'left',
                        color: '#f3f4f6',
                        fontWeight: 600,
                        borderBottom: '2px solid #374151',
                      }}
                    >
                      {header}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {selectedFile.content.rows.map((row: any[], rowIdx: number) => (
                  <tr
                    key={rowIdx}
                    style={{
                      backgroundColor: rowIdx % 2 === 0 ? '#111827' : '#1f2937',
                      borderBottom: '1px solid #374151',
                    }}
                  >
                    {row.map((cell, cellIdx) => (
                      <td
                        key={cellIdx}
                        style={{
                          padding: '0.75rem',
                          color: '#e5e7eb',
                        }}
                      >
                        {typeof cell === 'number' ? cell.toFixed(4) : String(cell)}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        );

      case 'json':
        return (
          <pre
            style={{
              backgroundColor: '#111827',
              padding: '1rem',
              borderRadius: '6px',
              overflow: 'auto',
              maxHeight: '500px',
              color: '#e5e7eb',
              fontSize: '0.875rem',
              fontFamily: 'monospace',
              lineHeight: 1.5,
            }}
          >
            {JSON.stringify(selectedFile.content, null, 2)}
          </pre>
        );

      default:
        return (
          <div style={{ padding: '2rem', textAlign: 'center', color: '#9ca3af' }}>
            <p>Unsupported file type: {selectedFile.type}</p>
          </div>
        );
    }
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
          maxWidth: '1400px',
          width: '90%',
          maxHeight: '90vh',
          overflow: 'auto',
          border: '2px solid #374151',
        }}
        onClick={(e) => e.stopPropagation()}
      >
        {/* Header */}
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1.5rem' }}>
          <h2 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.5rem', fontWeight: 700 }}>
            📤 {toolName} - Output Preview
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

        {/* File Selector */}
        {outputFiles.length > 1 && (
          <div style={{ marginBottom: '1.5rem' }}>
            <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem', fontWeight: 600 }}>
              Output Files ({outputFiles.length})
            </div>
            <div style={{ display: 'flex', gap: '0.5rem', flexWrap: 'wrap' }}>
              {outputFiles.map((file, idx) => (
                <button
                  key={idx}
                  onClick={() => setSelectedFileIndex(idx)}
                  style={{
                    padding: '0.5rem 1rem',
                    backgroundColor: idx === selectedFileIndex ? '#3b82f6' : '#374151',
                    color: '#f3f4f6',
                    border: 'none',
                    borderRadius: '6px',
                    cursor: 'pointer',
                    fontSize: '0.875rem',
                    fontWeight: idx === selectedFileIndex ? 600 : 400,
                  }}
                >
                  {file.name}
                  {file.size && <span style={{ marginLeft: '0.5rem', color: '#9ca3af' }}>({file.size})</span>}
                </button>
              ))}
            </div>
          </div>
        )}

        {/* File Info */}
        {selectedFile && (
          <div
            style={{
              marginBottom: '1rem',
              padding: '1rem',
              backgroundColor: '#111827',
              borderRadius: '8px',
              border: '1px solid #374151',
            }}
          >
            <div style={{ display: 'flex', gap: '2rem', flexWrap: 'wrap', fontSize: '0.875rem' }}>
              <div>
                <span style={{ color: '#9ca3af' }}>File: </span>
                <span style={{ color: '#f3f4f6', fontWeight: 600 }}>{selectedFile.name}</span>
              </div>
              <div>
                <span style={{ color: '#9ca3af' }}>Type: </span>
                <span style={{ color: '#10b981', fontWeight: 600 }}>{selectedFile.type.toUpperCase()}</span>
              </div>
              {selectedFile.size && (
                <div>
                  <span style={{ color: '#9ca3af' }}>Size: </span>
                  <span style={{ color: '#f3f4f6', fontWeight: 600 }}>{selectedFile.size}</span>
                </div>
              )}
              {selectedFile.timestamp && (
                <div>
                  <span style={{ color: '#9ca3af' }}>Modified: </span>
                  <span style={{ color: '#f3f4f6', fontWeight: 600 }}>{selectedFile.timestamp}</span>
                </div>
              )}
            </div>
          </div>
        )}

        {/* File Content */}
        <div>{renderFileContent()}</div>

        {/* Actions */}
        <div style={{ marginTop: '1.5rem', display: 'flex', gap: '1rem', justifyContent: 'flex-end' }}>
          <button
            onClick={() => {
              // TODO: 实现下载功能
              console.log('Download:', selectedFile.name);
            }}
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
            💾 Download
          </button>
          <button
            onClick={() => {
              // TODO: 实现导出功能
              console.log('Export:', selectedFile.name);
            }}
            style={{
              padding: '0.75rem 1.5rem',
              backgroundColor: '#8b5cf6',
              color: 'white',
              border: 'none',
              borderRadius: '6px',
              cursor: 'pointer',
              fontWeight: 600,
            }}
          >
            📊 Export to Data
          </button>
        </div>
      </div>
    </div>
  );
};

export default ToolOutputPreview;
