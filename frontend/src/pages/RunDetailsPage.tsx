import { useState } from 'react';
import RunExecutionVisualizer from '../components/RunExecutionVisualizer';
import ToolOutputPreview from '../components/ToolOutputPreview';

interface RunFile {
  id: string;
  name: string;
  type: string;
  size: string;
  path: string;
}

interface ToolExecution {
  id: string;
  name: string;
  tool: string;
  status: 'waiting' | 'running' | 'completed' | 'failed';
  startTime?: string;
  endTime?: string;
  duration?: number;
  progress?: number;
  cpuUsage?: number;
  memoryUsage?: number;
  inputFiles: RunFile[];
  outputFiles: RunFile[];
  logs: string[];
  error?: string;
}

interface RunDetailsPageProps {
  runId: string;
  pipelineName: string;
  status: 'running' | 'completed' | 'failed';
  startTime: string;
  endTime?: string;
  nodes: any[];
  edges: any[];
  toolExecutions: ToolExecution[];
}

const RunDetailsPage = ({
  runId,
  pipelineName,
  status,
  startTime,
  endTime,
  nodes,
  edges,
  toolExecutions,
}: RunDetailsPageProps) => {
  const [selectedView, setSelectedView] = useState<'overview' | 'timeline' | 'logs' | 'files'>('overview');
  const [selectedToolId, setSelectedToolId] = useState<string | null>(null);
  const [showOutputPreview, setShowOutputPreview] = useState(false);

  const selectedTool = toolExecutions.find((t) => t.id === selectedToolId);

  // 计算总体指标
  const totalDuration = toolExecutions.reduce((sum, t) => sum + (t.duration || 0), 0);
  const totalCpuUsage = toolExecutions.reduce((sum, t) => sum + (t.cpuUsage || 0), 0) / toolExecutions.length;
  const totalMemoryUsage = toolExecutions.reduce((sum, t) => sum + (t.memoryUsage || 0), 0);

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100vh', backgroundColor: '#0f172a' }}>
      {/* Header */}
      <div style={{ padding: '1.5rem', backgroundColor: '#1f2937', borderBottom: '1px solid #374151' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
          <div>
            <h1 style={{ margin: '0 0 0.5rem 0', color: '#f3f4f6', fontSize: '1.75rem', fontWeight: 700 }}>
              Run Details
            </h1>
            <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>
              {pipelineName} • Run ID: {runId}
            </div>
          </div>
          <div
            style={{
              padding: '0.75rem 1.5rem',
              borderRadius: '6px',
              fontWeight: 600,
              backgroundColor:
                status === 'completed' ? '#064e3b' : status === 'running' ? '#1e3a8a' : '#7f1d1d',
              color: status === 'completed' ? '#d1fae5' : status === 'running' ? '#dbeafe' : '#fecaca',
            }}
          >
            {status.toUpperCase()}
          </div>
        </div>

        {/* Summary Stats */}
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '1rem' }}>
          <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>Start Time</div>
            <div style={{ fontSize: '0.875rem', color: '#f3f4f6', fontWeight: 600 }}>{startTime}</div>
          </div>
          {endTime && (
            <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
              <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>End Time</div>
              <div style={{ fontSize: '0.875rem', color: '#f3f4f6', fontWeight: 600 }}>{endTime}</div>
            </div>
          )}
          <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>Total Duration</div>
            <div style={{ fontSize: '0.875rem', color: '#f3f4f6', fontWeight: 600 }}>{totalDuration.toFixed(1)}s</div>
          </div>
          <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>Avg CPU</div>
            <div style={{ fontSize: '0.875rem', color: '#f3f4f6', fontWeight: 600 }}>{totalCpuUsage.toFixed(1)}%</div>
          </div>
          <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>Total Memory</div>
            <div style={{ fontSize: '0.875rem', color: '#f3f4f6', fontWeight: 600 }}>{totalMemoryUsage.toFixed(0)} MB</div>
          </div>
        </div>

        {/* View Tabs */}
        <div style={{ display: 'flex', gap: '0.5rem', marginTop: '1rem' }}>
          {(['overview', 'timeline', 'logs', 'files'] as const).map((view) => (
            <button
              key={view}
              onClick={() => setSelectedView(view)}
              style={{
                padding: '0.75rem 1.5rem',
                backgroundColor: selectedView === view ? '#3b82f6' : '#374151',
                color: '#f3f4f6',
                border: 'none',
                borderRadius: '6px',
                cursor: 'pointer',
                fontWeight: 600,
                fontSize: '0.875rem',
              }}
            >
              {view.charAt(0).toUpperCase() + view.slice(1)}
            </button>
          ))}
        </div>
      </div>

      {/* Content */}
      <div style={{ flex: 1, overflow: 'auto' }}>
        {selectedView === 'overview' && (
          <div style={{ height: '100%' }}>
            <RunExecutionVisualizer
              runId={runId}
              pipelineName={pipelineName}
              nodes={nodes}
              edges={edges}
              toolExecutions={toolExecutions}
              onNodeClick={(nodeId) => setSelectedToolId(nodeId)}
            />
          </div>
        )}

        {selectedView === 'timeline' && (
          <div style={{ padding: '2rem' }}>
            <h2 style={{ color: '#f3f4f6', marginBottom: '1.5rem' }}>⏱️ Execution Timeline</h2>
            <div style={{ position: 'relative', paddingLeft: '2rem' }}>
              {/* Vertical line */}
              <div
                style={{
                  position: 'absolute',
                  left: '0.75rem',
                  top: 0,
                  bottom: 0,
                  width: '2px',
                  backgroundColor: '#374151',
                }}
              />

              {/* Timeline items */}
              {toolExecutions.map((execution, idx) => (
                <div key={execution.id} style={{ position: 'relative', marginBottom: '2rem' }}>
                  {/* Dot */}
                  <div
                    style={{
                      position: 'absolute',
                      left: '-1.375rem',
                      top: '0.5rem',
                      width: '1rem',
                      height: '1rem',
                      borderRadius: '50%',
                      backgroundColor:
                        execution.status === 'completed'
                          ? '#10b981'
                          : execution.status === 'running'
                          ? '#3b82f6'
                          : execution.status === 'failed'
                          ? '#dc2626'
                          : '#6b7280',
                      border: '3px solid #1f2937',
                    }}
                  />

                  {/* Content */}
                  <div
                    style={{
                      padding: '1.5rem',
                      backgroundColor: '#1f2937',
                      borderRadius: '8px',
                      border: '1px solid #374151',
                      cursor: 'pointer',
                    }}
                    onClick={() => setSelectedToolId(execution.id)}
                  >
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'start', marginBottom: '0.75rem' }}>
                      <div>
                        <h3 style={{ margin: '0 0 0.25rem 0', color: '#f3f4f6', fontSize: '1.125rem', fontWeight: 600 }}>
                          {execution.name}
                        </h3>
                        <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>{execution.tool}</div>
                      </div>
                      <div
                        style={{
                          padding: '0.25rem 0.75rem',
                          borderRadius: '4px',
                          fontSize: '0.75rem',
                          fontWeight: 600,
                          backgroundColor:
                            execution.status === 'completed'
                              ? '#064e3b'
                              : execution.status === 'running'
                              ? '#1e3a8a'
                              : execution.status === 'failed'
                              ? '#7f1d1d'
                              : '#374151',
                          color:
                            execution.status === 'completed'
                              ? '#d1fae5'
                              : execution.status === 'running'
                              ? '#dbeafe'
                              : execution.status === 'failed'
                              ? '#fecaca'
                              : '#9ca3af',
                        }}
                      >
                        {execution.status.toUpperCase()}
                      </div>
                    </div>

                    <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '1rem', fontSize: '0.875rem' }}>
                      {execution.startTime && (
                        <div>
                          <span style={{ color: '#9ca3af' }}>Start: </span>
                          <span style={{ color: '#f3f4f6' }}>{execution.startTime}</span>
                        </div>
                      )}
                      {execution.duration && (
                        <div>
                          <span style={{ color: '#9ca3af' }}>Duration: </span>
                          <span style={{ color: '#f3f4f6' }}>{execution.duration}s</span>
                        </div>
                      )}
                      <div>
                        <span style={{ color: '#9ca3af' }}>Input Files: </span>
                        <span style={{ color: '#f3f4f6' }}>{execution.inputFiles.length}</span>
                      </div>
                      <div>
                        <span style={{ color: '#9ca3af' }}>Output Files: </span>
                        <span style={{ color: '#f3f4f6' }}>{execution.outputFiles.length}</span>
                      </div>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {selectedView === 'logs' && (
          <div style={{ padding: '2rem' }}>
            <h2 style={{ color: '#f3f4f6', marginBottom: '1.5rem' }}>📋 Execution Logs</h2>
            {selectedToolId ? (
              <div>
                <div style={{ marginBottom: '1rem', fontSize: '0.875rem', color: '#9ca3af' }}>
                  Showing logs for: <span style={{ color: '#f3f4f6', fontWeight: 600 }}>{selectedTool?.name}</span>
                  <button
                    onClick={() => setSelectedToolId(null)}
                    style={{
                      marginLeft: '1rem',
                      padding: '0.25rem 0.75rem',
                      backgroundColor: '#374151',
                      color: '#f3f4f6',
                      border: 'none',
                      borderRadius: '4px',
                      cursor: 'pointer',
                      fontSize: '0.75rem',
                    }}
                  >
                    Show All
                  </button>
                </div>
                <pre
                  style={{
                    backgroundColor: '#111827',
                    padding: '1rem',
                    borderRadius: '6px',
                    overflow: 'auto',
                    maxHeight: '600px',
                    color: '#e5e7eb',
                    fontSize: '0.875rem',
                    fontFamily: 'monospace',
                    lineHeight: 1.6,
                  }}
                >
                  {selectedTool?.logs.join('\n') || 'No logs available'}
                </pre>
              </div>
            ) : (
              <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                {toolExecutions.map((execution) => (
                  <div
                    key={execution.id}
                    style={{
                      padding: '1rem',
                      backgroundColor: '#1f2937',
                      borderRadius: '8px',
                      border: '1px solid #374151',
                      cursor: 'pointer',
                    }}
                    onClick={() => setSelectedToolId(execution.id)}
                  >
                    <div style={{ fontWeight: 600, color: '#f3f4f6', marginBottom: '0.5rem' }}>{execution.name}</div>
                    <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>
                      {execution.logs.length} log entries
                    </div>
                  </div>
                ))}
              </div>
            )}
          </div>
        )}

        {selectedView === 'files' && (
          <div style={{ padding: '2rem' }}>
            <h2 style={{ color: '#f3f4f6', marginBottom: '1.5rem' }}>📁 Input/Output Files</h2>
            <div style={{ display: 'flex', flexDirection: 'column', gap: '2rem' }}>
              {toolExecutions.map((execution) => (
                <div
                  key={execution.id}
                  style={{
                    padding: '1.5rem',
                    backgroundColor: '#1f2937',
                    borderRadius: '8px',
                    border: '1px solid #374151',
                  }}
                >
                  <h3 style={{ margin: '0 0 1rem 0', color: '#f3f4f6', fontSize: '1.125rem', fontWeight: 600 }}>
                    {execution.name}
                  </h3>

                  <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1.5rem' }}>
                    {/* Input Files */}
                    <div>
                      <h4 style={{ margin: '0 0 0.75rem 0', color: '#9ca3af', fontSize: '0.875rem', fontWeight: 600 }}>
                        📥 Input Files ({execution.inputFiles.length})
                      </h4>
                      <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
                        {execution.inputFiles.map((file) => (
                          <div
                            key={file.id}
                            style={{
                              padding: '0.75rem',
                              backgroundColor: '#111827',
                              borderRadius: '6px',
                              border: '1px solid #374151',
                            }}
                          >
                            <div style={{ fontSize: '0.875rem', color: '#f3f4f6', fontWeight: 500 }}>{file.name}</div>
                            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginTop: '0.25rem' }}>
                              {file.type} • {file.size}
                            </div>
                          </div>
                        ))}
                      </div>
                    </div>

                    {/* Output Files */}
                    <div>
                      <h4 style={{ margin: '0 0 0.75rem 0', color: '#9ca3af', fontSize: '0.875rem', fontWeight: 600 }}>
                        📤 Output Files ({execution.outputFiles.length})
                      </h4>
                      <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
                        {execution.outputFiles.map((file) => (
                          <div
                            key={file.id}
                            style={{
                              padding: '0.75rem',
                              backgroundColor: '#111827',
                              borderRadius: '6px',
                              border: '1px solid #374151',
                              cursor: 'pointer',
                            }}
                            onClick={() => {
                              setSelectedToolId(execution.id);
                              setShowOutputPreview(true);
                            }}
                          >
                            <div style={{ fontSize: '0.875rem', color: '#f3f4f6', fontWeight: 500 }}>{file.name}</div>
                            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginTop: '0.25rem' }}>
                              {file.type} • {file.size}
                            </div>
                          </div>
                        ))}
                      </div>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}
      </div>

      {/* Tool Output Preview Modal */}
      {showOutputPreview && selectedTool && (
        <ToolOutputPreview
          isOpen={showOutputPreview}
          onClose={() => setShowOutputPreview(false)}
          toolName={selectedTool.name}
          outputFiles={selectedTool.outputFiles.map((f) => ({
            name: f.name,
            type: 'text', // TODO: 根据实际文件类型判断
            content: 'Sample output content...',
            size: f.size,
          }))}
        />
      )}
    </div>
  );
};

export default RunDetailsPage;
