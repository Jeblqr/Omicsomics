import { useState, useEffect } from 'react';
import ReactFlow, { Node, Edge, MarkerType, BackgroundVariant, Background, Controls } from 'reactflow';
import 'reactflow/dist/style.css';

interface ToolExecution {
  id: string;
  name: string;
  status: 'waiting' | 'running' | 'completed' | 'failed';
  startTime?: string;
  endTime?: string;
  duration?: number;
  progress?: number;
  cpuUsage?: number;
  memoryUsage?: number;
  error?: string;
}

interface RunExecutionVisualizerProps {
  runId: string;
  pipelineName: string;
  nodes: any[];
  edges: any[];
  toolExecutions: ToolExecution[];
  onNodeClick?: (nodeId: string) => void;
}

const RunExecutionVisualizer = ({
  runId,
  pipelineName,
  nodes: pipelineNodes,
  edges: pipelineEdges,
  toolExecutions,
  onNodeClick,
}: RunExecutionVisualizerProps) => {
  const [flowNodes, setFlowNodes] = useState<Node[]>([]);
  const [flowEdges, setFlowEdges] = useState<Edge[]>([]);

  useEffect(() => {
    // 转换节点，根据执行状态添加样式
    const convertedNodes: Node[] = pipelineNodes.map((node) => {
      const execution = toolExecutions.find((e) => e.id === node.id);
      const status = execution?.status || 'waiting';

      // 根据状态选择颜色
      const statusColors = {
        waiting: { bg: '#374151', border: '#4b5563', text: '#9ca3af' },
        running: { bg: '#1e3a8a', border: '#3b82f6', text: '#dbeafe' },
        completed: { bg: '#064e3b', border: '#10b981', text: '#d1fae5' },
        failed: { bg: '#7f1d1d', border: '#dc2626', text: '#fecaca' },
      };

      const colors = statusColors[status];

      return {
        id: node.id,
        type: 'default',
        position: node.position,
        data: {
          label: (
            <div style={{ textAlign: 'center' }}>
              <div style={{ fontWeight: 600, marginBottom: '0.25rem' }}>{node.data.label}</div>
              <div style={{ fontSize: '0.75rem', opacity: 0.8 }}>
                {status === 'running' && execution?.progress
                  ? `${execution.progress}%`
                  : status.charAt(0).toUpperCase() + status.slice(1)}
              </div>
            </div>
          ),
        },
        style: {
          backgroundColor: colors.bg,
          border: `2px solid ${colors.border}`,
          color: colors.text,
          padding: '12px',
          borderRadius: '8px',
          minWidth: '120px',
          boxShadow: status === 'running' ? `0 0 15px ${colors.border}` : 'none',
        },
      };
    });

    // 转换边，高亮已执行的路径
    const convertedEdges: Edge[] = pipelineEdges.map((edge) => {
      const sourceExecution = toolExecutions.find((e) => e.id === edge.source);
      const targetExecution = toolExecutions.find((e) => e.id === edge.target);

      const isActive =
        sourceExecution?.status === 'completed' ||
        sourceExecution?.status === 'running' ||
        targetExecution?.status === 'running';

      return {
        id: edge.id,
        source: edge.source,
        target: edge.target,
        type: 'smoothstep',
        animated: targetExecution?.status === 'running',
        style: {
          stroke: isActive ? '#3b82f6' : '#4b5563',
          strokeWidth: isActive ? 3 : 2,
        },
        markerEnd: {
          type: MarkerType.ArrowClosed,
          color: isActive ? '#3b82f6' : '#4b5563',
        },
      };
    });

    setFlowNodes(convertedNodes);
    setFlowEdges(convertedEdges);
  }, [pipelineNodes, pipelineEdges, toolExecutions]);

  // 计算总体统计
  const stats = {
    total: toolExecutions.length,
    waiting: toolExecutions.filter((e) => e.status === 'waiting').length,
    running: toolExecutions.filter((e) => e.status === 'running').length,
    completed: toolExecutions.filter((e) => e.status === 'completed').length,
    failed: toolExecutions.filter((e) => e.status === 'failed').length,
    progress: toolExecutions.length > 0 ? (toolExecutions.filter((e) => e.status === 'completed').length / toolExecutions.length) * 100 : 0,
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100%', backgroundColor: '#0f172a' }}>
      {/* Header with Stats */}
      <div style={{ padding: '1.5rem', backgroundColor: '#1f2937', borderBottom: '1px solid #374151' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
          <div>
            <h2 style={{ margin: '0 0 0.5rem 0', color: '#f3f4f6', fontSize: '1.5rem', fontWeight: 700 }}>
              ⚡ {pipelineName}
            </h2>
            <div style={{ fontSize: '0.875rem', color: '#9ca3af' }}>Run ID: {runId}</div>
          </div>
          <div
            style={{
              padding: '0.75rem 1.5rem',
              backgroundColor: '#111827',
              borderRadius: '8px',
              border: '1px solid #374151',
            }}
          >
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem', textTransform: 'uppercase' }}>
              Overall Progress
            </div>
            <div style={{ fontSize: '1.5rem', color: '#3b82f6', fontWeight: 700 }}>{stats.progress.toFixed(0)}%</div>
          </div>
        </div>

        {/* Progress Bar */}
        <div
          style={{
            height: '8px',
            backgroundColor: '#374151',
            borderRadius: '4px',
            overflow: 'hidden',
            marginBottom: '1rem',
          }}
        >
          <div
            style={{
              height: '100%',
              backgroundColor: '#3b82f6',
              width: `${stats.progress}%`,
              transition: 'width 0.3s ease',
            }}
          />
        </div>

        {/* Status Stats */}
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(120px, 1fr))', gap: '1rem' }}>
          <div style={{ padding: '0.75rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>Total</div>
            <div style={{ fontSize: '1.25rem', color: '#f3f4f6', fontWeight: 600 }}>{stats.total}</div>
          </div>
          <div style={{ padding: '0.75rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>⏳ Waiting</div>
            <div style={{ fontSize: '1.25rem', color: '#9ca3af', fontWeight: 600 }}>{stats.waiting}</div>
          </div>
          <div style={{ padding: '0.75rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>▶️ Running</div>
            <div style={{ fontSize: '1.25rem', color: '#3b82f6', fontWeight: 600 }}>{stats.running}</div>
          </div>
          <div style={{ padding: '0.75rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>✅ Completed</div>
            <div style={{ fontSize: '1.25rem', color: '#10b981', fontWeight: 600 }}>{stats.completed}</div>
          </div>
          <div style={{ padding: '0.75rem', backgroundColor: '#111827', borderRadius: '6px', border: '1px solid #374151' }}>
            <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginBottom: '0.25rem' }}>❌ Failed</div>
            <div style={{ fontSize: '1.25rem', color: '#dc2626', fontWeight: 600 }}>{stats.failed}</div>
          </div>
        </div>
      </div>

      {/* Pipeline Visualization */}
      <div style={{ flex: 1, position: 'relative' }}>
        <ReactFlow
          nodes={flowNodes}
          edges={flowEdges}
          onNodeClick={(_, node) => onNodeClick && onNodeClick(node.id)}
          fitView
          attributionPosition="bottom-left"
        >
          <Background variant={BackgroundVariant.Dots} gap={12} size={1} color="#374151" />
          <Controls />
        </ReactFlow>
      </div>

      {/* Tool Execution Details */}
      <div
        style={{
          maxHeight: '250px',
          overflowY: 'auto',
          backgroundColor: '#1f2937',
          borderTop: '1px solid #374151',
          padding: '1rem',
        }}
      >
        <h3 style={{ margin: '0 0 1rem 0', color: '#f3f4f6', fontSize: '1rem', fontWeight: 600 }}>
          🔧 Tool Execution Details
        </h3>
        <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
          {toolExecutions.map((execution) => (
            <div
              key={execution.id}
              style={{
                padding: '1rem',
                backgroundColor: '#111827',
                borderRadius: '6px',
                border: '1px solid #374151',
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
              }}
            >
              <div style={{ flex: 1 }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '0.75rem', marginBottom: '0.5rem' }}>
                  <div style={{ fontWeight: 600, color: '#f3f4f6' }}>{execution.name}</div>
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
                  {execution.status === 'running' && execution.progress !== undefined && (
                    <div style={{ fontSize: '0.875rem', color: '#3b82f6', fontWeight: 600 }}>
                      {execution.progress}%
                    </div>
                  )}
                </div>
                <div style={{ display: 'flex', gap: '1.5rem', fontSize: '0.75rem', color: '#9ca3af' }}>
                  {execution.startTime && (
                    <div>
                      Start: <span style={{ color: '#f3f4f6' }}>{execution.startTime}</span>
                    </div>
                  )}
                  {execution.duration && (
                    <div>
                      Duration: <span style={{ color: '#f3f4f6' }}>{execution.duration}s</span>
                    </div>
                  )}
                  {execution.cpuUsage !== undefined && (
                    <div>
                      CPU: <span style={{ color: '#f3f4f6' }}>{execution.cpuUsage}%</span>
                    </div>
                  )}
                  {execution.memoryUsage !== undefined && (
                    <div>
                      Memory: <span style={{ color: '#f3f4f6' }}>{execution.memoryUsage} MB</span>
                    </div>
                  )}
                </div>
                {execution.error && (
                  <div
                    style={{
                      marginTop: '0.5rem',
                      padding: '0.5rem',
                      backgroundColor: '#7f1d1d',
                      borderRadius: '4px',
                      fontSize: '0.75rem',
                      color: '#fecaca',
                    }}
                  >
                    ⚠️ {execution.error}
                  </div>
                )}
              </div>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
};

export default RunExecutionVisualizer;
