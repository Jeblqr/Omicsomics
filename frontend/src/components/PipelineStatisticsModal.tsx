/**
 * Pipeline Statistics Modal
 * 管道统计信息模态框 - 显示整个管道的概览和统计
 */

import React from 'react';
import { Node, Edge } from 'reactflow';
import { PipelineNodeData } from './pipelines/PipelineNode';

interface PipelineStatisticsProps {
  isOpen: boolean;
  onClose: () => void;
  pipelineName: string;
  nodes: Node<PipelineNodeData>[];
  edges: Edge[];
}

const PipelineStatisticsModal: React.FC<PipelineStatisticsProps> = ({
  isOpen,
  onClose,
  pipelineName,
  nodes,
  edges,
}) => {
  if (!isOpen) return null;

  // 按类型统计节点
  const nodesByType = nodes.reduce((acc, node) => {
    const type = node.data.nodeType || 'unknown';
    acc[type] = (acc[type] || 0) + 1;
    return acc;
  }, {} as Record<string, number>);

  // 检测孤立节点（没有连接的节点）
  const isolatedNodes = nodes.filter((node) => {
    const hasIncoming = edges.some((edge) => edge.target === node.id);
    const hasOutgoing = edges.some((edge) => edge.source === node.id);
    return !hasIncoming && !hasOutgoing;
  });

  // 入口节点（没有输入连接）
  const inputNodes = nodes.filter((node) => {
    return !edges.some((edge) => edge.target === node.id);
  });

  // 出口节点（没有输出连接）
  const outputNodes = nodes.filter((node) => {
    return !edges.some((edge) => edge.source === node.id);
  });

  // 计算管道深度（最长路径）
  const calculateMaxDepth = () => {
    const depths = new Map<string, number>();
    
    // 初始化入口节点深度为0
    inputNodes.forEach((node) => depths.set(node.id, 0));

    // DFS计算深度
    const visited = new Set<string>();
    const calculateDepth = (nodeId: string): number => {
      if (visited.has(nodeId)) return depths.get(nodeId) || 0;
      visited.add(nodeId);

      const incomingEdges = edges.filter((e) => e.target === nodeId);
      if (incomingEdges.length === 0) {
        depths.set(nodeId, 0);
        return 0;
      }

      const maxParentDepth = Math.max(
        ...incomingEdges.map((e) => calculateDepth(e.source))
      );
      const depth = maxParentDepth + 1;
      depths.set(nodeId, depth);
      return depth;
    };

    nodes.forEach((node) => calculateDepth(node.id));
    return Math.max(...Array.from(depths.values()), 0);
  };

  const maxDepth = calculateMaxDepth();

  // 按工具统计
  const toolCounts = nodes.reduce((acc, node) => {
    const tool = node.data.tool || 'Unknown';
    acc[tool] = (acc[tool] || 0) + 1;
    return acc;
  }, {} as Record<string, number>);

  // 类型图标映射
  const typeIcons: Record<string, string> = {
    input: '📥',
    output: '📤',
    processing: '⚙️',
    qc: '✅',
    alignment: '🧬',
    variant_calling: '🔬',
    visualization: '📊',
    analysis: '🔍',
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
        padding: '2rem',
      }}
      onClick={onClose}
    >
      <div
        style={{
          backgroundColor: '#1f2937',
          borderRadius: '12px',
          maxWidth: '900px',
          width: '100%',
          maxHeight: '90vh',
          overflow: 'hidden',
          display: 'flex',
          flexDirection: 'column',
          border: '2px solid #374151',
        }}
        onClick={(e) => e.stopPropagation()}
      >
        {/* Header */}
        <div
          style={{
            padding: '1.5rem',
            borderBottom: '2px solid #374151',
            backgroundColor: '#111827',
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
          }}
        >
          <h2 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.25rem', fontWeight: 600 }}>
            📈 Pipeline Statistics: {pipelineName}
          </h2>
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

        {/* Content */}
        <div style={{ flex: 1, overflowY: 'auto', padding: '1.5rem' }}>
          {/* Overall Statistics */}
          <div style={{ marginBottom: '2rem' }}>
            <h3 style={{ color: '#f3f4f6', fontSize: '1rem', fontWeight: 600, marginBottom: '1rem' }}>
              📊 Overall Statistics
            </h3>
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '1rem' }}>
              <div
                style={{
                  padding: '1.5rem',
                  backgroundColor: '#111827',
                  borderRadius: '8px',
                  border: '1px solid #374151',
                }}
              >
                <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem' }}>Total Nodes</div>
                <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#3b82f6' }}>{nodes.length}</div>
              </div>
              <div
                style={{
                  padding: '1.5rem',
                  backgroundColor: '#111827',
                  borderRadius: '8px',
                  border: '1px solid #374151',
                }}
              >
                <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem' }}>Total Connections</div>
                <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#10b981' }}>{edges.length}</div>
              </div>
              <div
                style={{
                  padding: '1.5rem',
                  backgroundColor: '#111827',
                  borderRadius: '8px',
                  border: '1px solid #374151',
                }}
              >
                <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.5rem' }}>Pipeline Depth</div>
                <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#8b5cf6' }}>{maxDepth + 1}</div>
              </div>
            </div>
          </div>

          {/* Node Types */}
          <div style={{ marginBottom: '2rem' }}>
            <h3 style={{ color: '#f3f4f6', fontSize: '1rem', fontWeight: 600, marginBottom: '1rem' }}>
              🏷️ Node Types
            </h3>
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '0.75rem' }}>
              {Object.entries(nodesByType).map(([type, count]) => (
                <div
                  key={type}
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
                  <span style={{ fontSize: '0.875rem', color: '#e5e7eb' }}>
                    {typeIcons[type] || '⚪'} {type}
                  </span>
                  <span
                    style={{
                      padding: '0.25rem 0.5rem',
                      backgroundColor: '#3b82f6',
                      color: 'white',
                      borderRadius: '4px',
                      fontSize: '0.75rem',
                      fontWeight: 600,
                    }}
                  >
                    {count}
                  </span>
                </div>
              ))}
            </div>
          </div>

          {/* Input/Output Nodes */}
          <div style={{ marginBottom: '2rem' }}>
            <h3 style={{ color: '#f3f4f6', fontSize: '1rem', fontWeight: 600, marginBottom: '1rem' }}>
              🔗 Entry & Exit Points
            </h3>
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(2, 1fr)', gap: '1rem' }}>
              <div
                style={{
                  padding: '1rem',
                  backgroundColor: '#111827',
                  borderRadius: '6px',
                  border: '1px solid #374151',
                }}
              >
                <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.75rem', fontWeight: 600 }}>
                  📥 Input Nodes ({inputNodes.length})
                </div>
                {inputNodes.length === 0 ? (
                  <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>No input nodes</div>
                ) : (
                  <ul style={{ margin: 0, padding: 0, listStyle: 'none' }}>
                    {inputNodes.map((node) => (
                      <li
                        key={node.id}
                        style={{
                          padding: '0.5rem',
                          fontSize: '0.75rem',
                          color: '#e5e7eb',
                          backgroundColor: '#1f2937',
                          borderRadius: '4px',
                          marginBottom: '0.25rem',
                        }}
                      >
                        {node.data.label}
                      </li>
                    ))}
                  </ul>
                )}
              </div>
              <div
                style={{
                  padding: '1rem',
                  backgroundColor: '#111827',
                  borderRadius: '6px',
                  border: '1px solid #374151',
                }}
              >
                <div style={{ fontSize: '0.875rem', color: '#9ca3af', marginBottom: '0.75rem', fontWeight: 600 }}>
                  📤 Output Nodes ({outputNodes.length})
                </div>
                {outputNodes.length === 0 ? (
                  <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>No output nodes</div>
                ) : (
                  <ul style={{ margin: 0, padding: 0, listStyle: 'none' }}>
                    {outputNodes.map((node) => (
                      <li
                        key={node.id}
                        style={{
                          padding: '0.5rem',
                          fontSize: '0.75rem',
                          color: '#e5e7eb',
                          backgroundColor: '#1f2937',
                          borderRadius: '4px',
                          marginBottom: '0.25rem',
                        }}
                      >
                        {node.data.label}
                      </li>
                    ))}
                  </ul>
                )}
              </div>
            </div>
          </div>

          {/* Tools Used */}
          <div style={{ marginBottom: '2rem' }}>
            <h3 style={{ color: '#f3f4f6', fontSize: '1rem', fontWeight: 600, marginBottom: '1rem' }}>
              🔧 Tools Used
            </h3>
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(200px, 1fr))', gap: '0.75rem' }}>
              {Object.entries(toolCounts)
                .sort((a, b) => b[1] - a[1])
                .map(([tool, count]) => (
                  <div
                    key={tool}
                    style={{
                      padding: '0.75rem',
                      backgroundColor: '#111827',
                      borderRadius: '6px',
                      border: '1px solid #374151',
                      display: 'flex',
                      justifyContent: 'space-between',
                      alignItems: 'center',
                    }}
                  >
                    <span style={{ fontSize: '0.875rem', color: '#e5e7eb' }}>{tool}</span>
                    <span
                      style={{
                        padding: '0.25rem 0.5rem',
                        backgroundColor: '#059669',
                        color: 'white',
                        borderRadius: '4px',
                        fontSize: '0.75rem',
                        fontWeight: 600,
                      }}
                    >
                      ×{count}
                    </span>
                  </div>
                ))}
            </div>
          </div>

          {/* Warnings */}
          {isolatedNodes.length > 0 && (
            <div
              style={{
                padding: '1rem',
                backgroundColor: '#7f1d1d',
                borderRadius: '6px',
                border: '1px solid #991b1b',
              }}
            >
              <div style={{ fontSize: '0.875rem', color: '#fecaca', fontWeight: 600, marginBottom: '0.5rem' }}>
                ⚠️ Warnings
              </div>
              <div style={{ fontSize: '0.75rem', color: '#fca5a5' }}>
                {isolatedNodes.length} isolated node(s) detected (not connected to the pipeline):
              </div>
              <ul style={{ margin: '0.5rem 0 0 0', padding: 0, listStyle: 'none' }}>
                {isolatedNodes.map((node) => (
                  <li key={node.id} style={{ padding: '0.25rem 0', fontSize: '0.75rem', color: '#fca5a5' }}>
                    • {node.data.label}
                  </li>
                ))}
              </ul>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default PipelineStatisticsModal;
