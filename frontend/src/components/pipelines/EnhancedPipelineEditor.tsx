import { useState, useCallback, useRef, DragEvent } from 'react';
import ReactFlow, {
  Node,
  Edge,
  Controls,
  Background,
  applyNodeChanges,
  applyEdgeChanges,
  addEdge,
  NodeChange,
  EdgeChange,
  Connection,
  MarkerType,
  BackgroundVariant,
  MiniMap,
} from 'reactflow';
import 'reactflow/dist/style.css';

import PipelineNode, { PipelineNodeData } from './PipelineNode';
import ToolboxV2 from './ToolboxV2';
import { ToolDefinition } from './ToolDefinitions';
import ConfigPanelV2 from './ConfigPanelV2';
import DataVisualizationModal from '../DataVisualizationModal';
import PipelineStatisticsModal from '../PipelineStatisticsModal';
import PipelinePreviewModal from '../PipelinePreviewModal';
import {
  validatePipeline,
  downloadPipelineJSON,
  importPipeline,
  getPipelineTemplates,
  type ValidationResult,
  type PipelineExport,
} from '../../utils/pipelineValidation';

interface EnhancedPipelineEditorProps {
  initialNodes?: Node<PipelineNodeData>[];
  initialEdges?: Edge[];
  onSave?: (nodes: Node<PipelineNodeData>[], edges: Edge[]) => void;
  readOnly?: boolean;
}

const nodeTypes = {
  pipelineNode: PipelineNode,
};

const EnhancedPipelineEditor = ({
  initialNodes = [],
  initialEdges = [],
  onSave,
  readOnly = false,
}: EnhancedPipelineEditorProps) => {
  const [nodes, setNodes] = useState<Node<PipelineNodeData>[]>(initialNodes);
  const [edges, setEdges] = useState<Edge[]>(initialEdges);
  const [selectedNode, setSelectedNode] = useState<Node<PipelineNodeData> | null>(null);
  const [showConfigPanel, setShowConfigPanel] = useState(false);
  const [validationResult, setValidationResult] = useState<ValidationResult | null>(null);
  const [showValidation, setShowValidation] = useState(false);
  const [showImportExport, setShowImportExport] = useState(false);
  const [showDataViz, setShowDataViz] = useState(false);
  const [showPipelineStats, setShowPipelineStats] = useState(false);
  const [showPreview, setShowPreview] = useState(false);
  const [vizData, setVizData] = useState<any[]>([]);
  const nodeIdCounter = useRef(initialNodes.length || 0);
  const reactFlowWrapper = useRef<HTMLDivElement>(null);
  const [reactFlowInstance, setReactFlowInstance] = useState<any>(null);

  const onNodesChange = useCallback(
    (changes: NodeChange[]) => {
      if (readOnly) return;
      setNodes((nds) => applyNodeChanges(changes, nds) as Node<PipelineNodeData>[]);
    },
    [readOnly]
  );

  const onEdgesChange = useCallback(
    (changes: EdgeChange[]) => {
      if (readOnly) return;
      setEdges((eds) => applyEdgeChanges(changes, eds));
    },
    [readOnly]
  );

  const onConnect = useCallback(
    (connection: Connection) => {
      if (readOnly) return;
      setEdges((eds) =>
        addEdge(
          {
            ...connection,
            animated: true,
            style: { stroke: '#3b82f6', strokeWidth: 2 },
            markerEnd: { type: MarkerType.ArrowClosed, color: '#3b82f6' },
          },
          eds
        )
      );
    },
    [readOnly]
  );

  const onNodeClick = useCallback((_event: React.MouseEvent, node: Node<PipelineNodeData>) => {
    setSelectedNode(node);
  }, []);

  const handleConfigureNode = useCallback((nodeId: string) => {
    const node = nodes.find((n) => n.id === nodeId);
    if (node) {
      setSelectedNode(node);
      setShowConfigPanel(true);
    }
  }, [nodes]);

  const handleDeleteNode = useCallback((nodeId: string) => {
    if (window.confirm('Delete this node?')) {
      setNodes((nds) => nds.filter((n) => n.id !== nodeId));
      setEdges((eds) => eds.filter((e) => e.source !== nodeId && e.target !== nodeId));
      if (selectedNode?.id === nodeId) {
        setSelectedNode(null);
        setShowConfigPanel(false);
      }
    }
  }, [selectedNode]);

  const handleUpdateNode = useCallback((nodeId: string, updates: Partial<PipelineNodeData>) => {
    setNodes((nds) =>
      nds.map((node) =>
        node.id === nodeId
          ? {
              ...node,
              data: {
                ...node.data,
                ...updates,
              },
            }
          : node
      )
    );
  }, []);

  const handleViewData = useCallback(() => {
    if (!selectedNode) return;
    
    // 生成示例数据用于演示
    const sampleData = Array.from({ length: 100 }, (_, i) => ({
      id: i + 1,
      sample_name: `Sample_${i + 1}`,
      value: Math.random() * 100,
      category: ['A', 'B', 'C'][Math.floor(Math.random() * 3)],
      score: Math.random() * 10,
    }));
    
    setVizData(sampleData);
    setShowDataViz(true);
  }, [selectedNode]);

  const onDragOver = useCallback((event: DragEvent) => {
    event.preventDefault();
    event.dataTransfer.dropEffect = 'move';
  }, []);

  const onDrop = useCallback(
    (event: DragEvent) => {
      event.preventDefault();

      if (!reactFlowInstance || !reactFlowWrapper.current) return;

      const toolData = event.dataTransfer.getData('application/reactflow');
      if (!toolData) return;

      const tool: ToolDefinition = JSON.parse(toolData);
      const reactFlowBounds = reactFlowWrapper.current.getBoundingClientRect();
      const position = reactFlowInstance.project({
        x: event.clientX - reactFlowBounds.left,
        y: event.clientY - reactFlowBounds.top,
      });

      const newNode: Node<PipelineNodeData> = {
        id: `node_${nodeIdCounter.current++}`,
        type: 'pipelineNode',
        position,
        data: {
          label: tool.name,
          tool: tool.tool,
          version: tool.version,
          parameters: tool.defaultParameters || {},
          nodeType: tool.type === 'visualization' ? 'output' : tool.type,
          onConfigure: handleConfigureNode,
          onDelete: handleDeleteNode,
        },
      };

      setNodes((nds) => [...nds, newNode]);
    },
    [reactFlowInstance, handleConfigureNode, handleDeleteNode]
  );

  const handleSave = () => {
    // Validate before saving
    const result = validatePipeline(nodes, edges);
    setValidationResult(result);
    
    if (!result.valid) {
      setShowValidation(true);
      if (!window.confirm(`Pipeline has ${result.errors.length} error(s). Save anyway?`)) {
        return;
      }
    }
    
    if (onSave) {
      onSave(nodes, edges);
    }
    
    if (result.valid || result.warnings.length === 0) {
      alert('✅ Pipeline saved successfully!');
    }
  };

  const handleValidate = () => {
    const result = validatePipeline(nodes, edges);
    setValidationResult(result);
    setShowValidation(true);
  };

  const handleExport = () => {
    const pipelineName = prompt('Enter pipeline name:');
    if (!pipelineName) return;
    
    downloadPipelineJSON(nodes, edges, {
      name: pipelineName,
      description: `Pipeline with ${nodes.length} nodes`,
      author: 'User',
    });
  };

  const handleImport = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    try {
      const imported = await importPipeline(file);
      setNodes(imported.nodes);
      setEdges(imported.edges);
      nodeIdCounter.current = imported.nodes.length;
      alert(`✅ Imported pipeline: ${imported.metadata.name}`);
    } catch (error) {
      alert(`❌ Failed to import pipeline: ${error}`);
    }
    
    // Reset file input
    event.target.value = '';
  };

  const handleLoadTemplate = (template: PipelineExport) => {
    if (nodes.length > 0 && !window.confirm('Replace current pipeline with template?')) {
      return;
    }
    setNodes(template.nodes);
    setEdges(template.edges);
    nodeIdCounter.current = template.nodes.length;
    setShowImportExport(false);
    alert(`✅ Loaded template: ${template.metadata.name}`);
  };

  const handleClear = () => {
    if (readOnly) return;
    if (window.confirm('Clear entire pipeline? This cannot be undone.')) {
      setNodes([]);
      setEdges([]);
      setSelectedNode(null);
      setShowConfigPanel(false);
      nodeIdCounter.current = 0;
    }
  };

  const handleAutoLayout = () => {
    // Simple vertical layout
    const layoutedNodes = nodes.map((node, index) => ({
      ...node,
      position: {
        x: 250 + (index % 3) * 250,
        y: 100 + Math.floor(index / 3) * 150,
      },
    }));
    setNodes(layoutedNodes);
  };

  return (
    <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column' }}>
      {/* Toolbar */}
      {!readOnly && (
        <div
          style={{
            padding: '1rem',
            backgroundColor: '#111827',
            borderBottom: '1px solid #374151',
            display: 'flex',
            gap: '0.75rem',
            alignItems: 'center',
            flexWrap: 'wrap',
          }}
        >
          <button
            onClick={handleSave}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#059669',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 600,
              fontSize: '0.875rem',
            }}
          >
            💾 Save
          </button>

          <button
            onClick={handleValidate}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: validationResult?.valid ? '#059669' : '#f59e0b',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
              fontSize: '0.875rem',
            }}
          >
            ✓ Validate
          </button>

          <button
            onClick={handleExport}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#8b5cf6',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
              fontSize: '0.875rem',
            }}
          >
            📤 Export
          </button>

          <label
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#8b5cf6',
              color: 'white',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
              fontSize: '0.875rem',
              display: 'inline-block',
            }}
          >
            📥 Import
            <input
              type="file"
              accept=".json"
              onChange={handleImport}
              style={{ display: 'none' }}
            />
          </label>

          <button
            onClick={() => setShowImportExport(!showImportExport)}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#6366f1',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
              fontSize: '0.875rem',
            }}
          >
            📚 Templates
          </button>

          <button
            onClick={handleAutoLayout}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#3b82f6',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
              fontSize: '0.875rem',
            }}
          >
            🔄 Layout
          </button>

          <button
            onClick={() => setShowPipelineStats(true)}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#10b981',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
              fontSize: '0.875rem',
            }}
          >
            📈 Statistics
          </button>

          <button
            onClick={() => setShowPreview(true)}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#f59e0b',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
              fontSize: '0.875rem',
            }}
          >
            🔍 Preview
          </button>

          <button
            onClick={handleClear}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#dc2626',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
              fontSize: '0.875rem',
            }}
          >
            🗑️ Clear
          </button>

          <div
            style={{
              marginLeft: 'auto',
              display: 'flex',
              gap: '1rem',
              fontSize: '0.875rem',
              color: '#9ca3af',
            }}
          >
            <span>
              <strong style={{ color: '#f3f4f6' }}>{nodes.length}</strong> nodes
            </span>
            <span>
              <strong style={{ color: '#f3f4f6' }}>{edges.length}</strong> connections
            </span>
          </div>
        </div>
      )}

      {/* Main Editor Area */}
      <div style={{ flex: 1, display: 'flex', overflow: 'hidden' }}>
        {/* Toolbox */}
        {!readOnly && <ToolboxV2 onDragStart={() => {}} />}

        {/* Canvas */}
        <div ref={reactFlowWrapper} style={{ flex: 1, backgroundColor: '#0f172a' }} onDrop={onDrop} onDragOver={onDragOver}>
          <ReactFlow
            nodes={nodes}
            edges={edges}
            onNodesChange={onNodesChange}
            onEdgesChange={onEdgesChange}
            onConnect={onConnect}
            onNodeClick={onNodeClick}
            onInit={setReactFlowInstance}
            nodeTypes={nodeTypes}
            fitView
            nodesDraggable={!readOnly}
            nodesConnectable={!readOnly}
            elementsSelectable={!readOnly}
            style={{ backgroundColor: '#0f172a' }}
          >
            <Background color="#1e293b" variant={BackgroundVariant.Dots} gap={16} size={1} />
            <Controls />
            <MiniMap
              nodeColor={(node) => {
                const data = node.data as PipelineNodeData;
                const colors: Record<string, string> = {
                  input: '#10b981',
                  process: '#3b82f6',
                  filter: '#f59e0b',
                  transform: '#8b5cf6',
                  analysis: '#ef4444',
                  output: '#6b7280',
                };
                return colors[data.nodeType] || '#3b82f6';
              }}
              style={{
                backgroundColor: '#1f2937',
                border: '1px solid #374151',
              }}
              maskColor="rgba(15, 23, 42, 0.8)"
            />
          </ReactFlow>
        </div>

        {/* Config Panel */}
        {!readOnly && showConfigPanel && (
          <ConfigPanelV2
            selectedNode={selectedNode}
            onUpdate={handleUpdateNode}
            onClose={() => setShowConfigPanel(false)}
          />
        )}

        {/* View Data Button */}
        {!readOnly && selectedNode && (
          <button
            onClick={handleViewData}
            style={{
              position: 'fixed',
              bottom: '2rem',
              right: '2rem',
              padding: '0.75rem 1.5rem',
              backgroundColor: '#8b5cf6',
              color: 'white',
              border: 'none',
              borderRadius: '8px',
              cursor: 'pointer',
              fontWeight: 600,
              fontSize: '0.875rem',
              boxShadow: '0 4px 6px rgba(0,0,0,0.3)',
              zIndex: 1000,
              display: 'flex',
              alignItems: 'center',
              gap: '0.5rem',
            }}
          >
            📊 View Data
          </button>
        )}
      </div>

      {/* Validation Results Dialog */}
      {showValidation && validationResult && (
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
          onClick={() => setShowValidation(false)}
        >
          <div
            style={{
              backgroundColor: '#1f2937',
              borderRadius: '8px',
              padding: '2rem',
              maxWidth: '600px',
              maxHeight: '80vh',
              overflow: 'auto',
              border: '2px solid #374151',
            }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 style={{ margin: '0 0 1rem 0', color: '#f3f4f6' }}>
              {validationResult.valid ? '✅ Pipeline Valid' : '❌ Validation Failed'}
            </h2>

            {validationResult.errors.length > 0 && (
              <div style={{ marginBottom: '1rem' }}>
                <h3 style={{ color: '#ef4444', fontSize: '1rem', marginBottom: '0.5rem' }}>
                  ⚠️ Errors ({validationResult.errors.length})
                </h3>
                {validationResult.errors.map((issue, idx) => (
                  <div
                    key={idx}
                    style={{
                      padding: '0.75rem',
                      backgroundColor: '#7f1d1d',
                      borderRadius: '4px',
                      marginBottom: '0.5rem',
                      border: '1px solid #991b1b',
                    }}
                  >
                    <div style={{ color: '#fecaca', fontSize: '0.875rem', marginBottom: '0.25rem' }}>
                      {issue.message}
                    </div>
                    {issue.suggestion && (
                      <div style={{ color: '#fca5a5', fontSize: '0.75rem', fontStyle: 'italic' }}>
                        💡 {issue.suggestion}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            )}

            {validationResult.warnings.length > 0 && (
              <div style={{ marginBottom: '1rem' }}>
                <h3 style={{ color: '#f59e0b', fontSize: '1rem', marginBottom: '0.5rem' }}>
                  ⚠️ Warnings ({validationResult.warnings.length})
                </h3>
                {validationResult.warnings.map((issue, idx) => (
                  <div
                    key={idx}
                    style={{
                      padding: '0.75rem',
                      backgroundColor: '#78350f',
                      borderRadius: '4px',
                      marginBottom: '0.5rem',
                      border: '1px solid #92400e',
                    }}
                  >
                    <div style={{ color: '#fde68a', fontSize: '0.875rem', marginBottom: '0.25rem' }}>
                      {issue.message}
                    </div>
                    {issue.suggestion && (
                      <div style={{ color: '#fcd34d', fontSize: '0.75rem', fontStyle: 'italic' }}>
                        💡 {issue.suggestion}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            )}

            {validationResult.valid && validationResult.warnings.length === 0 && (
              <div style={{ padding: '1rem', backgroundColor: '#064e3b', borderRadius: '4px', color: '#6ee7b7' }}>
                ✓ Pipeline structure is valid with no issues detected!
              </div>
            )}

            <button
              onClick={() => setShowValidation(false)}
              style={{
                marginTop: '1rem',
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

      {/* Templates Dialog */}
      {showImportExport && (
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
          onClick={() => setShowImportExport(false)}
        >
          <div
            style={{
              backgroundColor: '#1f2937',
              borderRadius: '8px',
              padding: '2rem',
              maxWidth: '600px',
              maxHeight: '80vh',
              overflow: 'auto',
              border: '2px solid #374151',
            }}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 style={{ margin: '0 0 1rem 0', color: '#f3f4f6' }}>📚 Pipeline Templates</h2>

            {getPipelineTemplates().map((template, idx) => (
              <div
                key={idx}
                style={{
                  padding: '1rem',
                  backgroundColor: '#111827',
                  borderRadius: '6px',
                  marginBottom: '1rem',
                  border: '1px solid #374151',
                }}
              >
                <h3 style={{ margin: '0 0 0.5rem 0', color: '#f3f4f6', fontSize: '1rem' }}>
                  {template.metadata.name}
                </h3>
                <p style={{ margin: '0 0 0.5rem 0', color: '#9ca3af', fontSize: '0.875rem' }}>
                  {template.metadata.description}
                </p>
                <div style={{ fontSize: '0.75rem', color: '#6b7280', marginBottom: '0.75rem' }}>
                  {template.nodes.length} nodes • {template.edges.length} connections
                </div>
                <button
                  onClick={() => handleLoadTemplate(template)}
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
                  Load Template
                </button>
              </div>
            ))}

            <button
              onClick={() => setShowImportExport(false)}
              style={{
                marginTop: '1rem',
                padding: '0.75rem 1.5rem',
                backgroundColor: '#374151',
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

      {/* Pipeline Statistics Modal */}
      <PipelineStatisticsModal
        isOpen={showPipelineStats}
        onClose={() => setShowPipelineStats(false)}
        pipelineName="Custom Pipeline"
        nodes={nodes}
        edges={edges}
      />

      {/* Data Visualization Modal */}
      <DataVisualizationModal
        isOpen={showDataViz}
        onClose={() => setShowDataViz(false)}
        title={selectedNode?.data.label || 'Node Data'}
        data={vizData}
      />

      {/* Pipeline Preview Modal */}
      <PipelinePreviewModal
        isOpen={showPreview}
        onClose={() => setShowPreview(false)}
        nodes={nodes}
        edges={edges}
        onExecute={() => {
          // TODO: 实际的管道执行逻辑
          console.log('Executing pipeline...');
        }}
      />
    </div>
  );
};

export default EnhancedPipelineEditor;
