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
    if (onSave) {
      onSave(nodes, edges);
    }
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
            💾 Save Pipeline
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
            🔄 Auto Layout
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
            🗑️ Clear All
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
      </div>
    </div>
  );
};

export default EnhancedPipelineEditor;
