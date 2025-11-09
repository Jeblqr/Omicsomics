import { useState, useCallback, useRef } from 'react';
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
} from 'reactflow';
import 'reactflow/dist/style.css';

interface PipelineEditorProps {
  initialNodes?: Node[];
  initialEdges?: Edge[];
  onSave?: (nodes: Node[], edges: Edge[]) => void;
  readOnly?: boolean;
}

// Available node types for the pipeline
const nodeTypes = [
  { type: 'input', label: 'Input', color: '#4caf50' },
  { type: 'process', label: 'Process', color: '#2196f3' },
  { type: 'filter', label: 'Filter', color: '#ff9800' },
  { type: 'transform', label: 'Transform', color: '#9c27b0' },
  { type: 'analysis', label: 'Analysis', color: '#f44336' },
  { type: 'output', label: 'Output', color: '#607d8b' },
];

const PipelineEditor = ({ 
  initialNodes = [], 
  initialEdges = [], 
  onSave,
  readOnly = false 
}: PipelineEditorProps) => {
  const [nodes, setNodes] = useState<Node[]>(initialNodes);
  const [edges, setEdges] = useState<Edge[]>(initialEdges);
  const [selectedNodeType, setSelectedNodeType] = useState('process');
  const nodeIdCounter = useRef(initialNodes.length);

  const onNodesChange = useCallback(
    (changes: NodeChange[]) => {
      if (readOnly) return;
      setNodes((nds) => applyNodeChanges(changes, nds));
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
            markerEnd: { type: MarkerType.ArrowClosed },
          },
          eds
        )
      );
    },
    [readOnly]
  );

  const addNode = () => {
    if (readOnly) return;
    
    const nodeTypeInfo = nodeTypes.find(nt => nt.type === selectedNodeType) || nodeTypes[1];
    const newNode: Node = {
      id: `node_${nodeIdCounter.current++}`,
      type: 'default',
      position: { x: 250, y: 50 + nodes.length * 100 },
      data: { 
        label: `${nodeTypeInfo.label} Step`,
        nodeType: selectedNodeType,
      },
      style: {
        background: nodeTypeInfo.color,
        color: 'white',
        border: '2px solid #222',
        borderRadius: '8px',
        padding: '10px',
      },
    };

    setNodes((nds) => [...nds, newNode]);
  };

  const handleSave = () => {
    if (onSave) {
      onSave(nodes, edges);
    }
  };

  const clearPipeline = () => {
    if (readOnly) return;
    if (window.confirm('Are you sure you want to clear the entire pipeline?')) {
      setNodes([]);
      setEdges([]);
      nodeIdCounter.current = 0;
    }
  };

  return (
    <div style={{ width: '100%', height: '600px', border: '1px solid #ddd', borderRadius: '8px', overflow: 'hidden' }}>
      {!readOnly && (
        <div style={{ 
          padding: '1rem', 
          backgroundColor: '#f5f5f5', 
          borderBottom: '1px solid #ddd',
          display: 'flex',
          gap: '1rem',
          alignItems: 'center',
          flexWrap: 'wrap',
        }}>
          <div style={{ display: 'flex', gap: '0.5rem', alignItems: 'center' }}>
            <label htmlFor="node-type" style={{ fontWeight: 500 }}>Node Type:</label>
            <select
              id="node-type"
              value={selectedNodeType}
              onChange={(e) => setSelectedNodeType(e.target.value)}
              style={{
                padding: '0.5rem',
                borderRadius: '4px',
                border: '1px solid #ccc',
              }}
            >
              {nodeTypes.map((nt) => (
                <option key={nt.type} value={nt.type}>
                  {nt.label}
                </option>
              ))}
            </select>
          </div>

          <button
            onClick={addNode}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#2196f3',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
            }}
          >
            + Add Node
          </button>

          <button
            onClick={handleSave}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#4caf50',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
            }}
          >
            Save Pipeline
          </button>

          <button
            onClick={clearPipeline}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#f44336',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: 500,
            }}
          >
            Clear
          </button>

          <div style={{ marginLeft: 'auto', fontSize: '0.9rem', color: '#666' }}>
            Nodes: {nodes.length} | Edges: {edges.length}
          </div>
        </div>
      )}

      <div style={{ width: '100%', height: readOnly ? '600px' : '500px' }}>
        <ReactFlow
          nodes={nodes}
          edges={edges}
          onNodesChange={onNodesChange}
          onEdgesChange={onEdgesChange}
          onConnect={onConnect}
          fitView
          nodesDraggable={!readOnly}
          nodesConnectable={!readOnly}
          elementsSelectable={!readOnly}
        >
          <Background />
          <Controls />
        </ReactFlow>
      </div>
    </div>
  );
};

export default PipelineEditor;
