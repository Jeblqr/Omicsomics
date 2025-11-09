import { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { Node, Edge } from 'reactflow';
import api from '../../lib/api';
import PipelineEditor from '../../components/PipelineEditor';
import { useProjectsContext } from '../../contexts/ProjectsContext';
import { ProjectSwitcher } from '../../components/ProjectSwitcher';

interface CustomPipeline {
  id: number;
  name: string;
  description: string | null;
  category: string | null;
  is_public: boolean;
  owner_id: number;
  definition: {
    nodes: Node[];
    edges: Edge[];
    parameters?: Record<string, any>;
  };
  created_at: string;
  updated_at: string;
}

const CustomPipelinesPage = () => {
  const navigate = useNavigate();
  const { currentProject } = useProjectsContext();
  const [pipelines, setPipelines] = useState<CustomPipeline[]>([]);
  const [loading, setLoading] = useState(true);
  const [showEditor, setShowEditor] = useState(false);
  const [editingPipeline, setEditingPipeline] = useState<CustomPipeline | null>(null);
  const [pipelineName, setPipelineName] = useState('');
  const [pipelineDescription, setPipelineDescription] = useState('');
  const [pipelineCategory, setPipelineCategory] = useState('custom');
  const [isPublic, setIsPublic] = useState(false);
  const [selectedPipelines, setSelectedPipelines] = useState<number[]>([]);

  useEffect(() => {
    if (currentProject) {
      loadPipelines();
    } else {
      setPipelines([]);
      setLoading(false);
    }
  }, [currentProject]);

  const loadPipelines = async () => {
    if (!currentProject) return;
    
    try {
      const response = await api.get(`/custom-pipelines/?project_id=${currentProject.id}`);
      setPipelines(response.data);
    } catch (error) {
      console.error('Failed to load pipelines:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleSavePipeline = async (nodes: Node[], edges: Edge[]) => {
    if (!pipelineName.trim()) {
      alert('Please enter a pipeline name');
      return;
    }

    if (!currentProject) {
      alert('Please select a project first');
      return;
    }

    try {
      const definition = {
        nodes: nodes.map((node) => ({
          id: node.id,
          type: node.data.nodeType || 'process',
          label: node.data.label,
          data: node.data || {},
          position: node.position || { x: 0, y: 0 },
        })),
        edges: edges.map((edge) => ({
          id: edge.id,
          source: edge.source,
          target: edge.target,
          sourceHandle: edge.sourceHandle || null,
          targetHandle: edge.targetHandle || null,
        })),
        parameters: {},
      };

      const payload = {
        name: pipelineName,
        description: pipelineDescription || 'Custom pipeline',
        category: pipelineCategory,
        is_public: isPublic,
        definition,
      };

      if (editingPipeline && editingPipeline.id > 0) {
        await api.put(`/custom-pipelines/${editingPipeline.id}`, payload);
        alert('Pipeline updated successfully ✅');
      } else {
        await api.post('/custom-pipelines/', payload);
        alert('Pipeline created successfully ✅');
      }

      setShowEditor(false);
      setEditingPipeline(null);
      setPipelineName('');
      setPipelineDescription('');
      setPipelineCategory('custom');
      setIsPublic(false);
      loadPipelines();
    } catch (error) {
      console.error('Failed to save pipeline:', error);
      alert('Failed to save pipeline');
    }
  };

  const handleDeletePipeline = async (id: number) => {
    if (!window.confirm('Are you sure you want to delete this pipeline?')) {
      return;
    }

    try {
      await api.delete(`/custom-pipelines/${id}`);
      alert('Pipeline deleted successfully');
      loadPipelines();
    } catch (error) {
      console.error('Failed to delete pipeline:', error);
      alert('Failed to delete pipeline');
    }
  };

  const handleEditPipeline = (pipeline: CustomPipeline) => {
    setEditingPipeline(pipeline);
    setPipelineName(pipeline.name);
    setPipelineDescription(pipeline.description || '');
    setPipelineCategory(pipeline.category || 'custom');
    setIsPublic(pipeline.is_public);
    setShowEditor(true);
  };

  const handleNewPipeline = () => {
    setEditingPipeline(null);
    setPipelineName('');
    setPipelineDescription('');
    setPipelineCategory('custom');
    setIsPublic(false);
    setShowEditor(true);
  };

  const togglePipelineSelection = (id: number) => {
    setSelectedPipelines((prev) =>
      prev.includes(id) ? prev.filter((pid) => pid !== id) : [...prev, id]
    );
  };

  const handleMergePipelines = async () => {
    if (selectedPipelines.length < 2) {
      alert('Please select at least 2 pipelines to merge');
      return;
    }

    try {
      const response = await api.post('/custom-pipelines/merge', {
        pipeline_ids: selectedPipelines,
      });

      const merged = response.data;
      setEditingPipeline({
        id: 0,
        name: 'Merged Pipeline',
        description: 'Merged from selected pipelines',
        category: 'custom',
        is_public: false,
        owner_id: 0,
        definition: merged,
        created_at: new Date().toISOString(),
        updated_at: new Date().toISOString(),
      });
      setPipelineName('Merged Pipeline');
      setPipelineDescription('Merged from selected pipelines');
      setPipelineCategory('custom');
      setIsPublic(false);
      setSelectedPipelines([]);
      setShowEditor(true);
    } catch (error) {
      console.error('Failed to merge pipelines:', error);
      alert('Failed to merge pipelines');
    }
  };

  const convertToReactFlowNodes = (definition: CustomPipeline['definition']): Node[] => {
    return definition.nodes.map((node: any) => ({
      id: node.id,
      type: 'default',
      position: node.position || { x: 0, y: 0 },
      data: {
        label: node.label,
        nodeType: node.type,
      },
      style: {
        background: getNodeColor(node.type),
        color: 'white',
        border: '2px solid #222',
        borderRadius: '8px',
        padding: '10px',
      },
    }));
  };

  const convertToReactFlowEdges = (definition: CustomPipeline['definition']): Edge[] => {
    return definition.edges.map((edge: any) => ({
      id: edge.id,
      source: edge.source,
      target: edge.target,
      animated: true,
    }));
  };

  const getNodeColor = (type: string): string => {
    const colors: Record<string, string> = {
      input: '#4caf50',
      process: '#2196f3',
      filter: '#ff9800',
      transform: '#9c27b0',
      analysis: '#f44336',
      output: '#607d8b',
    };
    return colors[type] || '#2196f3';
  };

  if (loading) {
    return (
      <div style={{ padding: '2rem', textAlign: 'center' }}>
        Loading pipelines...
      </div>
    );
  }

  if (showEditor) {
    return (
      <div style={{ padding: '2rem' }}>
        <div style={{ marginBottom: '2rem' }}>
          <button
            onClick={() => {
              setShowEditor(false);
              setEditingPipeline(null);
            }}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#666',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              marginBottom: '1rem',
            }}
          >
            ← Back to List
          </button>

          <h2 style={{ marginBottom: '1rem' }}>
            {editingPipeline && editingPipeline.id > 0 ? 'Edit Pipeline' : 'Create New Pipeline'}
          </h2>

          <div style={{ marginBottom: '1rem', display: 'flex', flexDirection: 'column', gap: '1rem' }}>
            <div>
              <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: 500 }}>
                Pipeline Name:
              </label>
              <input
                type="text"
                value={pipelineName}
                onChange={(e) => setPipelineName(e.target.value)}
                style={{
                  width: '100%',
                  padding: '0.5rem',
                  borderRadius: '4px',
                  border: '1px solid #ccc',
                }}
                placeholder="Enter pipeline name"
              />
            </div>

            <div>
              <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: 500 }}>
                Description:
              </label>
              <textarea
                value={pipelineDescription}
                onChange={(e) => setPipelineDescription(e.target.value)}
                style={{
                  width: '100%',
                  padding: '0.5rem',
                  borderRadius: '4px',
                  border: '1px solid #ccc',
                  minHeight: '80px',
                }}
                placeholder="Enter pipeline description"
              />
            </div>

            <div style={{ display: 'flex', gap: '2rem' }}>
              <div>
                <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: 500 }}>
                  Category:
                </label>
                <select
                  value={pipelineCategory}
                  onChange={(e) => setPipelineCategory(e.target.value)}
                  style={{
                    padding: '0.5rem',
                    borderRadius: '4px',
                    border: '1px solid #ccc',
                  }}
                >
                  <option value="custom">Custom</option>
                  <option value="genomics">Genomics</option>
                  <option value="proteomics">Proteomics</option>
                  <option value="metabolomics">Metabolomics</option>
                  <option value="multiomics">Multi-omics</option>
                </select>
              </div>

              <div>
                <label style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', cursor: 'pointer' }}>
                  <input
                    type="checkbox"
                    checked={isPublic}
                    onChange={(e) => setIsPublic(e.target.checked)}
                    style={{ width: '20px', height: '20px' }}
                  />
                  <span style={{ fontWeight: 500 }}>Make Public</span>
                </label>
              </div>
            </div>
          </div>
        </div>

        <PipelineEditor
          initialNodes={
            editingPipeline
              ? convertToReactFlowNodes(editingPipeline.definition)
              : []
          }
          initialEdges={
            editingPipeline
              ? convertToReactFlowEdges(editingPipeline.definition)
              : []
          }
          onSave={handleSavePipeline}
        />
      </div>
    );
  }

  return (
    <div>
      <ProjectSwitcher />
      <div style={{ padding: '2rem' }}>
        {!currentProject && (
          <div style={{
            padding: '1.5rem',
            backgroundColor: '#fff3cd',
            border: '1px solid #ffc107',
            borderRadius: '4px',
            marginBottom: '1rem',
          }}>
            <strong>⚠️ No project selected</strong>
            <p style={{ marginBottom: 0, marginTop: '0.5rem' }}>
              Please select a project from the dropdown above to view and create custom pipelines.
            </p>
          </div>
        )}

        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '2rem' }}>
          <h1>Custom Pipelines</h1>
        <div style={{ display: 'flex', gap: '1rem' }}>
          <button
            onClick={handleMergePipelines}
            disabled={selectedPipelines.length < 2}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: selectedPipelines.length >= 2 ? '#ff9800' : '#ccc',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: selectedPipelines.length >= 2 ? 'pointer' : 'not-allowed',
              fontWeight: 500,
            }}
          >
            Merge Selected ({selectedPipelines.length})
          </button>
          <button
            onClick={handleNewPipeline}
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
            + Create New Pipeline
          </button>
        </div>
      </div>

      {pipelines.length === 0 ? (
        <div style={{ textAlign: 'center', padding: '3rem', color: '#666' }}>
          <p>No custom pipelines yet. Create your first pipeline!</p>
        </div>
      ) : (
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(350px, 1fr))', gap: '1.5rem' }}>
          {pipelines.map((pipeline) => (
            <div
              key={pipeline.id}
              style={{
                border: '1px solid #ddd',
                borderRadius: '8px',
                padding: '1.5rem',
                backgroundColor: selectedPipelines.includes(pipeline.id) ? '#e3f2fd' : 'white',
                cursor: 'pointer',
              }}
            >
              <div style={{ display: 'flex', alignItems: 'flex-start', gap: '1rem', marginBottom: '1rem' }}>
                <input
                  type="checkbox"
                  checked={selectedPipelines.includes(pipeline.id)}
                  onChange={() => togglePipelineSelection(pipeline.id)}
                  style={{ marginTop: '0.25rem', width: '20px', height: '20px' }}
                />
                <div style={{ flex: 1 }}>
                  <h3 style={{ marginBottom: '0.5rem' }}>{pipeline.name}</h3>
                  {pipeline.description && (
                    <p style={{ color: '#666', fontSize: '0.9rem', marginBottom: '0.5rem' }}>
                      {pipeline.description}
                    </p>
                  )}
                  <div style={{ display: 'flex', gap: '0.5rem', flexWrap: 'wrap', fontSize: '0.85rem' }}>
                    {pipeline.category && (
                      <span
                        style={{
                          backgroundColor: '#e0e0e0',
                          padding: '0.25rem 0.5rem',
                          borderRadius: '4px',
                        }}
                      >
                        {pipeline.category}
                      </span>
                    )}
                    {pipeline.is_public && (
                      <span
                        style={{
                          backgroundColor: '#4caf50',
                          color: 'white',
                          padding: '0.25rem 0.5rem',
                          borderRadius: '4px',
                        }}
                      >
                        Public
                      </span>
                    )}
                    <span
                      style={{
                        backgroundColor: '#2196f3',
                        color: 'white',
                        padding: '0.25rem 0.5rem',
                        borderRadius: '4px',
                      }}
                    >
                      {pipeline.definition.nodes.length} nodes
                    </span>
                  </div>
                </div>
              </div>

              <div style={{ display: 'flex', gap: '0.5rem', marginTop: '1rem' }}>
                <button
                  onClick={() => handleEditPipeline(pipeline)}
                  style={{
                    flex: 1,
                    padding: '0.5rem',
                    backgroundColor: '#2196f3',
                    color: 'white',
                    border: 'none',
                    borderRadius: '4px',
                    cursor: 'pointer',
                  }}
                >
                  Edit
                </button>
                <button
                  onClick={() => handleDeletePipeline(pipeline.id)}
                  style={{
                    flex: 1,
                    padding: '0.5rem',
                    backgroundColor: '#f44336',
                    color: 'white',
                    border: 'none',
                    borderRadius: '4px',
                    cursor: 'pointer',
                  }}
                >
                  Delete
                </button>
              </div>
            </div>
          ))}
        </div>
      )}
      </div>
    </div>
  );
};

export default CustomPipelinesPage;
