import { useState, useEffect } from 'react';
import { ProjectSwitcher } from '../../components/ProjectSwitcher';
import { useProjectsContext } from '../../contexts/ProjectsContext';
import api from '../../lib/api';
import LoadingView from '../../components/LoadingView';
import { PipelineSelector, PipelineSelection } from '../../components/runs/PipelineSelector';
import { DataFileSelector } from '../../components/runs/DataFileSelector';
import { ParameterConfig } from '../../components/runs/ParameterConfig';

interface Run {
  id: number;
  name: string;
  description: string;
  status: string;
  pipeline_type?: string;
  pipeline_template_id?: string;
  custom_pipeline_id?: number;
  parameters?: Record<string, any>;
  input_files?: number[];
  progress?: number;
  project_id: number;
  owner_id: number;
  started_at: string | null;
  finished_at: string | null;
  created_at: string;
  updated_at: string;
}

const RunsPage = () => {
  const { currentProject } = useProjectsContext();
  const [runs, setRuns] = useState<Run[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [showForm, setShowForm] = useState(false);
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [pipelineSelection, setPipelineSelection] = useState<PipelineSelection>({ type: 'template' });
  const [selectedFiles, setSelectedFiles] = useState<number[]>([]);
  const [parameters, setParameters] = useState<Record<string, any>>({});
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [error, setError] = useState('');

  useEffect(() => {
    if (currentProject) {
      fetchRuns();
    } else {
      setRuns([]);
    }
  }, [currentProject]);

  const fetchRuns = async () => {
    if (!currentProject) return;
    
    setIsLoading(true);
    setError('');
    try {
      const response = await api.get<Run[]>(`/runs/?project_id=${currentProject.id}`);
      setRuns(response.data);
    } catch (err) {
      console.error('Failed to fetch runs:', err);
      const error = err as { response?: { data?: { detail?: string } }; message?: string };
      setError(error.response?.data?.detail || 'Failed to load runs');
    } finally {
      setIsLoading(false);
    }
  };

  const handleCreateRun = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!currentProject) return;

    // Validation
    if (!name.trim()) {
      setError('Please enter a run name');
      return;
    }
    
    if (!pipelineSelection.templateId && !pipelineSelection.customId) {
      setError('Please select a pipeline');
      return;
    }

    setIsSubmitting(true);
    setError('');

    try {
      const payload = {
        name,
        description,
        project_id: currentProject.id,
        pipeline_type: pipelineSelection.type,
        pipeline_template_id: pipelineSelection.templateId || null,
        custom_pipeline_id: pipelineSelection.customId || null,
        parameters,
        input_files: selectedFiles,
        input_mapping: {},
        auto_start: false,
        priority: 0,
      };

      await api.post('/runs/', payload);

      // Reset form
      setName('');
      setDescription('');
      setPipelineSelection({ type: 'template' });
      setSelectedFiles([]);
      setParameters({});
      setShowForm(false);
      fetchRuns();
    } catch (err) {
      const error = err as { response?: { data?: { detail?: string } }; message?: string };
      setError(error.response?.data?.detail || 'Failed to create run');
    } finally {
      setIsSubmitting(false);
    }
  };

  const getStatusBadge = (status: string) => {
    const colors: Record<string, string> = {
      pending: '#ffc107',
      running: '#007bff',
      completed: '#28a745',
      failed: '#dc3545',
    };
    
    return (
      <span style={{
        padding: '0.25rem 0.75rem',
        borderRadius: '12px',
        fontSize: '0.85rem',
        fontWeight: 500,
        backgroundColor: colors[status] || '#6c757d',
        color: 'white',
      }}>
        {status.toUpperCase()}
      </span>
    );
  };

  return (
    <section>
      <ProjectSwitcher />
      
      <div style={{ padding: '2rem' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
          <h2 style={{ color: '#212529' }}>Pipeline Runs</h2>
          {currentProject && (
            <button
              onClick={() => setShowForm(!showForm)}
              style={{
                padding: '0.5rem 1rem',
                background: showForm ? '#6c757d' : '#007bff',
                color: 'white',
                border: 'none',
                borderRadius: '4px',
                cursor: 'pointer',
              }}
            >
              {showForm ? 'Cancel' : '+ New Run'}
            </button>
          )}
        </div>
        
        <p style={{ color: '#6c757d' }}>Track workflow executions, logs, and outputs.</p>

        {!currentProject && (
          <div style={{
            padding: '1.5rem',
            backgroundColor: '#fff3cd',
            border: '1px solid #ffc107',
            borderRadius: '4px',
            marginTop: '1rem',
          }}>
            <strong>⚠️ No project selected</strong>
            <p style={{ marginBottom: 0, marginTop: '0.5rem' }}>
              Please select a project from the dropdown above to view and create runs.
            </p>
          </div>
        )}

        {showForm && currentProject && (
          <div style={{
            background: '#f8f9fa',
            padding: '1.5rem',
            borderRadius: '8px',
            marginBottom: '1.5rem',
            border: '1px solid #dee2e6',
          }}>
            <h3 style={{ marginTop: 0, color: '#212529' }}>Create New Run</h3>
            <form onSubmit={handleCreateRun}>
              <div style={{ marginBottom: '1.5rem' }}>
                <label htmlFor="run-name" style={{ display: 'block', marginBottom: '0.5rem', color: '#212529', fontWeight: 500 }}>
                  Run Name *
                </label>
                <input
                  id="run-name"
                  type="text"
                  value={name}
                  onChange={(e) => setName(e.target.value)}
                  required
                  placeholder="e.g., RNA-seq Analysis Run 1"
                  style={{
                    width: '100%',
                    padding: '0.5rem',
                    borderRadius: '4px',
                    border: '1px solid #ced4da',
                    color: '#212529',
                    backgroundColor: '#ffffff',
                    fontSize: '1rem',
                  }}
                />
              </div>

              <div style={{ marginBottom: '1.5rem' }}>
                <label htmlFor="run-desc" style={{ display: 'block', marginBottom: '0.5rem', color: '#212529', fontWeight: 500 }}>
                  Description
                </label>
                <textarea
                  id="run-desc"
                  value={description}
                  onChange={(e) => setDescription(e.target.value)}
                  rows={2}
                  placeholder="Optional description of this run"
                  style={{
                    width: '100%',
                    padding: '0.5rem',
                    borderRadius: '4px',
                    border: '1px solid #ced4da',
                    color: '#212529',
                    backgroundColor: '#ffffff',
                    fontSize: '1rem',
                    fontFamily: 'inherit',
                  }}
                />
              </div>

              <div style={{ marginBottom: '1.5rem' }}>
                <PipelineSelector
                  projectId={currentProject.id}
                  onSelect={setPipelineSelection}
                  value={pipelineSelection}
                />
              </div>

              <div style={{ marginBottom: '1.5rem' }}>
                <DataFileSelector
                  projectId={currentProject.id}
                  onSelect={setSelectedFiles}
                  value={selectedFiles}
                />
              </div>

              <div style={{ marginBottom: '1.5rem' }}>
                <ParameterConfig
                  parameters={parameters}
                  onChange={setParameters}
                />
              </div>

              {error && (
                <div style={{ color: '#dc3545', marginBottom: '1rem', padding: '0.75rem', backgroundColor: '#f8d7da', borderRadius: '4px', border: '1px solid #f5c6cb' }} role="alert">
                  {error}
                </div>
              )}

              <button
                type="submit"
                disabled={isSubmitting}
                style={{
                  padding: '0.5rem 1.5rem',
                  background: isSubmitting ? '#6c757d' : '#28a745',
                  color: 'white',
                  border: 'none',
                  borderRadius: '4px',
                  cursor: isSubmitting ? 'not-allowed' : 'pointer',
                  fontSize: '1rem',
                  fontWeight: 500,
                }}
              >
                {isSubmitting ? 'Creating...' : 'Create Run'}
              </button>
            </form>
          </div>
        )}

        {isLoading && <LoadingView />}
        
        {!isLoading && currentProject && (
          <>
            {error && !showForm && (
              <div style={{ color: '#dc3545', padding: '1rem', backgroundColor: '#f8d7da', borderRadius: '4px', marginBottom: '1rem' }}>
                {error}
              </div>
            )}

            {runs.length > 0 ? (
              <div style={{ 
                overflowX: 'auto',
                marginTop: '1rem',
                borderRadius: '8px',
                border: '1px solid #dee2e6',
                backgroundColor: '#ffffff',
              }}>
                <table style={{ 
                  width: '100%',
                  borderCollapse: 'collapse',
                  minWidth: '800px',
                }}>
                  <thead>
                    <tr style={{ backgroundColor: '#f8f9fa' }}>
                      <th style={{ 
                        padding: '1rem',
                        textAlign: 'left',
                        borderBottom: '2px solid #dee2e6',
                        fontWeight: 600,
                        color: '#212529',
                      }}>Run Name</th>
                      <th style={{ 
                        padding: '1rem',
                        textAlign: 'left',
                        borderBottom: '2px solid #dee2e6',
                        fontWeight: 600,
                        color: '#212529',
                      }}>Status</th>
                      <th style={{ 
                        padding: '1rem',
                        textAlign: 'left',
                        borderBottom: '2px solid #dee2e6',
                        fontWeight: 600,
                        color: '#212529',
                      }}>Description</th>
                      <th style={{ 
                        padding: '1rem',
                        textAlign: 'left',
                        borderBottom: '2px solid #dee2e6',
                        fontWeight: 600,
                        color: '#212529',
                      }}>Created</th>
                      <th style={{ 
                        padding: '1rem',
                        textAlign: 'left',
                        borderBottom: '2px solid #dee2e6',
                        fontWeight: 600,
                        color: '#212529',
                      }}>Progress</th>
                      <th style={{ 
                        padding: '1rem',
                        textAlign: 'center',
                        borderBottom: '2px solid #dee2e6',
                        fontWeight: 600,
                        color: '#212529',
                      }}>Actions</th>
                    </tr>
                  </thead>
                  <tbody>
                    {runs.map((run) => (
                      <tr key={run.id} style={{ 
                        borderBottom: '1px solid #dee2e6',
                        transition: 'background-color 0.2s',
                      }}>
                        <td style={{ padding: '1rem', color: '#212529' }}>
                          <div style={{ fontWeight: 600, marginBottom: '0.25rem' }}>{run.name}</div>
                          <div style={{ fontSize: '0.85rem', color: '#6c757d' }}>
                            ID: {run.id} | Pipeline: {run.pipeline_type || '—'}
                          </div>
                        </td>
                        <td style={{ padding: '1rem' }}>{getStatusBadge(run.status)}</td>
                        <td style={{ padding: '1rem', color: '#212529', maxWidth: '200px' }}>
                          {run.description || <span style={{ color: '#6c757d' }}>—</span>}
                        </td>
                        <td style={{ padding: '1rem', color: '#212529' }}>
                          <div>{new Date(run.created_at).toLocaleDateString()}</div>
                          <div style={{ fontSize: '0.85rem', color: '#6c757d' }}>
                            {new Date(run.created_at).toLocaleTimeString()}
                          </div>
                        </td>
                        <td style={{ padding: '1rem', color: '#212529' }}>
                          {run.progress !== undefined ? (
                            <div>
                              <div style={{ fontSize: '0.85rem', marginBottom: '0.25rem' }}>
                                {run.progress}%
                              </div>
                              <div style={{ 
                                width: '100px',
                                height: '6px',
                                backgroundColor: '#e9ecef',
                                borderRadius: '3px',
                                overflow: 'hidden',
                              }}>
                                <div style={{
                                  width: `${run.progress}%`,
                                  height: '100%',
                                  backgroundColor: run.status === 'completed' ? '#28a745' : '#007bff',
                                  transition: 'width 0.3s',
                                }}></div>
                              </div>
                            </div>
                          ) : (
                            <span style={{ color: '#6c757d' }}>—</span>
                          )}
                        </td>
                        <td style={{ padding: '1rem', textAlign: 'center' }}>
                          <div style={{ display: 'flex', gap: '0.5rem', justifyContent: 'center', flexWrap: 'wrap' }}>
                            {run.status === 'pending' && (
                              <button
                                onClick={async () => {
                                  try {
                                    await api.post(`/runs/${run.id}/start`);
                                    alert('Run started successfully! 🚀');
                                    fetchRuns();
                                  } catch (err) {
                                    console.error('Failed to start run:', err);
                                    alert('Failed to start run ❌');
                                  }
                                }}
                                style={{
                                  padding: '0.4rem 0.8rem',
                                  background: '#28a745',
                                  color: 'white',
                                  border: 'none',
                                  borderRadius: '4px',
                                  cursor: 'pointer',
                                  fontSize: '0.85rem',
                                  fontWeight: 500,
                                }}
                                title="Start run execution"
                              >
                                ▶️ Start
                              </button>
                            )}
                            {run.status === 'running' && (
                              <button
                                onClick={async () => {
                                  if (window.confirm(`Stop run "${run.name}"?`)) {
                                    try {
                                      await api.post(`/runs/${run.id}/stop`);
                                      alert('Run stopped ⏸️');
                                      fetchRuns();
                                    } catch (err) {
                                      console.error('Failed to stop run:', err);
                                      alert('Failed to stop run ❌');
                                    }
                                  }
                                }}
                                style={{
                                  padding: '0.4rem 0.8rem',
                                  background: '#ffc107',
                                  color: '#212529',
                                  border: 'none',
                                  borderRadius: '4px',
                                  cursor: 'pointer',
                                  fontSize: '0.85rem',
                                  fontWeight: 500,
                                }}
                                title="Stop run execution"
                              >
                                ⏸️ Stop
                              </button>
                            )}
                            <button
                              onClick={async () => {
                                try {
                                  const response = await api.get(`/runs/${run.id}`);
                                  const runDetails = response.data;
                                  const logs = runDetails.logs || 'No logs available yet.';
                                  const errorMsg = runDetails.error_message || '';
                                  
                                  const logContent = `
=== RUN LOGS ===
Run: ${run.name}
Status: ${run.status}
Created: ${new Date(run.created_at).toLocaleString()}
${run.started_at ? `Started: ${new Date(run.started_at).toLocaleString()}` : ''}
${run.finished_at ? `Finished: ${new Date(run.finished_at).toLocaleString()}` : ''}

${errorMsg ? `ERROR:\n${errorMsg}\n\n` : ''}LOGS:
${logs}
                                  `.trim();
                                  
                                  // Create modal
                                  const modal = document.createElement('div');
                                  modal.style.cssText = 'position: fixed; top: 0; left: 0; right: 0; bottom: 0; background: rgba(0,0,0,0.5); display: flex; align-items: center; justify-content: center; z-index: 9999;';
                                  
                                  const content = document.createElement('div');
                                  content.style.cssText = 'background: white; padding: 2rem; border-radius: 8px; max-width: 800px; max-height: 80vh; overflow: auto; color: #212529;';
                                  
                                  content.innerHTML = `
                                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 1rem;">
                                      <h3 style="margin: 0; color: #212529;">Run Logs - ${run.name}</h3>
                                      <button id="closeModal" style="background: #dc3545; color: white; border: none; padding: 0.5rem 1rem; border-radius: 4px; cursor: pointer;">✖ Close</button>
                                    </div>
                                    <pre style="background: #f8f9fa; padding: 1rem; border-radius: 4px; overflow-x: auto; white-space: pre-wrap; word-wrap: break-word; color: #212529; font-size: 0.9rem; line-height: 1.5;">${logContent}</pre>
                                  `;
                                  
                                  modal.appendChild(content);
                                  document.body.appendChild(modal);
                                  
                                  document.getElementById('closeModal')?.addEventListener('click', () => {
                                    document.body.removeChild(modal);
                                  });
                                  
                                  modal.addEventListener('click', (e) => {
                                    if (e.target === modal) {
                                      document.body.removeChild(modal);
                                    }
                                  });
                                } catch (err) {
                                  console.error('Failed to fetch logs:', err);
                                  alert('Failed to fetch logs ❌');
                                }
                              }}
                              style={{
                                padding: '0.4rem 0.8rem',
                                background: '#17a2b8',
                                color: 'white',
                                border: 'none',
                                borderRadius: '4px',
                                cursor: 'pointer',
                                fontSize: '0.85rem',
                                fontWeight: 500,
                              }}
                              title="View logs"
                            >
                              📄 Logs
                            </button>
                            <button
                              onClick={async () => {
                                if (window.confirm(`Delete run "${run.name}"?\n\nThis action cannot be undone.`)) {
                                  try {
                                    await api.delete(`/runs/${run.id}`);
                                    alert('Run deleted successfully ✅');
                                    fetchRuns();
                                  } catch (err) {
                                    console.error('Failed to delete run:', err);
                                    const error = err as { response?: { data?: { detail?: string } }; message?: string };
                                    alert(`Failed to delete run ❌\n${error.response?.data?.detail || error.message || 'Unknown error'}`);
                                  }
                                }
                              }}
                              style={{
                                padding: '0.4rem 0.8rem',
                                background: '#dc3545',
                                color: 'white',
                                border: 'none',
                                borderRadius: '4px',
                                cursor: 'pointer',
                                fontSize: '0.85rem',
                                fontWeight: 500,
                              }}
                              title="Delete run"
                            >
                              🗑️ Delete
                            </button>
                          </div>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            ) : (
              <div style={{
                textAlign: 'center',
                padding: '3rem',
                backgroundColor: '#f8f9fa',
                borderRadius: '8px',
                border: '2px dashed #dee2e6',
                marginTop: '1rem',
              }}>
                <div style={{ fontSize: '3rem', marginBottom: '1rem' }}>🚀</div>
                <p style={{ fontSize: '1.1rem', fontWeight: 500, marginBottom: '0.5rem', color: '#212529' }}>
                  No runs yet for this project
                </p>
                <p style={{ color: '#6c757d' }}>
                  Click &quot;+ New Run&quot; button above to create your first pipeline run.
                </p>
              </div>
            )}
          </>
        )}
      </div>
    </section>
  );
};

export default RunsPage;
