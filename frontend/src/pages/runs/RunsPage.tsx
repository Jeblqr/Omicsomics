import { useState, useEffect } from 'react';
import { ProjectSwitcher } from '../../components/ProjectSwitcher';
import { useProjectsContext } from '../../contexts/ProjectsContext';
import api from '../../lib/api';
import LoadingView from '../../components/LoadingView';

interface Run {
  id: number;
  name: string;
  description: string;
  status: string;
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

    setIsSubmitting(true);
    setError('');

    try {
      await api.post('/runs/', {
        name,
        description,
        project_id: currentProject.id,
      });

      setName('');
      setDescription('');
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
          <h2>Pipeline Runs</h2>
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
        
        <p>Track workflow executions, logs, and outputs.</p>

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
              <div style={{ marginBottom: '1rem' }}>
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

              <div style={{ marginBottom: '1rem' }}>
                <label htmlFor="run-desc" style={{ display: 'block', marginBottom: '0.5rem', color: '#212529', fontWeight: 500 }}>
                  Description
                </label>
                <textarea
                  id="run-desc"
                  value={description}
                  onChange={(e) => setDescription(e.target.value)}
                  rows={3}
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
              <table className="data-table" style={{ marginTop: '1rem' }}>
                <thead>
                  <tr>
                    <th>Run Name</th>
                    <th>Status</th>
                    <th>Description</th>
                    <th>Created</th>
                    <th>Started</th>
                    <th>Finished</th>
                  </tr>
                </thead>
                <tbody>
                  {runs.map((run) => (
                    <tr key={run.id}>
                      <td>
                        <strong>{run.name}</strong>
                      </td>
                      <td>{getStatusBadge(run.status)}</td>
                      <td>{run.description || '—'}</td>
                      <td>{new Date(run.created_at).toLocaleString()}</td>
                      <td>{run.started_at ? new Date(run.started_at).toLocaleString() : '—'}</td>
                      <td>{run.finished_at ? new Date(run.finished_at).toLocaleString() : '—'}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            ) : (
              <p>No runs yet for this project. Click &quot;New Run&quot; to create one.</p>
            )}
          </>
        )}
      </div>
    </section>
  );
};

export default RunsPage;
