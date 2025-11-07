import { useState } from 'react';
import LoadingView from '../../components/LoadingView';
import { useProjects } from '../../hooks/useProjects';
import api from '../../lib/api';

const ProjectsPage = () => {
  const { data, isLoading, isError, error, refetch } = useProjects();
  const [showForm, setShowForm] = useState(false);
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [submitError, setSubmitError] = useState('');

  const handleCreateProject = async (e: React.FormEvent) => {
    e.preventDefault();
    setIsSubmitting(true);
    setSubmitError('');

    try {
      await api.post('/projects/', {
        name,
        description,
      });
      
      setName('');
      setDescription('');
      setShowForm(false);
      refetch();
    } catch (err: any) {
      setSubmitError(err.response?.data?.detail || 'Failed to create project');
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <section>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
        <h2>Projects</h2>
        <button
          onClick={() => setShowForm(!showForm)}
          style={{
            padding: '0.5rem 1rem',
            background: '#007bff',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
          }}
        >
          {showForm ? 'Cancel' : '+ New Project'}
        </button>
      </div>
      
      <p>Manage projects, collaborators, and sample manifests.</p>

      {showForm && (
        <div style={{
          background: '#f8f9fa',
          padding: '1.5rem',
          borderRadius: '8px',
          marginBottom: '1.5rem',
        }}>
          <h3 style={{ marginTop: 0 }}>Create New Project</h3>
          <form onSubmit={handleCreateProject}>
            <div style={{ marginBottom: '1rem' }}>
              <label htmlFor="project-name" style={{ display: 'block', marginBottom: '0.5rem' }}>
                Project Name *
              </label>
              <input
                id="project-name"
                type="text"
                value={name}
                onChange={(e) => setName(e.target.value)}
                required
                style={{
                  width: '100%',
                  padding: '0.5rem',
                  borderRadius: '4px',
                  border: '1px solid #ced4da',
                }}
              />
            </div>

            <div style={{ marginBottom: '1rem' }}>
              <label htmlFor="project-desc" style={{ display: 'block', marginBottom: '0.5rem' }}>
                Description
              </label>
              <textarea
                id="project-desc"
                value={description}
                onChange={(e) => setDescription(e.target.value)}
                rows={3}
                style={{
                  width: '100%',
                  padding: '0.5rem',
                  borderRadius: '4px',
                  border: '1px solid #ced4da',
                }}
              />
            </div>

            {submitError && (
              <div style={{ color: '#dc3545', marginBottom: '1rem' }} role="alert">
                {submitError}
              </div>
            )}

            <button
              type="submit"
              disabled={isSubmitting}
              style={{
                padding: '0.5rem 1.5rem',
                background: '#28a745',
                color: 'white',
                border: 'none',
                borderRadius: '4px',
                cursor: isSubmitting ? 'not-allowed' : 'pointer',
              }}
            >
              {isSubmitting ? 'Creating...' : 'Create Project'}
            </button>
          </form>
        </div>
      )}

      {isLoading && <LoadingView />}
      {isError && <p role="alert">Failed to load projects: {error?.message}</p>}
      {data && data.length > 0 ? (
        <table className="data-table">
          <thead>
            <tr>
              <th>Name</th>
              <th>Description</th>
              <th>Last Updated</th>
            </tr>
          </thead>
          <tbody>
            {data.map((project) => (
              <tr key={project.id}>
                <td>{project.name}</td>
                <td>{project.description ?? '—'}</td>
                <td>{new Date(project.updated_at).toLocaleString()}</td>
              </tr>
            ))}
          </tbody>
        </table>
      ) : (
        !isLoading && <p>No projects yet. Click "New Project" to create one.</p>
      )}
    </section>
  );
};

export default ProjectsPage;
