import { useState } from 'react';
import LoadingView from '../../components/LoadingView';
import { useProjectsContext } from '../../contexts/ProjectsContext';

const ProjectsPage = () => {
  const { projects, isLoading, createProject, deleteProject, updateProject, setCurrentProject } = useProjectsContext();
  const [showForm, setShowForm] = useState(false);
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [submitError, setSubmitError] = useState('');
  const [editingId, setEditingId] = useState<number | null>(null);

  const handleCreateOrUpdate = async (e: React.FormEvent) => {
    e.preventDefault();
    setIsSubmitting(true);
    setSubmitError('');

    try {
      if (editingId) {
        await updateProject(editingId, name, description);
      } else {
        await createProject(name, description);
      }
      
      setName('');
      setDescription('');
      setShowForm(false);
      setEditingId(null);
    } catch (err) {
      const error = err as { response?: { data?: { detail?: string } }; message?: string };
      setSubmitError(error.response?.data?.detail || error.message || 'Failed to save project');
    } finally {
      setIsSubmitting(false);
    }
  };

  const handleEdit = (project: { id: number; name: string; description?: string }) => {
    setEditingId(project.id);
    setName(project.name);
    setDescription(project.description || '');
    setShowForm(true);
  };

  const handleDelete = async (id: number) => {
    if (!confirm('Are you sure you want to delete this project? This action cannot be undone.')) {
      return;
    }

    try {
      await deleteProject(id);
    } catch (err) {
      const error = err as { response?: { data?: { detail?: string } }; message?: string };
      alert(error.response?.data?.detail || error.message || 'Failed to delete project');
    }
  };

  const handleCancel = () => {
    setShowForm(false);
    setEditingId(null);
    setName('');
    setDescription('');
    setSubmitError('');
  };

  return (
    <section>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
        <h2>Projects</h2>
        <button
          onClick={() => {
            if (showForm) {
              handleCancel();
            } else {
              setShowForm(true);
            }
          }}
          style={{
            padding: '0.5rem 1rem',
            background: showForm ? '#6c757d' : '#007bff',
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
          border: '1px solid #dee2e6',
        }}>
          <h3 style={{ marginTop: 0, color: '#212529' }}>
            {editingId ? 'Edit Project' : 'Create New Project'}
          </h3>
          <form onSubmit={handleCreateOrUpdate}>
            <div style={{ marginBottom: '1rem' }}>
              <label htmlFor="project-name" style={{ display: 'block', marginBottom: '0.5rem', color: '#212529', fontWeight: 500 }}>
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
                  color: '#212529',
                  backgroundColor: '#ffffff',
                  fontSize: '1rem',
                }}
              />
            </div>

            <div style={{ marginBottom: '1rem' }}>
              <label htmlFor="project-desc" style={{ display: 'block', marginBottom: '0.5rem', color: '#212529', fontWeight: 500 }}>
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
                  color: '#212529',
                  backgroundColor: '#ffffff',
                  fontSize: '1rem',
                  fontFamily: 'inherit',
                }}
              />
            </div>

            {submitError && (
              <div style={{ color: '#dc3545', marginBottom: '1rem', padding: '0.75rem', backgroundColor: '#f8d7da', borderRadius: '4px', border: '1px solid #f5c6cb' }} role="alert">
                {submitError}
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
              {isSubmitting ? 'Creating...' : 'Create Project'}
            </button>
          </form>
        </div>
      )}

      {isLoading && <LoadingView />}
      {!isLoading && projects.length > 0 ? (
        <table className="data-table">
          <thead>
            <tr>
              <th>Name</th>
              <th>Description</th>
              <th>Last Updated</th>
              <th>Actions</th>
            </tr>
          </thead>
          <tbody>
            {projects.map((project) => (
              <tr key={project.id}>
                <td>
                  <strong>{project.name}</strong>
                </td>
                <td>{project.description || '—'}</td>
                <td>{new Date(project.updated_at).toLocaleString()}</td>
                <td>
                  <button
                    onClick={() => setCurrentProject(project)}
                    style={{
                      padding: '0.25rem 0.75rem',
                      marginRight: '0.5rem',
                      background: '#28a745',
                      color: 'white',
                      border: 'none',
                      borderRadius: '4px',
                      cursor: 'pointer',
                      fontSize: '0.85rem',
                    }}
                  >
                    Select
                  </button>
                  <button
                    onClick={() => handleEdit(project)}
                    style={{
                      padding: '0.25rem 0.75rem',
                      marginRight: '0.5rem',
                      background: '#ffc107',
                      color: '#212529',
                      border: 'none',
                      borderRadius: '4px',
                      cursor: 'pointer',
                      fontSize: '0.85rem',
                    }}
                  >
                    Edit
                  </button>
                  <button
                    onClick={() => handleDelete(project.id)}
                    style={{
                      padding: '0.25rem 0.75rem',
                      background: '#dc3545',
                      color: 'white',
                      border: 'none',
                      borderRadius: '4px',
                      cursor: 'pointer',
                      fontSize: '0.85rem',
                    }}
                  >
                    Delete
                  </button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      ) : (
        !isLoading &&           <p>No projects yet. Click &quot;Create New Project&quot; to get started.</p>
      )}
    </section>
  );
};

export default ProjectsPage;
