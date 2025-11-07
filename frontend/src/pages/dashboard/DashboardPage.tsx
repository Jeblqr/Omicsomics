import { useProjects } from '../../hooks/useProjects';
import { useAuth } from '../../contexts/AuthContext';
import LoadingView from '../../components/LoadingView';

const DashboardPage = () => {
  const { user } = useAuth();
  const { data: projects, isLoading } = useProjects();

  if (isLoading) {
    return <LoadingView />;
  }

  return (
    <section>
      <h2>Welcome back, {user?.full_name || user?.email}!</h2>
      <p>Overview of your omics analysis platform.</p>

      <div style={{
        display: 'grid',
        gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))',
        gap: '1.5rem',
        marginTop: '2rem',
      }}>
        <div style={{
          background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
          color: 'white',
          padding: '2rem',
          borderRadius: '12px',
          boxShadow: '0 4px 6px rgba(0,0,0,0.1)',
        }}>
          <h3 style={{ margin: '0 0 0.5rem 0', fontSize: '1rem', opacity: 0.9 }}>
            Total Projects
          </h3>
          <p style={{ fontSize: '2.5rem', fontWeight: 'bold', margin: 0 }}>
            {projects?.length || 0}
          </p>
        </div>

        <div style={{
          background: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
          color: 'white',
          padding: '2rem',
          borderRadius: '12px',
          boxShadow: '0 4px 6px rgba(0,0,0,0.1)',
        }}>
          <h3 style={{ margin: '0 0 0.5rem 0', fontSize: '1rem', opacity: 0.9 }}>
            Active Analyses
          </h3>
          <p style={{ fontSize: '2.5rem', fontWeight: 'bold', margin: 0 }}>
            0
          </p>
        </div>

        <div style={{
          background: 'linear-gradient(135deg, #4facfe 0%, #00f2fe 100%)',
          color: 'white',
          padding: '2rem',
          borderRadius: '12px',
          boxShadow: '0 4px 6px rgba(0,0,0,0.1)',
        }}>
          <h3 style={{ margin: '0 0 0.5rem 0', fontSize: '1rem', opacity: 0.9 }}>
            Data Files
          </h3>
          <p style={{ fontSize: '2.5rem', fontWeight: 'bold', margin: 0 }}>
            0
          </p>
        </div>

        <div style={{
          background: 'linear-gradient(135deg, #43e97b 0%, #38f9d7 100%)',
          color: 'white',
          padding: '2rem',
          borderRadius: '12px',
          boxShadow: '0 4px 6px rgba(0,0,0,0.1)',
        }}>
          <h3 style={{ margin: '0 0 0.5rem 0', fontSize: '1rem', opacity: 0.9 }}>
            Completed Runs
          </h3>
          <p style={{ fontSize: '2.5rem', fontWeight: 'bold', margin: 0 }}>
            0
          </p>
        </div>
      </div>

      <div style={{ marginTop: '3rem' }}>
        <h3>Recent Projects</h3>
        {projects && projects.length > 0 ? (
          <div style={{ marginTop: '1rem' }}>
            {projects.slice(0, 5).map((project) => (
              <div
                key={project.id}
                style={{
                  background: '#f8f9fa',
                  padding: '1rem',
                  marginBottom: '0.5rem',
                  borderRadius: '8px',
                  borderLeft: '4px solid #007bff',
                }}
              >
                <div style={{ fontWeight: 'bold' }}>{project.name}</div>
                <div style={{ fontSize: '0.9rem', color: '#6c757d', marginTop: '0.25rem' }}>
                  {project.description || 'No description'}
                </div>
                <div style={{ fontSize: '0.8rem', color: '#6c757d', marginTop: '0.5rem' }}>
                  Updated: {new Date(project.updated_at).toLocaleDateString()}
                </div>
              </div>
            ))}
          </div>
        ) : (
          <p>No projects yet. Create your first project to get started!</p>
        )}
      </div>

      <div style={{ marginTop: '3rem' }}>
        <h3>Quick Actions</h3>
        <div style={{ display: 'flex', gap: '1rem', marginTop: '1rem', flexWrap: 'wrap' }}>
          <button style={{
            padding: '0.75rem 1.5rem',
            background: '#007bff',
            color: 'white',
            border: 'none',
            borderRadius: '8px',
            cursor: 'pointer',
            fontSize: '1rem',
          }}>
            📊 New Analysis
          </button>
          <button style={{
            padding: '0.75rem 1.5rem',
            background: '#28a745',
            color: 'white',
            border: 'none',
            borderRadius: '8px',
            cursor: 'pointer',
            fontSize: '1rem',
          }}>
            📁 Upload Data
          </button>
          <button style={{
            padding: '0.75rem 1.5rem',
            background: '#17a2b8',
            color: 'white',
            border: 'none',
            borderRadius: '8px',
            cursor: 'pointer',
            fontSize: '1rem',
          }}>
            📈 View Reports
          </button>
        </div>
      </div>
    </section>
  );
};

export default DashboardPage;
