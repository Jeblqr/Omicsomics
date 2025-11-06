import LoadingView from '../../components/LoadingView';
import { useProjects } from '../../hooks/useProjects';

const ProjectsPage = () => {
  const { data, isLoading, isError, error } = useProjects();

  return (
    <section>
      <h2>Projects</h2>
      <p>Manage projects, collaborators, and sample manifests.</p>
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
        !isLoading && <p>No projects yet. Create one via the API or upcoming UI.</p>
      )}
    </section>
  );
};

export default ProjectsPage;
