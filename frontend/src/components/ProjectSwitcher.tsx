import React from 'react';
import { useProjectsContext } from '../contexts/ProjectsContext';

export const ProjectSwitcher: React.FC = () => {
  const { currentProject, projects, setCurrentProject, isLoading } = useProjectsContext();

  if (isLoading) {
    return (
      <div style={{
        padding: '1rem',
        borderBottom: '1px solid #dee2e6',
        backgroundColor: '#f8f9fa',
      }}>
        <span style={{ color: '#6c757d' }}>Loading projects...</span>
      </div>
    );
  }

  return (
    <div style={{
      padding: '1rem',
      borderBottom: '2px solid #007bff',
      backgroundColor: '#ffffff',
      display: 'flex',
      alignItems: 'center',
      gap: '1rem',
    }}>
      <label style={{
        display: 'flex',
        alignItems: 'center',
        gap: '0.5rem',
        fontSize: '1rem',
        fontWeight: 500,
        color: '#212529',
      }}>
        <span>📁 Current Project:</span>
        <select
          value={currentProject?.id || ''}
          onChange={(e) => {
            const id = parseInt(e.target.value);
            const project = projects.find(p => p.id === id);
            setCurrentProject(project || null);
          }}
          style={{
            padding: '0.5rem 1rem',
            fontSize: '0.95rem',
            borderRadius: '4px',
            border: '1px solid #ced4da',
            backgroundColor: '#ffffff',
            color: '#212529',
            cursor: 'pointer',
            minWidth: '200px',
          }}
        >
          <option value="">-- Select Project --</option>
          {projects.map((proj) => (
            <option key={proj.id} value={proj.id}>
              {proj.name}
            </option>
          ))}
        </select>
      </label>
      
      {currentProject && (
        <div style={{
          flex: 1,
          fontSize: '0.9rem',
          color: '#6c757d',
          fontStyle: 'italic',
        }}>
          {currentProject.description || 'No description'}
        </div>
      )}

      {!currentProject && projects.length > 0 && (
        <div style={{
          flex: 1,
          fontSize: '0.9rem',
          color: '#dc3545',
          fontWeight: 500,
        }}>
          ⚠️ Please select a project to continue
        </div>
      )}
    </div>
  );
};
