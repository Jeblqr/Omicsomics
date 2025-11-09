import React, { createContext, useContext, useState, useEffect } from 'react';
import api from '../lib/api';

export interface Project {
  id: number;
  name: string;
  description: string;
  created_at: string;
  updated_at: string;
  owner_id: number;
}

interface ProjectsContextType {
  currentProject: Project | null;
  projects: Project[];
  setCurrentProject: (project: Project | null) => void;
  fetchProjects: () => Promise<void>;
  createProject: (name: string, description: string) => Promise<Project>;
  deleteProject: (id: number) => Promise<void>;
  updateProject: (id: number, name: string, description: string) => Promise<Project>;
  isLoading: boolean;
}

const ProjectsContext = createContext<ProjectsContextType | undefined>(undefined);

export const ProjectsProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const [currentProject, setCurrentProjectState] = useState<Project | null>(null);
  const [projects, setProjects] = useState<Project[]>([]);
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    fetchProjects();
  }, []);

  useEffect(() => {
    // 从 localStorage 恢复当前项目
    const savedProjectId = localStorage.getItem('current_project_id');
    if (savedProjectId && projects.length > 0) {
      const id = parseInt(savedProjectId);
      const found = projects.find(p => p.id === id);
      if (found) {
        setCurrentProjectState(found);
      } else {
        // 如果保存的项目不存在，清除并选择第一个
        localStorage.removeItem('current_project_id');
        if (projects.length > 0) {
          setCurrentProject(projects[0]);
        }
      }
    }
  }, [projects]);

  const setCurrentProject = (project: Project | null) => {
    setCurrentProjectState(project);
    if (project) {
      localStorage.setItem('current_project_id', String(project.id));
    } else {
      localStorage.removeItem('current_project_id');
    }
  };

  const fetchProjects = async () => {
    setIsLoading(true);
    try {
      const response = await api.get<Project[]>('/projects/');
      setProjects(response.data);
    } catch (error) {
      console.error('Failed to fetch projects:', error);
      setProjects([]);
    } finally {
      setIsLoading(false);
    }
  };

  const createProject = async (name: string, description: string): Promise<Project> => {
    const response = await api.post<Project>('/projects/', { name, description });
    const newProject = response.data;
    setProjects([...projects, newProject]);
    // 自动切换到新创建的项目
    setCurrentProject(newProject);
    return newProject;
  };

  const deleteProject = async (id: number) => {
    await api.delete(`/projects/${id}`);
    setProjects(projects.filter(p => p.id !== id));
    if (currentProject?.id === id) {
      setCurrentProject(null);
    }
  };

  const updateProject = async (id: number, name: string, description: string): Promise<Project> => {
    const response = await api.put<Project>(`/projects/${id}`, { name, description });
    const updated = response.data;
    setProjects(projects.map(p => p.id === id ? updated : p));
    if (currentProject?.id === id) {
      setCurrentProject(updated);
    }
    return updated;
  };

  return (
    <ProjectsContext.Provider
      value={{
        currentProject,
        projects,
        setCurrentProject,
        fetchProjects,
        createProject,
        deleteProject,
        updateProject,
        isLoading,
      }}
    >
      {children}
    </ProjectsContext.Provider>
  );
};

export const useProjectsContext = () => {
  const context = useContext(ProjectsContext);
  if (!context) {
    throw new Error('useProjectsContext must be used within ProjectsProvider');
  }
  return context;
};
