import { useQuery } from '@tanstack/react-query';
import api from '../lib/api';

export interface Project {
  id: number;
  name: string;
  description: string | null;
  created_at: string;
  updated_at: string;
}

const fetchProjects = async (): Promise<Project[]> => {
  const response = await api.get<Project[]>('/projects');
  return response.data;
};

export const useProjects = () =>
  useQuery({
    queryKey: ['projects'],
    queryFn: fetchProjects,
  });
