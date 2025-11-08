import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import api from '../lib/api';

export enum WorkflowStatus {
  PENDING = 'pending',
  RUNNING = 'running',
  COMPLETED = 'completed',
  FAILED = 'failed',
  CANCELLED = 'cancelled',
}

export interface Workflow {
  id: number;
  name: string;
  workflow_type: string;
  status: WorkflowStatus;
  sample_id: number;
  input_files: Record<string, any>;
  output_files: Record<string, any>;
  parameters: Record<string, any>;
  logs: string;
  error_message: string | null;
  started_at: string | null;
  completed_at: string | null;
  created_at: string;
  updated_at: string;
}

export interface WorkflowCreate {
  name: string;
  workflow_type: string;
  sample_id: number;
  input_files?: Record<string, any>;
  parameters?: Record<string, any>;
}

const fetchWorkflows = async (sampleId?: number, status?: WorkflowStatus): Promise<Workflow[]> => {
  const params = new URLSearchParams();
  if (sampleId) params.append('sample_id', sampleId.toString());
  if (status) params.append('status', status);
  
  const response = await api.get<Workflow[]>(`/workflows?${params.toString()}`);
  return response.data;
};

const createWorkflow = async (workflow: WorkflowCreate): Promise<Workflow> => {
  const response = await api.post<Workflow>('/workflows/', workflow);
  return response.data;
};

const cancelWorkflow = async (workflowId: number): Promise<Workflow> => {
  const response = await api.put<Workflow>(`/workflows/${workflowId}`, {
    status: WorkflowStatus.CANCELLED,
  });
  return response.data;
};

export const useWorkflows = (sampleId?: number, status?: WorkflowStatus) =>
  useQuery({
    queryKey: ['workflows', sampleId, status],
    queryFn: () => fetchWorkflows(sampleId, status),
  });

export const useCreateWorkflow = () => {
  const queryClient = useQueryClient();
  
  return useMutation({
    mutationFn: createWorkflow,
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['workflows'] });
    },
  });
};

export const useCancelWorkflow = () => {
  const queryClient = useQueryClient();
  
  return useMutation({
    mutationFn: cancelWorkflow,
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['workflows'] });
    },
  });
};
