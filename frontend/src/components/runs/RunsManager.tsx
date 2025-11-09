/**
 * Pipeline Runs Management Component
 * 
 * Displays and manages pipeline execution runs:
 * - List all runs with status
 * - Create new runs from templates
 * - Monitor run progress
 * - View logs and results
 * - Start/stop/delete runs
 */

import React, { useState } from 'react';
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';

interface Run {
  id: number;
  name: string;
  description: string;
  status: string;
  progress: number;
  pipeline_type: string | null;
  pipeline_template_id: string | null;
  input_files: number[];
  output_files: number[];
  started_at: string | null;
  finished_at: string | null;
  created_at: string;
  project_id: number;
}

interface RunsManagerProps {
  projectId: number;
}

const RunsManager: React.FC<RunsManagerProps> = ({ projectId }) => {
  const [selectedRun, setSelectedRun] = useState<Run | null>(null);
  const [logsDialogOpen, setLogsDialogOpen] = useState(false);
  const queryClient = useQueryClient();

  // Fetch runs for the project
  const { data: runs, isLoading } = useQuery<Run[]>({
    queryKey: ['runs', projectId],
    queryFn: async () => {
      const response = await fetch(`/api/v1/runs/?project_id=${projectId}`, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });
      if (!response.ok) throw new Error('Failed to fetch runs');
      return response.json();
    },
    refetchInterval: 5000, // Auto-refresh every 5 seconds
  });

  // Start run mutation
  const startRunMutation = useMutation({
    mutationFn: async (runId: number) => {
      const response = await fetch(`/api/v1/runs/${runId}/start`, {
        method: 'POST',
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });
      if (!response.ok) throw new Error('Failed to start run');
      return response.json();
    },
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['runs', projectId] });
    },
  });

  // Stop run mutation
  const stopRunMutation = useMutation({
    mutationFn: async (runId: number) => {
      const response = await fetch(`/api/v1/runs/${runId}/stop`, {
        method: 'POST',
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });
      if (!response.ok) throw new Error('Failed to stop run');
      return response.json();
    },
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['runs', projectId] });
    },
  });

  // Delete run mutation
  const deleteRunMutation = useMutation({
    mutationFn: async (runId: number) => {
      const response = await fetch(`/api/v1/runs/${runId}`, {
        method: 'DELETE',
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });
      if (!response.ok) throw new Error('Failed to delete run');
    },
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['runs', projectId] });
    },
  });

  const getStatusColor = (status: string): string => {
    const colors: Record<string, string> = {
      pending: 'bg-gray-200 text-gray-800',
      running: 'bg-blue-500 text-white animate-pulse',
      completed: 'bg-green-500 text-white',
      failed: 'bg-red-500 text-white',
      cancelled: 'bg-yellow-500 text-white',
    };
    return colors[status] || 'bg-gray-200 text-gray-800';
  };

  const getStatusIcon = (status: string): string => {
    const icons: Record<string, string> = {
      pending: '⏳',
      running: '▶️',
      completed: '✅',
      failed: '❌',
      cancelled: '⏹️',
    };
    return icons[status] || '●';
  };

  const formatDuration = (startTime: string | null, endTime: string | null): string => {
    if (!startTime) return '-';
    const start = new Date(startTime).getTime();
    const end = endTime ? new Date(endTime).getTime() : Date.now();
    const duration = Math.floor((end - start) / 1000);
    
    if (duration < 60) return `${duration}s`;
    if (duration < 3600) return `${Math.floor(duration / 60)}m ${duration % 60}s`;
    return `${Math.floor(duration / 3600)}h ${Math.floor((duration % 3600) / 60)}m`;
  };

  const handleViewLogs = async (run: Run) => {
    setSelectedRun(run);
    setLogsDialogOpen(true);
  };

  if (isLoading) {
    return (
      <div className="flex justify-center items-center min-h-[400px]">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-500"></div>
      </div>
    );
  }

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex justify-between items-center">
        <div>
          <h2 className="text-2xl font-bold text-gray-900">Pipeline Runs</h2>
          <p className="text-gray-600">Monitor and manage pipeline executions</p>
        </div>
        <button className="px-4 py-2 bg-blue-600 hover:bg-blue-700 text-white rounded-lg font-medium transition-colors">
          + New Run
        </button>
      </div>

      {/* Runs List */}
      <div className="bg-white rounded-lg shadow overflow-hidden">
        <table className="min-w-full divide-y divide-gray-200">
          <thead className="bg-gray-50">
            <tr>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Name
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Pipeline
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Status
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Progress
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Duration
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Actions
              </th>
            </tr>
          </thead>
          <tbody className="bg-white divide-y divide-gray-200">
            {runs && runs.length > 0 ? (
              runs.map((run) => (
                <tr key={run.id} className="hover:bg-gray-50">
                  <td className="px-6 py-4 whitespace-nowrap">
                    <div className="text-sm font-medium text-gray-900">{run.name}</div>
                    <div className="text-sm text-gray-500">{run.description}</div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500">
                    {run.pipeline_template_id || 'Custom Pipeline'}
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap">
                    <span className={`px-3 py-1 inline-flex text-xs leading-5 font-semibold rounded-full ${getStatusColor(run.status)}`}>
                      {getStatusIcon(run.status)} {run.status}
                    </span>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap">
                    <div className="flex items-center">
                      <div className="w-full bg-gray-200 rounded-full h-2 mr-2">
                        <div
                          className="bg-blue-600 h-2 rounded-full transition-all duration-300"
                          style={{ width: `${run.progress}%` }}
                        ></div>
                      </div>
                      <span className="text-sm text-gray-700">{run.progress.toFixed(0)}%</span>
                    </div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500">
                    {formatDuration(run.started_at, run.finished_at)}
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm font-medium space-x-2">
                    {run.status === 'pending' && (
                      <button
                        onClick={() => startRunMutation.mutate(run.id)}
                        className="text-green-600 hover:text-green-900"
                        disabled={startRunMutation.isPending}
                      >
                        ▶️ Start
                      </button>
                    )}
                    {run.status === 'running' && (
                      <button
                        onClick={() => stopRunMutation.mutate(run.id)}
                        className="text-red-600 hover:text-red-900"
                        disabled={stopRunMutation.isPending}
                      >
                        ⏹️ Stop
                      </button>
                    )}
                    <button
                      onClick={() => handleViewLogs(run)}
                      className="text-blue-600 hover:text-blue-900"
                    >
                      📋 Logs
                    </button>
                    {(run.status === 'completed' || run.status === 'failed' || run.status === 'cancelled') && (
                      <button
                        onClick={() => {
                          if (confirm(`Delete run "${run.name}"?`)) {
                            deleteRunMutation.mutate(run.id);
                          }
                        }}
                        className="text-red-600 hover:text-red-900"
                        disabled={deleteRunMutation.isPending}
                      >
                        🗑️ Delete
                      </button>
                    )}
                  </td>
                </tr>
              ))
            ) : (
              <tr>
                <td colSpan={6} className="px-6 py-12 text-center text-gray-500">
                  <div className="text-lg">No runs yet</div>
                  <div className="text-sm mt-2">Create your first pipeline run to get started</div>
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>

      {/* Logs Modal */}
      {logsDialogOpen && selectedRun && (
        <LogsDialog
          run={selectedRun}
          onClose={() => setLogsDialogOpen(false)}
        />
      )}
    </div>
  );
};

// Logs Dialog Component
interface LogsDialogProps {
  run: Run;
  onClose: () => void;
}

const LogsDialog: React.FC<LogsDialogProps> = ({ run, onClose }) => {
  const { data: logsData } = useQuery({
    queryKey: ['run-logs', run.id],
    queryFn: async () => {
      const response = await fetch(`/api/v1/runs/${run.id}/logs`, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });
      if (!response.ok) throw new Error('Failed to fetch logs');
      return response.json();
    },
    refetchInterval: run.status === 'running' ? 2000 : false, // Auto-refresh while running
  });

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50 p-4">
      <div className="bg-white rounded-lg shadow-xl max-w-4xl w-full max-h-[90vh] flex flex-col">
        {/* Modal Header */}
        <div className="sticky top-0 bg-white border-b px-6 py-4 flex justify-between items-center">
          <h3 className="text-xl font-bold text-gray-900">Run Logs: {run.name}</h3>
          <button
            onClick={onClose}
            className="text-gray-500 hover:text-gray-700 text-2xl"
          >
            ×
          </button>
        </div>

        {/* Logs Content */}
        <div className="flex-1 overflow-y-auto p-6">
          <pre className="bg-gray-900 text-green-400 p-4 rounded-lg text-sm font-mono overflow-x-auto">
            {logsData?.logs || 'No logs available yet...'}
          </pre>
          
          {logsData?.error_message && (
            <div className="mt-4 bg-red-50 border border-red-200 rounded-lg p-4">
              <h4 className="text-red-800 font-semibold mb-2">Error:</h4>
              <pre className="text-red-700 text-sm whitespace-pre-wrap">
                {logsData.error_message}
              </pre>
            </div>
          )}
        </div>

        {/* Modal Footer */}
        <div className="sticky bottom-0 bg-gray-50 border-t px-6 py-4 flex justify-end">
          <button
            onClick={onClose}
            className="px-4 py-2 border border-gray-300 rounded-lg hover:bg-gray-100 transition-colors"
          >
            Close
          </button>
        </div>
      </div>
    </div>
  );
};

export default RunsManager;
