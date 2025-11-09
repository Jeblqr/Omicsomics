/**
 * Create Run from Template Component
 * 
 * Wizard-style interface for creating a new pipeline run from a template:
 * - Select template
 * - Configure parameters
 * - Select input files
 * - Review and create
 */

import React, { useState } from 'react';
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';

interface PipelineTemplate {
  id: string;
  name: string;
  description: string;
  category: string;
  parameters: Record<string, {
    type: string;
    required?: boolean;
    default?: unknown;
    description: string;
  }>;
  inputs: string[];
}

interface DataFile {
  id: number;
  filename: string;
  file_type: string;
  omics_type: string;
}

interface CreateRunDialogProps {
  projectId: number;
  onClose: () => void;
  onSuccess?: () => void;
}

const CreateRunDialog: React.FC<CreateRunDialogProps> = ({
  projectId,
  onClose,
  onSuccess,
}) => {
  const [step, setStep] = useState(1);
  const [selectedTemplate, setSelectedTemplate] = useState<PipelineTemplate | null>(null);
  const [runName, setRunName] = useState('');
  const [runDescription, setRunDescription] = useState('');
  const [parameters, setParameters] = useState<Record<string, unknown>>({});
  const [selectedFiles, setSelectedFiles] = useState<number[]>([]);
  const queryClient = useQueryClient();

  // Fetch templates
  const { data: templates } = useQuery<PipelineTemplate[]>({
    queryKey: ['pipeline-templates'],
    queryFn: async () => {
      const response = await fetch('/api/v1/pipelines/', {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });
      if (!response.ok) throw new Error('Failed to fetch templates');
      return response.json();
    },
  });

  // Fetch data files
  const { data: dataFiles } = useQuery<DataFile[]>({
    queryKey: ['datafiles', projectId],
    queryFn: async () => {
      const response = await fetch(`/api/v1/data/?project_id=${projectId}`, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });
      if (!response.ok) throw new Error('Failed to fetch data files');
      return response.json();
    },
  });

  // Create run mutation
  const createRunMutation = useMutation({
    mutationFn: async (runData: {
      name: string;
      description: string;
      project_id: number;
      pipeline_type: string;
      pipeline_template_id: string;
      parameters: Record<string, unknown>;
      input_files: number[];
      auto_start: boolean;
    }) => {
      const response = await fetch('/api/v1/runs/', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
        body: JSON.stringify(runData),
      });
      if (!response.ok) throw new Error('Failed to create run');
      return response.json();
    },
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['runs', projectId] });
      if (onSuccess) onSuccess();
      onClose();
    },
  });

  const handleTemplateSelect = (template: PipelineTemplate) => {
    setSelectedTemplate(template);
    setRunName(`${template.name} - ${new Date().toLocaleDateString()}`);
    // Set default parameters
    const defaults: Record<string, unknown> = {};
    Object.entries(template.parameters).forEach(([key, param]) => {
      if (param.default !== undefined) {
        defaults[key] = param.default;
      }
    });
    setParameters(defaults);
    setStep(2);
  };

  const handleParameterChange = (key: string, value: unknown) => {
    setParameters((prev) => ({ ...prev, [key]: value }));
  };

  const handleFileToggle = (fileId: number) => {
    setSelectedFiles((prev) =>
      prev.includes(fileId)
        ? prev.filter((id) => id !== fileId)
        : [...prev, fileId]
    );
  };

  const handleCreate = () => {
    if (!selectedTemplate) return;

    createRunMutation.mutate({
      name: runName,
      description: runDescription,
      project_id: projectId,
      pipeline_type: 'template',
      pipeline_template_id: selectedTemplate.id,
      parameters,
      input_files: selectedFiles,
      auto_start: false,
    });
  };

  const isValidStep2 = runName.trim().length > 0;
  const isValidStep3 = Object.entries(selectedTemplate?.parameters || {})
    .filter(([, param]) => param.required)
    .every(([key]) => parameters[key] !== undefined && parameters[key] !== '');
  const isValidStep4 = selectedFiles.length > 0;

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50 p-4">
      <div className="bg-white rounded-lg shadow-xl max-w-4xl w-full max-h-[90vh] flex flex-col">
        {/* Header */}
        <div className="sticky top-0 bg-white border-b px-6 py-4">
          <div className="flex justify-between items-center">
            <h2 className="text-2xl font-bold text-gray-900">Create New Run</h2>
            <button
              onClick={onClose}
              className="text-gray-500 hover:text-gray-700 text-2xl"
            >
              ×
            </button>
          </div>
          
          {/* Progress Steps */}
          <div className="flex items-center mt-4 space-x-2">
            {[1, 2, 3, 4].map((i) => (
              <React.Fragment key={i}>
                <div
                  className={`flex items-center justify-center w-8 h-8 rounded-full ${
                    step >= i ? 'bg-blue-600 text-white' : 'bg-gray-200 text-gray-600'
                  }`}
                >
                  {i}
                </div>
                {i < 4 && (
                  <div
                    className={`flex-1 h-1 ${
                      step > i ? 'bg-blue-600' : 'bg-gray-200'
                    }`}
                  ></div>
                )}
              </React.Fragment>
            ))}
          </div>
          <div className="flex justify-between mt-2 text-sm text-gray-600">
            <span>Template</span>
            <span>Details</span>
            <span>Parameters</span>
            <span>Files</span>
          </div>
        </div>

        {/* Content */}
        <div className="flex-1 overflow-y-auto p-6">
          {/* Step 1: Select Template */}
          {step === 1 && (
            <div className="space-y-4">
              <h3 className="text-lg font-semibold text-gray-900 mb-4">
                Select Pipeline Template
              </h3>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {templates?.map((template) => (
                  <div
                    key={template.id}
                    onClick={() => handleTemplateSelect(template)}
                    className="border border-gray-200 rounded-lg p-4 hover:border-blue-500 hover:shadow-md cursor-pointer transition-all"
                  >
                    <h4 className="font-semibold text-gray-900 mb-2">{template.name}</h4>
                    <p className="text-sm text-gray-600 mb-2">{template.description}</p>
                    <span className="inline-block px-2 py-1 bg-blue-100 text-blue-800 text-xs rounded">
                      {template.category}
                    </span>
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Step 2: Run Details */}
          {step === 2 && (
            <div className="space-y-4">
              <h3 className="text-lg font-semibold text-gray-900 mb-4">Run Details</h3>
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Run Name *
                </label>
                <input
                  type="text"
                  value={runName}
                  onChange={(e) => setRunName(e.target.value)}
                  className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                  placeholder="Enter run name"
                />
              </div>
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Description
                </label>
                <textarea
                  value={runDescription}
                  onChange={(e) => setRunDescription(e.target.value)}
                  className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                  rows={4}
                  placeholder="Enter run description (optional)"
                />
              </div>
              <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
                <h4 className="font-medium text-blue-900 mb-2">Selected Template:</h4>
                <p className="text-blue-800">{selectedTemplate?.name}</p>
              </div>
            </div>
          )}

          {/* Step 3: Configure Parameters */}
          {step === 3 && selectedTemplate && (
            <div className="space-y-4">
              <h3 className="text-lg font-semibold text-gray-900 mb-4">
                Configure Parameters
              </h3>
              {Object.entries(selectedTemplate.parameters).map(([key, param]) => (
                <div key={key}>
                  <label className="block text-sm font-medium text-gray-700 mb-2">
                    {key}
                    {param.required && <span className="text-red-500 ml-1">*</span>}
                  </label>
                  <p className="text-sm text-gray-500 mb-2">{param.description}</p>
                  {param.type === 'boolean' ? (
                    <input
                      type="checkbox"
                      checked={Boolean(parameters[key])}
                      onChange={(e) => handleParameterChange(key, e.target.checked)}
                      className="h-4 w-4 text-blue-600 border-gray-300 rounded"
                    />
                  ) : param.type === 'integer' ? (
                    <input
                      type="number"
                      value={parameters[key] as number || ''}
                      onChange={(e) => handleParameterChange(key, parseInt(e.target.value))}
                      className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                    />
                  ) : param.type === 'float' ? (
                    <input
                      type="number"
                      step="0.01"
                      value={parameters[key] as number || ''}
                      onChange={(e) => handleParameterChange(key, parseFloat(e.target.value))}
                      className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                    />
                  ) : (
                    <input
                      type="text"
                      value={parameters[key] as string || ''}
                      onChange={(e) => handleParameterChange(key, e.target.value)}
                      className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                      placeholder={param.default !== undefined ? `Default: ${param.default}` : ''}
                    />
                  )}
                </div>
              ))}
            </div>
          )}

          {/* Step 4: Select Input Files */}
          {step === 4 && (
            <div className="space-y-4">
              <h3 className="text-lg font-semibold text-gray-900 mb-4">
                Select Input Files
              </h3>
              <div className="space-y-2">
                {dataFiles?.map((file) => (
                  <label
                    key={file.id}
                    className="flex items-center p-4 border border-gray-200 rounded-lg hover:bg-gray-50 cursor-pointer"
                  >
                    <input
                      type="checkbox"
                      checked={selectedFiles.includes(file.id)}
                      onChange={() => handleFileToggle(file.id)}
                      className="h-4 w-4 text-blue-600 border-gray-300 rounded mr-3"
                    />
                    <div className="flex-1">
                      <div className="font-medium text-gray-900">{file.filename}</div>
                      <div className="text-sm text-gray-500">
                        {file.omics_type} • {file.file_type}
                      </div>
                    </div>
                  </label>
                ))}
              </div>
              {dataFiles && dataFiles.length === 0 && (
                <div className="text-center py-8 text-gray-500">
                  No data files available in this project
                </div>
              )}
            </div>
          )}
        </div>

        {/* Footer */}
        <div className="sticky bottom-0 bg-gray-50 border-t px-6 py-4 flex justify-between">
          <button
            onClick={() => step > 1 && setStep(step - 1)}
            disabled={step === 1}
            className="px-4 py-2 border border-gray-300 rounded-lg hover:bg-gray-100 transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
          >
            ← Back
          </button>
          <div className="flex gap-2">
            <button
              onClick={onClose}
              className="px-4 py-2 border border-gray-300 rounded-lg hover:bg-gray-100 transition-colors"
            >
              Cancel
            </button>
            {step < 4 ? (
              <button
                onClick={() => setStep(step + 1)}
                disabled={
                  (step === 2 && !isValidStep2) ||
                  (step === 3 && !isValidStep3)
                }
                className="px-6 py-2 bg-blue-600 hover:bg-blue-700 text-white rounded-lg font-medium transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
              >
                Next →
              </button>
            ) : (
              <button
                onClick={handleCreate}
                disabled={!isValidStep4 || createRunMutation.isPending}
                className="px-6 py-2 bg-green-600 hover:bg-green-700 text-white rounded-lg font-medium transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
              >
                {createRunMutation.isPending ? 'Creating...' : '✓ Create Run'}
              </button>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default CreateRunDialog;
