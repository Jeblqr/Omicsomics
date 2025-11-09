import React, { useEffect, useState } from 'react';
import api from '../../lib/api';

interface PipelineSelectorProps {
  projectId: number;
  onSelect: (selection: PipelineSelection) => void;
  value?: PipelineSelection;
}

export interface PipelineSelection {
  type: 'template' | 'custom' | 'merged';
  templateId?: string;
  customId?: number;
}

interface PipelineTemplate {
  id: string;
  name: string;
  description: string;
  category: string;
}

interface CustomPipeline {
  id: number;
  name: string;
  description: string;
}

export function PipelineSelector({ projectId, onSelect, value }: PipelineSelectorProps) {
  const [templates, setTemplates] = useState<PipelineTemplate[]>([]);
  const [customPipelines, setCustomPipelines] = useState<CustomPipeline[]>([]);
  const [loading, setLoading] = useState(true);
  const [selectedType, setSelectedType] = useState<'template' | 'custom' | 'merged'>(value?.type || 'template');

  useEffect(() => {
    loadPipelines();
  }, [projectId]);

  const loadPipelines = async () => {
    try {
      setLoading(true);
      // Load pipeline templates
      const templatesRes = await api.get('/pipelines/');
      setTemplates(templatesRes.data);

      // Load custom pipelines for this project
      const customRes = await api.get(`/custom-pipelines/?project_id=${projectId}`);
      setCustomPipelines(customRes.data);
    } catch (error) {
      console.error('Failed to load pipelines:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleTypeChange = (type: 'template' | 'custom' | 'merged') => {
    setSelectedType(type);
    onSelect({ type });
  };

  const handleTemplateSelect = (templateId: string) => {
    onSelect({
      type: selectedType,
      templateId,
      customId: value?.customId,
    });
  };

  const handleCustomSelect = (customId: number) => {
    onSelect({
      type: selectedType,
      customId,
      templateId: value?.templateId,
    });
  };

  if (loading) {
    return <div className="text-gray-500">Loading pipelines...</div>;
  }

  return (
    <div className="space-y-4">
      <div>
        <label className="block text-sm font-medium mb-2">Pipeline Type</label>
        <div className="flex gap-4">
          <label className="flex items-center">
            <input
              type="radio"
              name="pipeline-type"
              value="template"
              checked={selectedType === 'template'}
              onChange={() => handleTypeChange('template')}
              className="mr-2"
            />
            Template Pipeline
          </label>
          <label className="flex items-center">
            <input
              type="radio"
              name="pipeline-type"
              value="custom"
              checked={selectedType === 'custom'}
              onChange={() => handleTypeChange('custom')}
              className="mr-2"
            />
            Custom Pipeline
          </label>
          <label className="flex items-center">
            <input
              type="radio"
              name="pipeline-type"
              value="merged"
              checked={selectedType === 'merged'}
              onChange={() => handleTypeChange('merged')}
              className="mr-2"
            />
            Merged Workflow
          </label>
        </div>
      </div>

      {(selectedType === 'template' || selectedType === 'merged') && (
        <div>
          <label className="block text-sm font-medium mb-2">Select Template</label>
          <select
            className="w-full border rounded p-2"
            value={value?.templateId || ''}
            onChange={(e) => handleTemplateSelect(e.target.value)}
          >
            <option value="">Choose a template...</option>
            {templates.map((tmpl) => (
              <option key={tmpl.id} value={tmpl.id}>
                {tmpl.name} - {tmpl.description}
              </option>
            ))}
          </select>
        </div>
      )}

      {(selectedType === 'custom' || selectedType === 'merged') && (
        <div>
          <label className="block text-sm font-medium mb-2">
            {selectedType === 'merged' ? 'Select Additional Custom Pipeline' : 'Select Custom Pipeline'}
          </label>
          <select
            className="w-full border rounded p-2"
            value={value?.customId || ''}
            onChange={(e) => handleCustomSelect(Number(e.target.value))}
          >
            <option value="">Choose a custom pipeline...</option>
            {customPipelines.map((pipe) => (
              <option key={pipe.id} value={pipe.id}>
                {pipe.name} - {pipe.description}
              </option>
            ))}
          </select>
        </div>
      )}

      {selectedType === 'merged' && (
        <div className="text-sm text-gray-600 bg-blue-50 p-3 rounded">
          <strong>Merged Workflow:</strong> Combines a template pipeline with your custom pipeline steps.
          This allows you to extend standard workflows with custom processing.
        </div>
      )}
    </div>
  );
}
