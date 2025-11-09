/**
 * Pipeline Template Library Component
 * 
 * Displays available pipeline templates and allows users to:
 * - Browse templates by category
 * - Search templates
 * - View template details
 * - Create pipelines from templates
 */

import React, { useState, useMemo } from 'react';
import { useQuery } from '@tanstack/react-query';

interface PipelineStep {
  name: string;
  tool: string;
  version: string;
  parameters: Record<string, unknown>;
}

interface PipelineParameter {
  type: string;
  required?: boolean;
  default?: unknown;
  description: string;
}

interface PipelineTemplate {
  id: string;
  name: string;
  description: string;
  category: string;
  steps: PipelineStep[];
  parameters: Record<string, PipelineParameter>;
  inputs: string[];
  outputs: string[];
}

interface PipelineTemplateLibraryProps {
  onCreateFromTemplate?: (templateId: string) => void;
}

const PipelineTemplateLibrary: React.FC<PipelineTemplateLibraryProps> = ({
  onCreateFromTemplate,
}) => {
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedCategory, setSelectedCategory] = useState('all');
  const [selectedTemplate, setSelectedTemplate] = useState<PipelineTemplate | null>(null);
  const [detailsDialogOpen, setDetailsDialogOpen] = useState(false);

  // Fetch templates
  const {
    data: templates,
    isLoading,
    error,
  } = useQuery<PipelineTemplate[]>({
    queryKey: ['pipeline-templates', selectedCategory],
    queryFn: async () => {
      const url = selectedCategory === 'all'
        ? '/api/v1/pipelines/'
        : `/api/v1/pipelines/?category=${selectedCategory}`;
      
      const response = await fetch(url, {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });
      
      if (!response.ok) {
        throw new Error('Failed to fetch templates');
      }
      
      return response.json();
    },
  });

  // Fetch categories
  const { data: categoriesData } = useQuery<{ categories: string[] }>({
    queryKey: ['pipeline-categories'],
    queryFn: async () => {
      const response = await fetch('/api/v1/pipelines/categories', {
        headers: {
          Authorization: `Bearer ${localStorage.getItem('token')}`,
        },
      });
      
      if (!response.ok) {
        throw new Error('Failed to fetch categories');
      }
      
      return response.json();
    },
  });

  // Filter templates by search query
  const filteredTemplates = useMemo(() => {
    if (!templates) return [];
    
    if (!searchQuery) return templates;
    
    const query = searchQuery.toLowerCase();
    return templates.filter(
      (template) =>
        template.name.toLowerCase().includes(query) ||
        template.description.toLowerCase().includes(query)
    );
  }, [templates, searchQuery]);

  const handleViewDetails = (template: PipelineTemplate) => {
    setSelectedTemplate(template);
    setDetailsDialogOpen(true);
  };

  const handleCreateFromTemplate = () => {
    if (selectedTemplate && onCreateFromTemplate) {
      onCreateFromTemplate(selectedTemplate.id);
      setDetailsDialogOpen(false);
    }
  };

  const getCategoryColor = (category: string): string => {
    const colors: Record<string, string> = {
      transcriptomics: '#4caf50',
      genomics: '#2196f3',
      proteomics: '#ff9800',
      metabolomics: '#9c27b0',
      epigenomics: '#f44336',
      singlecell: '#00bcd4',
      multiomics: '#3f51b5',
      gwas: '#e91e63',
      quality_control: '#607d8b',
      enrichment: '#8bc34a',
    };
    return colors[category] || '#757575';
  };

  if (isLoading) {
    return (
      <div className="flex justify-center items-center min-h-[400px]">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-500"></div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="bg-red-50 border border-red-200 rounded-lg p-4 text-red-800">
        <strong>Error:</strong> Failed to load pipeline templates: {(error as Error).message}
      </div>
    );
  }

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="mb-6">
        <h1 className="text-3xl font-bold text-gray-900 mb-2">Pipeline Template Library</h1>
        <p className="text-gray-600">Browse and use pre-configured analysis workflows</p>
      </div>

      {/* Filters */}
      <div className="bg-white rounded-lg shadow p-4 mb-6">
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div>
            <input
              type="text"
              placeholder="🔍 Search templates..."
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
            />
          </div>
          <div>
            <select
              value={selectedCategory}
              onChange={(e) => setSelectedCategory(e.target.value)}
              className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
            >
              <option value="all">All Categories</option>
              {categoriesData?.categories.map((cat) => (
                <option key={cat} value={cat}>
                  {cat.replace(/_/g, ' ').replace(/\b\w/g, (l) => l.toUpperCase())}
                </option>
              ))}
            </select>
          </div>
        </div>
      </div>

      {/* Template Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
        {filteredTemplates.map((template) => (
          <div
            key={template.id}
            className="bg-white rounded-lg shadow hover:shadow-lg transition-all duration-200 hover:-translate-y-1 flex flex-col"
          >
            <div className="p-6 flex-grow">
              <div className="flex justify-between items-start mb-3">
                <h3 className="text-xl font-semibold text-gray-900">{template.name}</h3>
                <span
                  className="px-3 py-1 rounded-full text-xs font-medium text-white"
                  style={{ backgroundColor: getCategoryColor(template.category) }}
                >
                  {template.category}
                </span>
              </div>

              <p className="text-gray-600 text-sm mb-4">{template.description}</p>

              <div className="mb-3">
                <span className="text-xs text-gray-500">Steps: {template.steps.length}</span>
              </div>

              <div className="flex flex-wrap gap-2">
                {template.steps.slice(0, 3).map((step, idx) => (
                  <span
                    key={idx}
                    className="px-2 py-1 bg-gray-100 border border-gray-300 rounded text-xs"
                  >
                    {step.tool}
                  </span>
                ))}
                {template.steps.length > 3 && (
                  <span className="px-2 py-1 bg-gray-100 border border-gray-300 rounded text-xs">
                    +{template.steps.length - 3} more
                  </span>
                )}
              </div>
            </div>

            <div className="px-6 py-4 bg-gray-50 flex justify-between items-center border-t">
              <button
                onClick={() => handleViewDetails(template)}
                className="text-blue-600 hover:text-blue-800 text-sm font-medium"
              >
                ℹ️ Details
              </button>
              <button
                onClick={() => handleViewDetails(template)}
                className="px-4 py-2 bg-blue-600 hover:bg-blue-700 text-white rounded-lg text-sm font-medium transition-colors"
              >
                ▶️ Use Template
              </button>
            </div>
          </div>
        ))}
      </div>

      {filteredTemplates.length === 0 && (
        <div className="text-center py-12">
          <p className="text-xl text-gray-500">No templates found</p>
        </div>
      )}

      {/* Template Details Modal */}
      {detailsDialogOpen && selectedTemplate && (
        <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50 p-4">
          <div className="bg-white rounded-lg shadow-xl max-w-4xl w-full max-h-[90vh] overflow-y-auto">
            {/* Modal Header */}
            <div className="sticky top-0 bg-white border-b px-6 py-4 flex justify-between items-center">
              <div className="flex items-center gap-3">
                <h2 className="text-2xl font-bold text-gray-900">{selectedTemplate.name}</h2>
                <span
                  className="px-3 py-1 rounded-full text-xs font-medium text-white"
                  style={{ backgroundColor: getCategoryColor(selectedTemplate.category) }}
                >
                  {selectedTemplate.category}
                </span>
              </div>
              <button
                onClick={() => setDetailsDialogOpen(false)}
                className="text-gray-500 hover:text-gray-700 text-2xl"
              >
                ×
              </button>
            </div>

            {/* Modal Content */}
            <div className="p-6 space-y-6">
              <p className="text-gray-700">{selectedTemplate.description}</p>

              <hr className="border-gray-200" />

              <div>
                <h3 className="text-lg font-semibold text-gray-900 mb-3">Pipeline Steps</h3>
                <ul className="space-y-2">
                  {selectedTemplate.steps.map((step, idx) => (
                    <li key={idx} className="flex items-start">
                      <span className="text-blue-600 font-medium mr-2">{idx + 1}.</span>
                      <div>
                        <p className="font-medium text-gray-900">{step.name}</p>
                        <p className="text-sm text-gray-600">
                          Tool: {step.tool} (v{step.version})
                        </p>
                      </div>
                    </li>
                  ))}
                </ul>
              </div>

              <hr className="border-gray-200" />

              <div>
                <h3 className="text-lg font-semibold text-gray-900 mb-3">Parameters</h3>
                <ul className="space-y-3">
                  {Object.entries(selectedTemplate.parameters).map(([key, param]) => (
                    <li key={key} className="border-l-4 border-blue-500 pl-4">
                      <div className="flex items-center gap-2 mb-1">
                        <span className="font-medium text-gray-900">{key}</span>
                        {param.required && (
                          <span className="px-2 py-0.5 bg-red-100 text-red-800 text-xs rounded">
                            Required
                          </span>
                        )}
                      </div>
                      <p className="text-sm text-gray-600">
                        {param.description}
                        {param.default !== undefined && ` (Default: ${param.default})`}
                      </p>
                    </li>
                  ))}
                </ul>
              </div>

              <hr className="border-gray-200" />

              <div className="grid grid-cols-2 gap-6">
                <div>
                  <h3 className="text-lg font-semibold text-gray-900 mb-3">Inputs</h3>
                  <ul className="space-y-1">
                    {selectedTemplate.inputs.map((input) => (
                      <li key={input} className="text-sm text-gray-700">
                        • {input}
                      </li>
                    ))}
                  </ul>
                </div>
                <div>
                  <h3 className="text-lg font-semibold text-gray-900 mb-3">Outputs</h3>
                  <ul className="space-y-1">
                    {selectedTemplate.outputs.map((output) => (
                      <li key={output} className="text-sm text-gray-700">
                        • {output}
                      </li>
                    ))}
                  </ul>
                </div>
              </div>
            </div>

            {/* Modal Footer */}
            <div className="sticky bottom-0 bg-gray-50 border-t px-6 py-4 flex justify-end gap-3">
              <button
                onClick={() => setDetailsDialogOpen(false)}
                className="px-4 py-2 border border-gray-300 rounded-lg hover:bg-gray-100 transition-colors"
              >
                Close
              </button>
              <button
                onClick={handleCreateFromTemplate}
                className="px-6 py-2 bg-blue-600 hover:bg-blue-700 text-white rounded-lg font-medium transition-colors"
              >
                ▶️ Create Pipeline from Template
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default PipelineTemplateLibrary;
