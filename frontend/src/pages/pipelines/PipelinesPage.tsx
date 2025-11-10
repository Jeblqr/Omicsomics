import { useState, useEffect } from 'react';
import api from '../../lib/api';
import LoadingView from '../../components/LoadingView';

interface PipelineStep {
  name: string;
  tool: string;
  version: string;
  parameters: Record<string, any>;
}

interface PipelineTemplate {
  id: string;
  name: string;
  description: string;
  category: string;
  steps: PipelineStep[];
  parameters: Record<string, any>;
  inputs: string[];
  outputs: string[];
}

const PipelinesPage = () => {
  const [templates, setTemplates] = useState<PipelineTemplate[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');
  const [selectedCategory, setSelectedCategory] = useState<string>('all');

  useEffect(() => {
    fetchTemplates();
  }, [selectedCategory]);

  const fetchTemplates = async () => {
    setIsLoading(true);
    setError('');
    try {
      const url = selectedCategory === 'all' 
        ? '/pipelines/' 
        : `/pipelines/?category=${selectedCategory}`;
      const response = await api.get<PipelineTemplate[]>(url);
      setTemplates(response.data);
    } catch (err) {
      console.error('Failed to fetch pipeline templates:', err);
      const error = err as { response?: { data?: { detail?: string }; status?: number }; message?: string };
      if (error.response?.status === 401) {
        setError('⚠️ Please log in to view pipeline templates');
      } else {
        setError(error.response?.data?.detail || 'Failed to load pipeline templates');
      }
    } finally {
      setIsLoading(false);
    }
  };

  const categories = [
    { value: 'all', label: 'All Pipelines' },
    { value: 'transcriptomics', label: 'Transcriptomics' },
    { value: 'genomics', label: 'Genomics' },
    { value: 'epigenomics', label: 'Epigenomics' },
    { value: 'singlecell', label: 'Single Cell' },
    { value: 'proteomics', label: 'Proteomics' },
    { value: 'metabolomics', label: 'Metabolomics' },
    { value: 'gwas', label: 'GWAS' },
    { value: 'multiomics', label: 'Multi-omics' },
  ];

  const getCategoryBadgeColor = (category: string): string => {
    const colors: Record<string, string> = {
      transcriptomics: '#17a2b8',
      genomics: '#007bff',
      epigenomics: '#6f42c1',
      singlecell: '#e83e8c',
      proteomics: '#fd7e14',
      metabolomics: '#20c997',
      gwas: '#ffc107',
      multiomics: '#6c757d',
    };
    return colors[category] || '#6c757d';
  };

  return (
    <section>
      <div style={{ padding: '2rem' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
          <h2 style={{ color: '#f3f4f6' }}>Common Pipelines</h2>
        </div>
        
        <p style={{ color: '#9ca3af' }}>Pre-configured analysis pipelines available to all users.</p>

        <div style={{ marginBottom: '1.5rem' }}>
          <label htmlFor="category-filter" style={{ display: 'block', marginBottom: '0.5rem', fontWeight: 500, color: '#f3f4f6' }}>
            Filter by Category:
          </label>
          <select
            id="category-filter"
            value={selectedCategory}
            onChange={(e) => setSelectedCategory(e.target.value)}
            style={{
              padding: '0.5rem',
              borderRadius: '4px',
              border: '1px solid #ced4da',
              fontSize: '1rem',
              minWidth: '200px',
            }}
          >
            {categories.map((cat) => (
              <option key={cat.value} value={cat.value}>
                {cat.label}
              </option>
            ))}
          </select>
        </div>

        {error && (
          <div style={{ color: '#fcd34d', padding: '1rem', backgroundColor: '#422006', border: '1px solid #f59e0b', borderRadius: '4px', marginBottom: '1rem' }}>
            <strong style={{ color: '#fbbf24' }}>{error}</strong>
          </div>
        )}

        {isLoading && <LoadingView />}
        
        {!isLoading && templates.length > 0 && (
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(350px, 1fr))', gap: '1.5rem' }}>
            {templates.map((template) => (
              <div
                key={template.id}
                style={{
                  border: '1px solid #dee2e6',
                  borderRadius: '8px',
                  padding: '1.5rem',
                  backgroundColor: '#ffffff',
                  boxShadow: '0 2px 4px rgba(0,0,0,0.1)',
                }}
              >
                <div style={{ marginBottom: '1rem' }}>
                  <h3 style={{ marginTop: 0, marginBottom: '0.5rem', color: '#212529' }}>
                    {template.name}
                  </h3>
                  <span
                    style={{
                      padding: '0.25rem 0.75rem',
                      borderRadius: '12px',
                      fontSize: '0.85rem',
                      fontWeight: 500,
                      backgroundColor: getCategoryBadgeColor(template.category),
                      color: 'white',
                    }}
                  >
                    {template.category}
                  </span>
                </div>

                <p style={{ color: '#6c757d', marginBottom: '1rem', fontSize: '0.95rem' }}>
                  {template.description}
                </p>

                <div style={{ marginBottom: '1rem' }}>
                  <strong style={{ fontSize: '0.9rem', color: '#495057' }}>Steps:</strong>
                  <div style={{ marginTop: '0.5rem', display: 'flex', flexWrap: 'wrap', gap: '0.5rem' }}>
                    {template.steps.map((step, idx) => (
                      <span
                        key={idx}
                        style={{
                          padding: '0.25rem 0.5rem',
                          backgroundColor: '#e9ecef',
                          borderRadius: '4px',
                          fontSize: '0.85rem',
                          color: '#495057',
                        }}
                        title={`${step.tool} ${step.version}`}
                      >
                        {step.name}
                      </span>
                    ))}
                  </div>
                </div>

                <div style={{ marginBottom: '0.75rem' }}>
                  <strong style={{ fontSize: '0.85rem', color: '#495057' }}>Inputs:</strong>
                  <span style={{ fontSize: '0.85rem', color: '#6c757d', marginLeft: '0.5rem' }}>
                    {template.inputs.join(', ')}
                  </span>
                </div>

                <div>
                  <strong style={{ fontSize: '0.85rem', color: '#495057' }}>Outputs:</strong>
                  <span style={{ fontSize: '0.85rem', color: '#6c757d', marginLeft: '0.5rem' }}>
                    {template.outputs.join(', ')}
                  </span>
                </div>
              </div>
            ))}
          </div>
        )}

        {!isLoading && templates.length === 0 && !error && (
          <p style={{ color: '#9ca3af' }}>No pipeline templates found for the selected category.</p>
        )}
      </div>
    </section>
  );
};

export default PipelinesPage;
