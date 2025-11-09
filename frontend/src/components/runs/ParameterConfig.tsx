import React from 'react';

interface ParameterConfigProps {
  parameters: Record<string, any>;
  onChange: (params: Record<string, any>) => void;
}

export function ParameterConfig({ parameters, onChange }: ParameterConfigProps) {
  const handleChange = (key: string, value: any) => {
    onChange({
      ...parameters,
      [key]: value,
    });
  };

  const handleAddParameter = () => {
    const key = prompt('Enter parameter name:');
    if (key && !parameters[key]) {
      handleChange(key, '');
    }
  };

  const handleRemoveParameter = (key: string) => {
    const newParams = { ...parameters };
    delete newParams[key];
    onChange(newParams);
  };

  return (
    <div className="space-y-3">
      <div className="flex justify-between items-center">
        <label className="block text-sm font-medium">Pipeline Parameters</label>
        <button
          type="button"
          onClick={handleAddParameter}
          className="text-sm text-blue-600 hover:text-blue-800"
        >
          + Add Parameter
        </button>
      </div>

      {Object.keys(parameters).length === 0 ? (
        <div className="text-gray-500 bg-gray-50 p-4 rounded text-sm">
          No parameters configured. Click "Add Parameter" to add custom parameters.
        </div>
      ) : (
        <div className="space-y-2">
          {Object.entries(parameters).map(([key, value]) => (
            <div key={key} className="flex gap-2 items-start">
              <div className="flex-1 grid grid-cols-2 gap-2">
                <input
                  type="text"
                  value={key}
                  disabled
                  className="border rounded p-2 bg-gray-50"
                  placeholder="Parameter name"
                />
                <input
                  type="text"
                  value={String(value)}
                  onChange={(e) => handleChange(key, e.target.value)}
                  className="border rounded p-2"
                  placeholder="Parameter value"
                />
              </div>
              <button
                type="button"
                onClick={() => handleRemoveParameter(key)}
                className="text-red-600 hover:text-red-800 p-2"
                title="Remove parameter"
              >
                ✕
              </button>
            </div>
          ))}
        </div>
      )}

      <div className="text-xs text-gray-500 bg-blue-50 p-2 rounded">
        <strong>Tip:</strong> Parameters will be passed to the pipeline execution.
        Common parameters include: threads, memory, output_dir, quality_threshold, etc.
      </div>
    </div>
  );
}
