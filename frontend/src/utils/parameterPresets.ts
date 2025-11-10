/**
 * Parameter Preset Management
 * Allows users to save and load parameter configurations for tools
 */

export interface ParameterPreset {
  id: string;
  name: string;
  description?: string;
  toolId: string;
  parameters: Record<string, any>;
  createdAt: string;
  updatedAt: string;
}

const STORAGE_KEY = 'pipeline_parameter_presets';

// Get all presets from localStorage
export const getAllPresets = (): ParameterPreset[] => {
  try {
    const data = localStorage.getItem(STORAGE_KEY);
    return data ? JSON.parse(data) : [];
  } catch (error) {
    console.error('Failed to load presets:', error);
    return [];
  }
};

// Get presets for a specific tool
export const getPresetsForTool = (toolId: string): ParameterPreset[] => {
  return getAllPresets().filter((preset) => preset.toolId === toolId);
};

// Save a new preset
export const savePreset = (preset: Omit<ParameterPreset, 'id' | 'createdAt' | 'updatedAt'>): ParameterPreset => {
  const allPresets = getAllPresets();
  const newPreset: ParameterPreset = {
    ...preset,
    id: `preset_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
    createdAt: new Date().toISOString(),
    updatedAt: new Date().toISOString(),
  };
  
  allPresets.push(newPreset);
  localStorage.setItem(STORAGE_KEY, JSON.stringify(allPresets));
  return newPreset;
};

// Update an existing preset
export const updatePreset = (id: string, updates: Partial<ParameterPreset>): ParameterPreset | null => {
  const allPresets = getAllPresets();
  const index = allPresets.findIndex((p) => p.id === id);
  
  if (index === -1) return null;
  
  allPresets[index] = {
    ...allPresets[index],
    ...updates,
    updatedAt: new Date().toISOString(),
  };
  
  localStorage.setItem(STORAGE_KEY, JSON.stringify(allPresets));
  return allPresets[index];
};

// Delete a preset
export const deletePreset = (id: string): boolean => {
  const allPresets = getAllPresets();
  const filtered = allPresets.filter((p) => p.id !== id);
  
  if (filtered.length === allPresets.length) return false;
  
  localStorage.setItem(STORAGE_KEY, JSON.stringify(filtered));
  return true;
};

// Export presets to JSON file
export const exportPresets = (presets: ParameterPreset[]): void => {
  const data = JSON.stringify(presets, null, 2);
  const blob = new Blob([data], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = `parameter-presets-${Date.now()}.json`;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};

// Import presets from JSON file
export const importPresets = (file: File): Promise<ParameterPreset[]> => {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    
    reader.onload = (e) => {
      try {
        const imported = JSON.parse(e.target?.result as string) as ParameterPreset[];
        
        // Validate structure
        if (!Array.isArray(imported)) {
          throw new Error('Invalid preset file format');
        }
        
        // Merge with existing presets (avoid duplicates by ID)
        const existing = getAllPresets();
        const existingIds = new Set(existing.map((p) => p.id));
        const newPresets = imported.filter((p) => !existingIds.has(p.id));
        
        const merged = [...existing, ...newPresets];
        localStorage.setItem(STORAGE_KEY, JSON.stringify(merged));
        
        resolve(newPresets);
      } catch (error) {
        reject(error);
      }
    };
    
    reader.onerror = () => reject(reader.error);
    reader.readAsText(file);
  });
};

// Default presets for common QC tools
export const getDefaultPresets = (): Omit<ParameterPreset, 'id' | 'createdAt' | 'updatedAt'>[] => {
  return [
    {
      name: 'Standard Quality (Q20)',
      description: 'Standard quality trimming with Q20 threshold',
      toolId: 'genomics-trim-galore',
      parameters: {
        quality: 20,
        length: 20,
        stringency: 1,
        error_rate: 0.1,
        clip_r1: 0,
        clip_r2: 0,
        three_prime_clip_r1: 0,
      },
    },
    {
      name: 'High Quality (Q30)',
      description: 'Stringent quality trimming with Q30 threshold',
      toolId: 'genomics-trim-galore',
      parameters: {
        quality: 30,
        length: 30,
        stringency: 2,
        error_rate: 0.05,
        clip_r1: 0,
        clip_r2: 0,
        three_prime_clip_r1: 0,
      },
    },
    {
      name: 'Illumina Nextera',
      description: 'Optimized for Illumina Nextera adapters',
      toolId: 'genomics-trim-galore',
      parameters: {
        quality: 20,
        length: 20,
        stringency: 1,
        error_rate: 0.1,
        clip_r1: 0,
        clip_r2: 0,
        three_prime_clip_r1: 0,
      },
    },
    {
      name: 'Low Input RNA-seq',
      description: 'Gentle trimming for low-input RNA-seq',
      toolId: 'genomics-trim-galore',
      parameters: {
        quality: 15,
        length: 15,
        stringency: 1,
        error_rate: 0.15,
        clip_r1: 0,
        clip_r2: 0,
        three_prime_clip_r1: 0,
      },
    },
    {
      name: 'Standard FastQC',
      description: 'Default FastQC parameters',
      toolId: 'genomics-fastqc',
      parameters: {
        threads: 4,
        kmers: 7,
      },
    },
    {
      name: 'High Throughput FastQC',
      description: 'FastQC with more threads for large datasets',
      toolId: 'genomics-fastqc',
      parameters: {
        threads: 16,
        kmers: 7,
      },
    },
  ];
};

// Initialize default presets if storage is empty
export const initializeDefaultPresets = (): void => {
  const existing = getAllPresets();
  if (existing.length === 0) {
    const defaults = getDefaultPresets();
    defaults.forEach((preset) => savePreset(preset));
  }
};
