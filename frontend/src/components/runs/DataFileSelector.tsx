import React, { useEffect, useState } from 'react';
import api from '../../lib/api';

interface DataFileSelectorProps {
  projectId: number;
  onSelect: (fileIds: number[]) => void;
  value?: number[];
}

interface DataFile {
  id: number;
  filename: string;
  object_key: string;
  metadata_: Record<string, unknown>;
  size: number;
  checksum: string;
  project_id: number;
  run_id: number | null;
  uploaded_by_id: number;
  created_at: string;
  updated_at: string;
}

export function DataFileSelector({ projectId, onSelect, value = [] }: DataFileSelectorProps) {
  const [files, setFiles] = useState<DataFile[]>([]);
  const [loading, setLoading] = useState(true);
  const [selectedFiles, setSelectedFiles] = useState<Set<number>>(new Set(value));

  useEffect(() => {
    loadDataFiles();
  }, [projectId]);

  useEffect(() => {
    setSelectedFiles(new Set(value));
  }, [value]);

  const loadDataFiles = async () => {
    try {
      setLoading(true);
      const res = await api.get(`/data/?project_id=${projectId}`);
      setFiles(res.data);
    } catch (error) {
      console.error('Failed to load data files:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleToggle = (fileId: number) => {
    const newSelected = new Set(selectedFiles);
    if (newSelected.has(fileId)) {
      newSelected.delete(fileId);
    } else {
      newSelected.add(fileId);
    }
    setSelectedFiles(newSelected);
    onSelect(Array.from(newSelected));
  };

  const handleSelectAll = () => {
    if (selectedFiles.size === files.length) {
      setSelectedFiles(new Set());
      onSelect([]);
    } else {
      const allIds = files.map(f => f.id);
      setSelectedFiles(new Set(allIds));
      onSelect(allIds);
    }
  };

  const formatFileSize = (bytes: number) => {
    if (bytes < 1024) return `${bytes} B`;
    if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
    if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
    return `${(bytes / (1024 * 1024 * 1024)).toFixed(1)} GB`;
  };

  if (loading) {
    return <div className="text-gray-500">Loading data files...</div>;
  }

  if (files.length === 0) {
    return (
      <div className="text-gray-500 bg-gray-50 p-4 rounded">
        No data files found in this project. Please upload data files first.
      </div>
    );
  }

  return (
    <div className="space-y-3">
      <div className="flex justify-between items-center">
        <label className="block text-sm font-medium">
          Select Input Data Files ({selectedFiles.size} selected)
        </label>
        <button
          type="button"
          onClick={handleSelectAll}
          className="text-sm text-blue-600 hover:text-blue-800"
        >
          {selectedFiles.size === files.length ? 'Deselect All' : 'Select All'}
        </button>
      </div>

      <div className="max-h-64 overflow-y-auto border rounded divide-y">
        {files.map((file) => (
          <label
            key={file.id}
            className={`flex items-start p-3 cursor-pointer hover:bg-gray-50 ${
              selectedFiles.has(file.id) ? 'bg-blue-50' : ''
            }`}
          >
            <input
              type="checkbox"
              checked={selectedFiles.has(file.id)}
              onChange={() => handleToggle(file.id)}
              className="mt-1 mr-3"
            />
            <div className="flex-1 min-w-0">
              <div className="font-medium truncate">{file.filename}</div>
              <div className="text-sm text-gray-500 flex gap-3">
                <span>{formatFileSize(file.size)}</span>
                <span>•</span>
                <span>{file.checksum.substring(0, 8)}...</span>
              </div>
              <div className="text-xs text-gray-400">
                Uploaded: {new Date(file.created_at).toLocaleString()}
              </div>
            </div>
          </label>
        ))}
      </div>
    </div>
  );
}
