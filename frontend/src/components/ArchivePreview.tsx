import React, { useState } from 'react';
import { useMutation, useQuery } from '@tanstack/react-query';
import { api } from '../lib/api';
import { Button } from './ui/Button';
import { Spinner } from './ui/Spinner';

interface ArchiveFile {
  name: string;
  path: string;
  size: number;
  mime_type: string;
  is_dir: boolean;
}

interface ArchivePreviewData {
  is_archive: boolean;
  files?: ArchiveFile[];
  total_files?: number;
  total_size?: number;
  archive_filename?: string;
  message?: string;
}

interface ArchivePreviewProps {
  datafileId: number;
  onFileProcessed?: (datafileId: number) => void;
}

/**
 * Component for previewing and processing files from archives.
 * 
 * Features:
 * - Display archive contents in a tree structure
 * - Select individual files for processing
 * - Multi-select support
 * - File metadata display (size, type)
 * - Batch processing
 */
export function ArchivePreview({ datafileId, onFileProcessed }: ArchivePreviewProps) {
  const [selectedFiles, setSelectedFiles] = useState<Set<string>>(new Set());
  const [processingFile, setProcessingFile] = useState<string | null>(null);

  // Query archive preview
  const {
    data: previewData,
    isLoading,
    error,
  } = useQuery<ArchivePreviewData>({
    queryKey: ['archive-preview', datafileId],
    queryFn: async () => {
      const response = await api.get(`/data/${datafileId}/preview`);
      return response.data;
    },
  });

  // Process file from archive mutation
  const processFileMutation = useMutation({
    mutationFn: async (filePath: string) => {
      const formData = new FormData();
      formData.append('file_path', filePath);
      
      const response = await api.post(
        `/data/${datafileId}/process-from-archive`,
        formData
      );
      return response.data;
    },
    onSuccess: (data, filePath) => {
      setProcessingFile(null);
      setSelectedFiles((prev) => {
        const next = new Set(prev);
        next.delete(filePath);
        return next;
      });
      
      if (onFileProcessed) {
        onFileProcessed(data.datafile_id);
      }
    },
    onError: (error: any) => {
      setProcessingFile(null);
      console.error('Failed to process file:', error);
      alert(`Failed to process file: ${error.response?.data?.detail || error.message}`);
    },
  });

  // Toggle file selection
  const toggleFileSelection = (filePath: string) => {
    setSelectedFiles((prev) => {
      const next = new Set(prev);
      if (next.has(filePath)) {
        next.delete(filePath);
      } else {
        next.add(filePath);
      }
      return next;
    });
  };

  // Select all files
  const selectAll = () => {
    if (!previewData?.files) return;
    setSelectedFiles(new Set(previewData.files.map((f) => f.path)));
  };

  // Deselect all files
  const deselectAll = () => {
    setSelectedFiles(new Set());
  };

  // Process selected file
  const processFile = (filePath: string) => {
    setProcessingFile(filePath);
    processFileMutation.mutate(filePath);
  };

  // Format file size
  const formatSize = (bytes: number): string => {
    if (bytes === 0) return '0 B';
    const k = 1024;
    const sizes = ['B', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return `${(bytes / Math.pow(k, i)).toFixed(2)} ${sizes[i]}`;
  };

  // Get icon for mime type
  const getMimeIcon = (mimeType: string): string => {
    if (mimeType.startsWith('text/')) return '📄';
    if (mimeType.startsWith('image/')) return '🖼️';
    if (mimeType.includes('csv')) return '📊';
    if (mimeType.includes('json')) return '📋';
    if (mimeType.includes('xml')) return '📰';
    if (mimeType.includes('vcf')) return '🧬';
    return '📦';
  };

  if (isLoading) {
    return (
      <div className="flex items-center justify-center p-8">
        <Spinner size="lg" />
        <span className="ml-3 text-gray-600">Loading archive contents...</span>
      </div>
    );
  }

  if (error) {
    return (
      <div className="p-4 bg-red-50 border border-red-200 rounded-lg">
        <h3 className="text-red-800 font-semibold mb-2">Failed to load archive</h3>
        <p className="text-red-600 text-sm">
          {(error as any).response?.data?.detail || (error as Error).message}
        </p>
      </div>
    );
  }

  if (!previewData?.is_archive) {
    return (
      <div className="p-4 bg-yellow-50 border border-yellow-200 rounded-lg">
        <p className="text-yellow-800">{previewData?.message || 'Not an archive file'}</p>
      </div>
    );
  }

  const files = previewData.files || [];

  return (
    <div className="space-y-4">
      {/* Archive Info */}
      <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
        <h3 className="text-blue-900 font-semibold mb-2">Archive Contents</h3>
        <div className="text-sm text-blue-800 space-y-1">
          <div>
            <strong>Filename:</strong> {previewData.archive_filename}
          </div>
          <div>
            <strong>Total Files:</strong> {previewData.total_files}
          </div>
          <div>
            <strong>Total Size:</strong> {formatSize(previewData.total_size || 0)}
          </div>
        </div>
      </div>

      {/* Selection Controls */}
      <div className="flex items-center justify-between">
        <div className="text-sm text-gray-600">
          {selectedFiles.size > 0 ? (
            <span>{selectedFiles.size} file(s) selected</span>
          ) : (
            <span>No files selected</span>
          )}
        </div>
        <div className="flex gap-2">
          <Button
            variant="secondary"
            size="sm"
            onClick={selectAll}
            disabled={files.length === 0}
          >
            Select All
          </Button>
          <Button
            variant="secondary"
            size="sm"
            onClick={deselectAll}
            disabled={selectedFiles.size === 0}
          >
            Deselect All
          </Button>
        </div>
      </div>

      {/* File List */}
      <div className="border border-gray-200 rounded-lg overflow-hidden">
        <div className="max-h-96 overflow-y-auto">
          {files.length === 0 ? (
            <div className="p-8 text-center text-gray-500">No files in archive</div>
          ) : (
            <table className="w-full">
              <thead className="bg-gray-50 border-b border-gray-200 sticky top-0">
                <tr>
                  <th className="w-10 p-3"></th>
                  <th className="text-left p-3 text-sm font-medium text-gray-700">Name</th>
                  <th className="text-left p-3 text-sm font-medium text-gray-700">Path</th>
                  <th className="text-right p-3 text-sm font-medium text-gray-700">Size</th>
                  <th className="text-center p-3 text-sm font-medium text-gray-700">Type</th>
                  <th className="text-center p-3 text-sm font-medium text-gray-700">Action</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-200">
                {files.map((file) => {
                  const isSelected = selectedFiles.has(file.path);
                  const isProcessing = processingFile === file.path;

                  return (
                    <tr
                      key={file.path}
                      className={`hover:bg-gray-50 transition-colors ${
                        isSelected ? 'bg-blue-50' : ''
                      }`}
                    >
                      <td className="p-3">
                        <input
                          type="checkbox"
                          checked={isSelected}
                          onChange={() => toggleFileSelection(file.path)}
                          disabled={isProcessing}
                          className="w-4 h-4 text-blue-600 rounded focus:ring-blue-500"
                        />
                      </td>
                      <td className="p-3">
                        <div className="flex items-center gap-2">
                          <span>{getMimeIcon(file.mime_type)}</span>
                          <span className="text-sm font-medium text-gray-900">
                            {file.name}
                          </span>
                        </div>
                      </td>
                      <td className="p-3">
                        <span className="text-sm text-gray-600 font-mono">{file.path}</span>
                      </td>
                      <td className="p-3 text-right">
                        <span className="text-sm text-gray-600">{formatSize(file.size)}</span>
                      </td>
                      <td className="p-3 text-center">
                        <span className="inline-block px-2 py-1 text-xs font-mono bg-gray-100 text-gray-700 rounded">
                          {file.mime_type}
                        </span>
                      </td>
                      <td className="p-3 text-center">
                        <Button
                          variant="primary"
                          size="sm"
                          onClick={() => processFile(file.path)}
                          disabled={isProcessing || processFileMutation.isPending}
                        >
                          {isProcessing ? (
                            <>
                              <Spinner size="sm" className="mr-2" />
                              Processing...
                            </>
                          ) : (
                            'Process'
                          )}
                        </Button>
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          )}
        </div>
      </div>

      {/* Batch Process (Future Enhancement) */}
      {selectedFiles.size > 1 && (
        <div className="p-4 bg-gray-50 border border-gray-200 rounded-lg">
          <p className="text-sm text-gray-600 mb-2">
            Batch processing is not yet implemented. Please process files individually.
          </p>
          <p className="text-xs text-gray-500">
            Future versions will support processing multiple selected files at once.
          </p>
        </div>
      )}
    </div>
  );
}
