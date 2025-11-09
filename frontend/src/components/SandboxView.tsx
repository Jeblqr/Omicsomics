import React, { useState, useEffect } from 'react';
import { useProjects } from '../contexts/ProjectsContext';
import api from '../lib/api';

interface DataFile {
  id: number;
  filename: string;
  object_key: string;
  metadata: any;
}

interface Run {
  id: number;
  name: string;
  status: string;
  created_at: string;
}

export const SandboxView: React.FC = () => {
  const { currentProject } = useProjects();
  const [files, setFiles] = useState<DataFile[]>([]);
  const [runs, setRuns] = useState<Run[]>([]);
  const [uploadFile, setUploadFile] = useState<File | null>(null);
  const [runName, setRunName] = useState('');
  const [isUploading, setIsUploading] = useState(false);
  const [isCreatingRun, setIsCreatingRun] = useState(false);

  useEffect(() => {
    if (currentProject) {
      fetchData();
      fetchRuns();
    }
  }, [currentProject]);

  const fetchData = async () => {
    if (!currentProject) return;
    try {
      const response = await api.get<DataFile[]>(`/data/?project_id=${currentProject.id}`);
      setFiles(response.data);
    } catch (error) {
      console.error('Failed to fetch files:', error);
    }
  };

  const fetchRuns = async () => {
    if (!currentProject) return;
    try {
      const response = await api.get<Run[]>(`/runs/?project_id=${currentProject.id}`);
      setRuns(response.data);
    } catch (error) {
      console.error('Failed to fetch runs:', error);
    }
  };

  const handleUpload = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!uploadFile || !currentProject) return;

    setIsUploading(true);
    try {
      const formData = new FormData();
      formData.append('file', uploadFile);
      formData.append('project_id', String(currentProject.id));

      await api.post('/data/upload', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });

      alert('File uploaded successfully!');
      setUploadFile(null);
      fetchData();
    } catch (error) {
      console.error('Upload failed:', error);
      alert('Upload failed');
    } finally {
      setIsUploading(false);
    }
  };

  const handleCreateRun = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!runName || !currentProject) return;

    setIsCreatingRun(true);
    try {
      await api.post('/runs/', {
        name: runName,
        description: '',
        project_id: currentProject.id,
      });

      alert('Run created!');
      setRunName('');
      fetchRuns();
    } catch (error) {
      console.error('Create run failed:', error);
      alert('Failed to create run');
    } finally {
      setIsCreatingRun(false);
    }
  };

  const handleDownload = async (fileId: number, filename: string) => {
    try {
      const response = await api.get(`/data/${fileId}/download-decrypted`, {
        responseType: 'blob',
      });
      const url = window.URL.createObjectURL(new Blob([response.data]));
      const link = document.createElement('a');
      link.href = url;
      link.setAttribute('download', filename);
      document.body.appendChild(link);
      link.click();
      link.remove();
    } catch (error) {
      console.error('Download failed:', error);
      alert('Download failed');
    }
  };

  if (!currentProject) {
    return <div style={{ padding: '2rem' }}>Please select a project above.</div>;
  }

  return (
    <div style={{ padding: '2rem' }}>
      <h2>Sandbox for {currentProject.name}</h2>

      {/* Upload Section */}
      <section style={{ marginTop: '2rem', padding: '1rem', border: '1px solid #ddd' }}>
        <h3>Upload File</h3>
        <form onSubmit={handleUpload} style={{ display: 'flex', gap: '1rem', alignItems: 'center' }}>
          <input
            type="file"
            onChange={(e) => setUploadFile(e.target.files?.[0] || null)}
            required
          />
          <button type="submit" disabled={isUploading || !uploadFile}>
            {isUploading ? 'Uploading...' : 'Upload'}
          </button>
        </form>
      </section>

      {/* Files List */}
      <section style={{ marginTop: '2rem' }}>
        <h3>Data Files ({files.length})</h3>
        {files.length === 0 ? (
          <p>No files uploaded yet.</p>
        ) : (
          <table style={{ width: '100%', borderCollapse: 'collapse' }}>
            <thead>
              <tr style={{ borderBottom: '2px solid #333' }}>
                <th style={{ padding: '0.5rem', textAlign: 'left' }}>Filename</th>
                <th style={{ padding: '0.5rem', textAlign: 'left' }}>Object Key</th>
                <th style={{ padding: '0.5rem', textAlign: 'left' }}>Actions</th>
              </tr>
            </thead>
            <tbody>
              {files.map((file) => (
                <tr key={file.id} style={{ borderBottom: '1px solid #ddd' }}>
                  <td style={{ padding: '0.5rem' }}>{file.filename}</td>
                  <td style={{ padding: '0.5rem', fontSize: '0.85rem', color: '#666' }}>
                    {file.object_key}
                  </td>
                  <td style={{ padding: '0.5rem' }}>
                    <button onClick={() => handleDownload(file.id, file.filename)} style={{ marginRight: '0.5rem' }}>
                      Download (Decrypted)
                    </button>
                    <button
                      onClick={async () => {
                        if (window.confirm(`Delete file "${file.filename}"?`)) {
                          try {
                            await api.delete(`/data/${file.id}`);
                            fetchData();
                          } catch (err) {
                            console.error('Failed to delete file:', err);
                            alert('Failed to delete file');
                          }
                        }
                      }}
                      style={{ background: '#dc3545', color: 'white', border: 'none', padding: '0.25rem 0.75rem', borderRadius: '4px', cursor: 'pointer' }}
                    >
                      Delete
                    </button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
      </section>

      {/* Runs Section */}
      <section style={{ marginTop: '2rem' }}>
        <h3>Runs</h3>
        <form onSubmit={handleCreateRun} style={{ display: 'flex', gap: '1rem', marginBottom: '1rem' }}>
          <input
            type="text"
            placeholder="Run name"
            value={runName}
            onChange={(e) => setRunName(e.target.value)}
            required
            style={{ padding: '0.5rem', flex: 1 }}
          />
          <button type="submit" disabled={isCreatingRun || !runName}>
            {isCreatingRun ? 'Creating...' : 'Create Run'}
          </button>
        </form>

        {runs.length === 0 ? (
          <p>No runs yet.</p>
        ) : (
          <ul>
            {runs.map((run) => (
              <li key={run.id} style={{ padding: '0.5rem', borderBottom: '1px solid #eee' }}>
                <strong>{run.name}</strong> - Status: {run.status} (Created: {new Date(run.created_at).toLocaleString()})
              </li>
            ))}
          </ul>
        )}
      </section>
    </div>
  );
};
