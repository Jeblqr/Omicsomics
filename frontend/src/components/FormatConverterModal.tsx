import React, { useState, useEffect } from 'react';
import axios from 'axios';

interface FormatInfo {
  format_id: string;
  extensions: string[];
  mime_type: string;
  can_convert_to: string[];
  can_convert_from: string[];
}

interface ConversionEstimate {
  estimated_seconds: number;
  estimated_minutes: number;
  file_size_mb: number;
}

interface FormatConverterModalProps {
  isOpen: boolean;
  onClose: () => void;
  initialFilePath?: string;
  onConversionComplete?: (targetPath: string) => void;
}

const FormatConverterModal: React.FC<FormatConverterModalProps> = ({
  isOpen,
  onClose,
  initialFilePath,
  onConversionComplete,
}) => {
  const [sourcePath, setSourcePath] = useState(initialFilePath || '');
  const [targetPath, setTargetPath] = useState('');
  const [sourceFormat, setSourceFormat] = useState('');
  const [targetFormat, setTargetFormat] = useState('');
  const [supportedFormats, setSupportedFormats] = useState<FormatInfo[]>([]);
  const [conversionPath, setConversionPath] = useState<string[]>([]);
  const [estimate, setEstimate] = useState<ConversionEstimate | null>(null);
  const [isConverting, setIsConverting] = useState(false);
  const [error, setError] = useState('');
  const [success, setSuccess] = useState('');

  useEffect(() => {
    if (isOpen) {
      loadSupportedFormats();
      if (initialFilePath) {
        detectSourceFormat(initialFilePath);
      }
    }
  }, [isOpen, initialFilePath]);

  useEffect(() => {
    if (sourceFormat && targetFormat) {
      loadConversionPath();
      if (sourcePath) {
        loadEstimate();
      }
    }
  }, [sourceFormat, targetFormat, sourcePath]);

  const loadSupportedFormats = async () => {
    try {
      const response = await axios.get('/api/formats/supported');
      setSupportedFormats(response.data);
    } catch (err) {
      console.error('Failed to load formats:', err);
    }
  };

  const detectSourceFormat = async (filePath: string) => {
    try {
      const response = await axios.post('/api/formats/detect', {
        file_path: filePath,
      });
      setSourceFormat(response.data.format);
    } catch (err) {
      console.error('Failed to detect format:', err);
    }
  };

  const loadConversionPath = async () => {
    try {
      const response = await axios.post('/api/formats/conversion-path', {
        from_format: sourceFormat,
        to_format: targetFormat,
      });
      setConversionPath(response.data.conversion_path);
    } catch (err) {
      setConversionPath([]);
      console.error('Failed to get conversion path:', err);
    }
  };

  const loadEstimate = async () => {
    try {
      const response = await axios.post('/api/formats/estimate', {
        file_path: sourcePath,
        from_format: sourceFormat,
        to_format: targetFormat,
      });
      setEstimate(response.data);
    } catch (err) {
      setEstimate(null);
      console.error('Failed to get estimate:', err);
    }
  };

  const handleConvert = async () => {
    setError('');
    setSuccess('');
    setIsConverting(true);

    try {
      const response = await axios.post('/api/formats/convert', {
        source_path: sourcePath,
        target_path: targetPath,
        from_format: sourceFormat,
        to_format: targetFormat,
      });

      setSuccess(`Conversion started! ID: ${response.data.id}`);
      
      // Poll for completion
      pollConversionStatus(response.data.id);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Conversion failed');
      setIsConverting(false);
    }
  };

  const pollConversionStatus = async (conversionId: number) => {
    const maxAttempts = 60; // 5 minutes max
    let attempts = 0;

    const poll = setInterval(async () => {
      try {
        const response = await axios.get(`/api/formats/conversions/${conversionId}`);
        const status = response.data.status;

        if (status === 'completed') {
          clearInterval(poll);
          setSuccess('Conversion completed successfully!');
          setIsConverting(false);
          if (onConversionComplete) {
            onConversionComplete(targetPath);
          }
        } else if (status === 'failed') {
          clearInterval(poll);
          setError(`Conversion failed: ${response.data.error_message}`);
          setIsConverting(false);
        }

        attempts++;
        if (attempts >= maxAttempts) {
          clearInterval(poll);
          setError('Conversion timeout - check status later');
          setIsConverting(false);
        }
      } catch (err) {
        clearInterval(poll);
        setError('Failed to check conversion status');
        setIsConverting(false);
      }
    }, 5000); // Check every 5 seconds
  };

  const getAvailableTargetFormats = () => {
    if (!sourceFormat) return [];
    const format = supportedFormats.find((f) => f.format_id === sourceFormat);
    return format?.can_convert_to || [];
  };

  const autoGenerateTargetPath = () => {
    if (!sourcePath || !targetFormat) return;
    
    const sourceFormatInfo = supportedFormats.find((f) => f.format_id === sourceFormat);
    const targetFormatInfo = supportedFormats.find((f) => f.format_id === targetFormat);
    
    if (!sourceFormatInfo || !targetFormatInfo) return;
    
    // Remove source extension and add target extension
    let newPath = sourcePath;
    for (const ext of sourceFormatInfo.extensions) {
      if (newPath.endsWith(ext)) {
        newPath = newPath.slice(0, -ext.length);
        break;
      }
    }
    newPath += targetFormatInfo.extensions[0];
    setTargetPath(newPath);
  };

  useEffect(() => {
    autoGenerateTargetPath();
  }, [targetFormat]);

  if (!isOpen) return null;

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
      <div className="bg-gray-800 rounded-lg shadow-xl w-full max-w-2xl max-h-[90vh] overflow-y-auto">
        {/* Header */}
        <div className="flex items-center justify-between p-6 border-b border-gray-700">
          <h2 className="text-2xl font-bold text-white">格式转换器 (Format Converter)</h2>
          <button
            onClick={onClose}
            className="text-gray-400 hover:text-white transition-colors"
          >
            ✕
          </button>
        </div>

        {/* Content */}
        <div className="p-6 space-y-6">
          {/* Source File */}
          <div>
            <label className="block text-sm font-medium text-gray-300 mb-2">
              源文件路径 (Source File Path)
            </label>
            <input
              type="text"
              value={sourcePath}
              onChange={(e) => {
                setSourcePath(e.target.value);
                detectSourceFormat(e.target.value);
              }}
              className="w-full px-3 py-2 bg-gray-700 border border-gray-600 rounded-md text-white"
              placeholder="/path/to/source/file.csv"
            />
          </div>

          {/* Source Format */}
          <div>
            <label className="block text-sm font-medium text-gray-300 mb-2">
              源格式 (Source Format)
            </label>
            <select
              value={sourceFormat}
              onChange={(e) => setSourceFormat(e.target.value)}
              className="w-full px-3 py-2 bg-gray-700 border border-gray-600 rounded-md text-white"
            >
              <option value="">选择格式 (Select Format)</option>
              {supportedFormats.map((format) => (
                <option key={format.format_id} value={format.format_id}>
                  {format.format_id.toUpperCase()} - {format.extensions.join(', ')}
                </option>
              ))}
            </select>
          </div>

          {/* Target Format */}
          <div>
            <label className="block text-sm font-medium text-gray-300 mb-2">
              目标格式 (Target Format)
            </label>
            <select
              value={targetFormat}
              onChange={(e) => setTargetFormat(e.target.value)}
              className="w-full px-3 py-2 bg-gray-700 border border-gray-600 rounded-md text-white"
              disabled={!sourceFormat}
            >
              <option value="">选择格式 (Select Format)</option>
              {getAvailableTargetFormats().map((format) => (
                <option key={format} value={format}>
                  {format.toUpperCase()}
                </option>
              ))}
            </select>
          </div>

          {/* Target File */}
          <div>
            <label className="block text-sm font-medium text-gray-300 mb-2">
              目标文件路径 (Target File Path)
            </label>
            <input
              type="text"
              value={targetPath}
              onChange={(e) => setTargetPath(e.target.value)}
              className="w-full px-3 py-2 bg-gray-700 border border-gray-600 rounded-md text-white"
              placeholder="/path/to/target/file.xlsx"
            />
          </div>

          {/* Conversion Path */}
          {conversionPath.length > 0 && (
            <div className="bg-gray-700 rounded-lg p-4">
              <h3 className="text-sm font-medium text-gray-300 mb-2">
                转换路径 (Conversion Path)
              </h3>
              <div className="flex items-center space-x-2 text-white">
                {conversionPath.map((format, index) => (
                  <React.Fragment key={index}>
                    <span className="px-3 py-1 bg-blue-600 rounded-md text-sm font-medium">
                      {format.toUpperCase()}
                    </span>
                    {index < conversionPath.length - 1 && (
                      <span className="text-gray-400">→</span>
                    )}
                  </React.Fragment>
                ))}
              </div>
              {conversionPath.length > 2 && (
                <p className="text-sm text-yellow-400 mt-2">
                  ⚠️ 需要 {conversionPath.length - 1} 步转换
                </p>
              )}
            </div>
          )}

          {/* Estimate */}
          {estimate && (
            <div className="bg-gray-700 rounded-lg p-4">
              <h3 className="text-sm font-medium text-gray-300 mb-2">
                预估时间 (Estimated Time)
              </h3>
              <div className="grid grid-cols-2 gap-4 text-white">
                <div>
                  <p className="text-sm text-gray-400">文件大小</p>
                  <p className="text-lg font-semibold">
                    {estimate.file_size_mb.toFixed(2)} MB
                  </p>
                </div>
                <div>
                  <p className="text-sm text-gray-400">预计时间</p>
                  <p className="text-lg font-semibold">
                    {estimate.estimated_seconds < 60
                      ? `${estimate.estimated_seconds.toFixed(1)}s`
                      : `${estimate.estimated_minutes.toFixed(1)}min`}
                  </p>
                </div>
              </div>
            </div>
          )}

          {/* Error Message */}
          {error && (
            <div className="bg-red-900/50 border border-red-700 rounded-lg p-4">
              <p className="text-red-200">{error}</p>
            </div>
          )}

          {/* Success Message */}
          {success && (
            <div className="bg-green-900/50 border border-green-700 rounded-lg p-4">
              <p className="text-green-200">{success}</p>
            </div>
          )}
        </div>

        {/* Footer */}
        <div className="flex justify-end space-x-3 p-6 border-t border-gray-700">
          <button
            onClick={onClose}
            className="px-4 py-2 text-gray-300 hover:text-white transition-colors"
            disabled={isConverting}
          >
            取消 (Cancel)
          </button>
          <button
            onClick={handleConvert}
            disabled={
              !sourcePath ||
              !targetPath ||
              !sourceFormat ||
              !targetFormat ||
              isConverting
            }
            className="px-6 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors disabled:bg-gray-600 disabled:cursor-not-allowed"
          >
            {isConverting ? '转换中... (Converting...)' : '开始转换 (Convert)'}
          </button>
        </div>
      </div>
    </div>
  );
};

export default FormatConverterModal;
