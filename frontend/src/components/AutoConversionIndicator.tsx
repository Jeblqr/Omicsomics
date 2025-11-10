import React from 'react';

interface ConversionInfo {
  edge_id?: string;
  source_node: string;
  target_node: string;
  source_format: string;
  target_format: string;
  conversion_path: string[];
  estimated_seconds: number;
  num_steps: number;
}

interface AutoConversionIndicatorProps {
  conversions: ConversionInfo[];
  warnings: string[];
  totalEstimatedTime: number;
  onApplyConversions?: () => void;
  onDismiss?: () => void;
}

const AutoConversionIndicator: React.FC<AutoConversionIndicatorProps> = ({
  conversions,
  warnings,
  totalEstimatedTime,
  onApplyConversions,
  onDismiss,
}) => {
  if (conversions.length === 0 && warnings.length === 0) {
    return null;
  }

  return (
    <div className="bg-amber-900/90 border border-amber-600 rounded-lg p-4 mb-4">
      {/* Header */}
      <div className="flex items-start justify-between mb-3">
        <div className="flex items-center space-x-2">
          <span className="text-2xl">🔄</span>
          <div>
            <h3 className="text-lg font-semibold text-amber-100">
              自动格式转换 (Auto Format Conversion)
            </h3>
            <p className="text-sm text-amber-200">
              检测到 {conversions.length} 个格式不匹配，需要自动转换
            </p>
          </div>
        </div>
        {onDismiss && (
          <button
            onClick={onDismiss}
            className="text-amber-300 hover:text-amber-100 transition-colors"
          >
            ✕
          </button>
        )}
      </div>

      {/* Conversions List */}
      <div className="space-y-3 mb-4">
        {conversions.map((conversion, idx) => (
          <div
            key={idx}
            className="bg-gray-800/50 rounded-md p-3 border border-amber-700/50"
          >
            <div className="flex items-center justify-between mb-2">
              <div className="flex items-center space-x-2">
                <span className="text-amber-400 font-medium">#{idx + 1}</span>
                <span className="text-gray-300 text-sm">
                  {conversion.source_node} → {conversion.target_node}
                </span>
              </div>
              <span className="text-xs text-gray-400">
                ~{conversion.estimated_seconds.toFixed(1)}s
              </span>
            </div>

            {/* Conversion Path */}
            <div className="flex items-center space-x-1 flex-wrap">
              {conversion.conversion_path.map((format, formatIdx) => (
                <React.Fragment key={formatIdx}>
                  <span
                    className={`px-2 py-1 rounded text-xs font-medium ${
                      formatIdx === 0
                        ? 'bg-blue-600 text-white'
                        : formatIdx === conversion.conversion_path.length - 1
                        ? 'bg-green-600 text-white'
                        : 'bg-amber-600 text-white'
                    }`}
                  >
                    {format.toUpperCase()}
                  </span>
                  {formatIdx < conversion.conversion_path.length - 1 && (
                    <span className="text-amber-400 text-xs mx-1">→</span>
                  )}
                </React.Fragment>
              ))}
            </div>

            {/* Multi-step Warning */}
            {conversion.num_steps > 1 && (
              <div className="mt-2 flex items-center space-x-1 text-xs text-amber-300">
                <span>⚠️</span>
                <span>需要 {conversion.num_steps} 步转换</span>
              </div>
            )}
          </div>
        ))}
      </div>

      {/* Warnings */}
      {warnings.length > 0 && (
        <div className="mb-4 space-y-2">
          {warnings.map((warning, idx) => (
            <div
              key={idx}
              className="bg-red-900/30 border border-red-700/50 rounded-md p-2 text-sm text-red-200"
            >
              <span className="mr-2">⚠️</span>
              {warning}
            </div>
          ))}
        </div>
      )}

      {/* Summary */}
      <div className="bg-gray-800/50 rounded-md p-3 mb-4">
        <div className="grid grid-cols-2 gap-4 text-sm">
          <div>
            <p className="text-gray-400 mb-1">转换数量</p>
            <p className="text-xl font-semibold text-amber-100">
              {conversions.length}
            </p>
          </div>
          <div>
            <p className="text-gray-400 mb-1">预计总时间</p>
            <p className="text-xl font-semibold text-amber-100">
              {totalEstimatedTime < 60
                ? `${totalEstimatedTime.toFixed(1)}s`
                : `${(totalEstimatedTime / 60).toFixed(1)}min`}
            </p>
          </div>
        </div>
      </div>

      {/* Actions */}
      <div className="flex items-center space-x-3">
        {onApplyConversions && (
          <button
            onClick={onApplyConversions}
            className="flex-1 px-4 py-2 bg-amber-600 text-white rounded-md hover:bg-amber-700 transition-colors font-medium"
          >
            自动插入转换节点 (Auto Insert Conversion Nodes)
          </button>
        )}
        {onDismiss && (
          <button
            onClick={onDismiss}
            className="px-4 py-2 text-gray-300 hover:text-white transition-colors"
          >
            忽略 (Ignore)
          </button>
        )}
      </div>

      {/* Info */}
      <div className="mt-3 text-xs text-amber-200/70">
        💡 提示：转换节点将自动插入到管道中，并在执行时自动进行格式转换
      </div>
    </div>
  );
};

export default AutoConversionIndicator;
