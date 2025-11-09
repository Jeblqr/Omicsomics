import React from 'react';
import { ScatterPlot, ScatterPoint } from './ScatterPlot';

export interface VolcanoPoint {
  label: string;  // Gene/Protein ID
  logFC: number;  // Log2 fold change
  pValue: number;  // P-value (will be converted to -log10)
  significant?: boolean;  // Pre-computed significance
}

export interface VolcanoPlotProps {
  data: VolcanoPoint[];
  width?: number;
  height?: number;
  foldChangeThreshold?: number;  // Fold change threshold (default: 1)
  pValueThreshold?: number;  // P-value threshold (default: 0.05)
  title?: string;
  onPointClick?: (point: VolcanoPoint) => void;
}

/**
 * Volcano plot for differential expression analysis.
 * 
 * X-axis: Log2 fold change
 * Y-axis: -Log10 p-value
 * 
 * Points are colored by significance:
 * - Red: upregulated (logFC > threshold, p < threshold)
 * - Blue: downregulated (logFC < -threshold, p < threshold)
 * - Gray: not significant
 */
export function VolcanoPlot({
  data,
  width = 700,
  height = 600,
  foldChangeThreshold = 1,
  pValueThreshold = 0.05,
  title = 'Volcano Plot',
  onPointClick,
}: VolcanoPlotProps) {
  // Convert to scatter plot format
  const scatterData: ScatterPoint[] = data.map(point => {
    const negLogP = -Math.log10(point.pValue);
    const isUpregulated = point.logFC >= foldChangeThreshold && point.pValue < pValueThreshold;
    const isDownregulated = point.logFC <= -foldChangeThreshold && point.pValue < pValueThreshold;
    
    let color = '#888888';  // Gray for not significant
    if (isUpregulated) color = '#dc2626';  // Red for upregulated
    if (isDownregulated) color = '#2563eb';  // Blue for downregulated
    
    return {
      x: point.logFC,
      y: negLogP,
      label: point.label,
      color,
      size: (isUpregulated || isDownregulated) ? 5 : 3,
    };
  });

  // Custom click handler
  const handlePointClick = (scatterPoint: ScatterPoint, index: number) => {
    if (onPointClick) {
      onPointClick(data[index]);
    }
  };

  return (
    <div className="relative">
      {/* Plot */}
      <ScatterPlot
        data={scatterData}
        width={width}
        height={height}
        xLabel="Log2 Fold Change"
        yLabel="-Log10 P-value"
        title={title}
        showGrid={true}
        onPointClick={handlePointClick}
      />
      
      {/* Legend */}
      <div className="mt-4 flex justify-center gap-6 text-sm">
        <div className="flex items-center gap-2">
          <div className="w-3 h-3 rounded-full bg-red-600"></div>
          <span>Upregulated ({scatterData.filter(p => p.color === '#dc2626').length})</span>
        </div>
        <div className="flex items-center gap-2">
          <div className="w-3 h-3 rounded-full bg-blue-600"></div>
          <span>Downregulated ({scatterData.filter(p => p.color === '#2563eb').length})</span>
        </div>
        <div className="flex items-center gap-2">
          <div className="w-3 h-3 rounded-full bg-gray-500"></div>
          <span>Not significant ({scatterData.filter(p => p.color === '#888888').length})</span>
        </div>
      </div>
      
      {/* Thresholds info */}
      <div className="mt-2 text-center text-xs text-gray-600">
        Thresholds: |Log2FC| ≥ {foldChangeThreshold}, P-value &lt; {pValueThreshold}
      </div>
      
      {/* Statistics */}
      <div className="mt-4 p-3 bg-gray-50 rounded border border-gray-200 text-sm">
        <div className="font-semibold mb-2">Statistics</div>
        <div className="grid grid-cols-3 gap-4">
          <div>
            <div className="text-gray-600">Total features</div>
            <div className="text-lg font-bold">{data.length}</div>
          </div>
          <div>
            <div className="text-gray-600">Significant</div>
            <div className="text-lg font-bold text-green-600">
              {scatterData.filter(p => p.color !== '#888888').length}
            </div>
          </div>
          <div>
            <div className="text-gray-600">% Significant</div>
            <div className="text-lg font-bold text-blue-600">
              {((scatterData.filter(p => p.color !== '#888888').length / data.length) * 100).toFixed(1)}%
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
