import React from 'react';
import { ScatterPlot, ScatterPoint } from './ScatterPlot';

export interface MAPoint {
  label: string;  // Gene/Protein ID
  baseMean: number;  // Average expression
  logFC: number;  // Log2 fold change
  pValue: number;  // P-value
}

export interface MAPlotProps {
  data: MAPoint[];
  width?: number;
  height?: number;
  foldChangeThreshold?: number;  // Fold change threshold (default: 1)
  pValueThreshold?: number;  // P-value threshold (default: 0.05)
  title?: string;
  onPointClick?: (point: MAPoint) => void;
}

/**
 * MA Plot for differential expression analysis.
 * 
 * X-axis: Log2 mean expression (A = average)
 * Y-axis: Log2 fold change (M = minus)
 * 
 * Points are colored by significance and direction:
 * - Red: significantly upregulated
 * - Blue: significantly downregulated
 * - Gray: not significant
 */
export function MAPlot({
  data,
  width = 700,
  height = 600,
  foldChangeThreshold = 1,
  pValueThreshold = 0.05,
  title = 'MA Plot',
  onPointClick,
}: MAPlotProps) {
  // Convert to scatter plot format
  const scatterData: ScatterPoint[] = data.map(point => {
    const logBaseMean = Math.log2(point.baseMean + 1);  // Add 1 to avoid log(0)
    const isUpregulated = point.logFC >= foldChangeThreshold && point.pValue < pValueThreshold;
    const isDownregulated = point.logFC <= -foldChangeThreshold && point.pValue < pValueThreshold;
    
    let color = '#888888';  // Gray for not significant
    if (isUpregulated) color = '#dc2626';  // Red for upregulated
    if (isDownregulated) color = '#2563eb';  // Blue for downregulated
    
    return {
      x: logBaseMean,
      y: point.logFC,
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
        xLabel="Log2 Mean Expression"
        yLabel="Log2 Fold Change"
        title={title}
        showGrid={true}
        onPointClick={handlePointClick}
      />
      
      {/* Threshold lines indicator */}
      <div className="mt-2 text-center text-xs text-gray-600">
        Horizontal lines indicate fold change thresholds: ±{foldChangeThreshold}
      </div>
      
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
      
      {/* Statistics */}
      <div className="mt-4 p-3 bg-gray-50 rounded border border-gray-200 text-sm">
        <div className="font-semibold mb-2">Expression Statistics</div>
        <div className="grid grid-cols-2 gap-4">
          <div>
            <div className="text-gray-600">Mean expression range</div>
            <div className="font-mono text-sm">
              {Math.min(...data.map(p => p.baseMean)).toFixed(2)} - {Math.max(...data.map(p => p.baseMean)).toFixed(2)}
            </div>
          </div>
          <div>
            <div className="text-gray-600">Fold change range</div>
            <div className="font-mono text-sm">
              {Math.min(...data.map(p => p.logFC)).toFixed(2)} - {Math.max(...data.map(p => p.logFC)).toFixed(2)}
            </div>
          </div>
          <div>
            <div className="text-gray-600">Significant features</div>
            <div className="text-lg font-bold text-green-600">
              {scatterData.filter(p => p.color !== '#888888').length} / {data.length}
            </div>
          </div>
          <div>
            <div className="text-gray-600">Up/Down ratio</div>
            <div className="text-lg font-bold text-purple-600">
              {(scatterData.filter(p => p.color === '#dc2626').length / 
                Math.max(scatterData.filter(p => p.color === '#2563eb').length, 1)).toFixed(2)}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
