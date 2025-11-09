import React, { useRef, useEffect } from 'react';

export interface HeatmapData {
  rows: string[];  // Row labels (e.g., genes)
  columns: string[];  // Column labels (e.g., samples)
  values: number[][];  // Matrix values
  min?: number;  // Manual min value for color scale
  max?: number;  // Manual max value for color scale
}

export interface HeatmapProps {
  data: HeatmapData;
  width?: number;
  height?: number;
  cellSize?: number;
  colorScheme?: 'reds' | 'blues' | 'greens' | 'diverging';
  showLabels?: boolean;
  onCellClick?: (row: number, col: number, value: number) => void;
}

/**
 * Heatmap visualization component for displaying matrix data.
 * 
 * Common use cases:
 * - Gene expression matrices
 * - Protein abundance across samples
 * - Correlation matrices
 * - Distance matrices
 */
export function Heatmap({
  data,
  width = 800,
  height = 600,
  cellSize = 20,
  colorScheme = 'diverging',
  showLabels = true,
  onCellClick,
}: HeatmapProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [hoveredCell, setHoveredCell] = React.useState<{
    row: number;
    col: number;
    value: number;
  } | null>(null);

  // Calculate color scale
  const getColor = (value: number, min: number, max: number): string => {
    // Normalize value to 0-1
    const normalized = (value - min) / (max - min);
    
    switch (colorScheme) {
      case 'reds':
        return `rgb(${Math.floor(255 * normalized)}, 0, 0)`;
      
      case 'blues':
        return `rgb(0, 0, ${Math.floor(255 * normalized)})`;
      
      case 'greens':
        return `rgb(0, ${Math.floor(255 * normalized)}, 0)`;
      
      case 'diverging':
      default:
        // Blue-White-Red for diverging data (centered at 0)
        if (value < 0) {
          const intensity = Math.abs(value) / Math.abs(min);
          return `rgb(${Math.floor(255 * (1 - intensity))}, ${Math.floor(255 * (1 - intensity))}, 255)`;
        } else {
          const intensity = value / max;
          return `rgb(255, ${Math.floor(255 * (1 - intensity))}, ${Math.floor(255 * (1 - intensity))})`;
        }
    }
  };

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !data.values.length) return;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    // Calculate dimensions
    const numRows = data.rows.length;
    const numCols = data.columns.length;
    const labelWidth = showLabels ? 100 : 0;
    const labelHeight = showLabels ? 50 : 0;
    
    const plotWidth = Math.min(numCols * cellSize, width - labelWidth);
    const plotHeight = Math.min(numRows * cellSize, height - labelHeight);
    const actualCellWidth = plotWidth / numCols;
    const actualCellHeight = plotHeight / numRows;

    // Clear canvas
    ctx.clearRect(0, 0, width, height);

    // Find min/max values
    let minVal = data.min !== undefined ? data.min : Infinity;
    let maxVal = data.max !== undefined ? data.max : -Infinity;
    
    if (data.min === undefined || data.max === undefined) {
      data.values.forEach(row => {
        row.forEach(val => {
          if (val < minVal) minVal = val;
          if (val > maxVal) maxVal = val;
        });
      });
    }

    // Draw heatmap cells
    data.values.forEach((row, rowIdx) => {
      row.forEach((value, colIdx) => {
        const x = labelWidth + colIdx * actualCellWidth;
        const y = labelHeight + rowIdx * actualCellHeight;
        
        ctx.fillStyle = getColor(value, minVal, maxVal);
        ctx.fillRect(x, y, actualCellWidth, actualCellHeight);
        
        // Draw cell border
        ctx.strokeStyle = '#ddd';
        ctx.lineWidth = 0.5;
        ctx.strokeRect(x, y, actualCellWidth, actualCellHeight);
      });
    });

    // Draw row labels
    if (showLabels) {
      ctx.fillStyle = '#333';
      ctx.font = '12px sans-serif';
      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      
      data.rows.forEach((label, idx) => {
        const y = labelHeight + idx * actualCellHeight + actualCellHeight / 2;
        ctx.fillText(
          label.length > 15 ? label.substring(0, 12) + '...' : label,
          labelWidth - 5,
          y
        );
      });
    }

    // Draw column labels
    if (showLabels) {
      ctx.save();
      ctx.fillStyle = '#333';
      ctx.font = '12px sans-serif';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';
      
      data.columns.forEach((label, idx) => {
        const x = labelWidth + idx * actualCellWidth + actualCellWidth / 2;
        const y = labelHeight - 5;
        
        ctx.translate(x, y);
        ctx.rotate(-Math.PI / 4);
        ctx.fillText(
          label.length > 15 ? label.substring(0, 12) + '...' : label,
          0,
          0
        );
        ctx.setTransform(1, 0, 0, 1, 0, 0);
      });
      
      ctx.restore();
    }

    // Draw color scale legend
    const legendWidth = 200;
    const legendHeight = 20;
    const legendX = width - legendWidth - 20;
    const legendY = height - legendHeight - 30;
    
    // Draw gradient
    for (let i = 0; i < legendWidth; i++) {
      const value = minVal + (i / legendWidth) * (maxVal - minVal);
      ctx.fillStyle = getColor(value, minVal, maxVal);
      ctx.fillRect(legendX + i, legendY, 1, legendHeight);
    }
    
    // Draw legend labels
    ctx.fillStyle = '#333';
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(minVal.toFixed(2), legendX, legendY + legendHeight + 15);
    ctx.fillText(maxVal.toFixed(2), legendX + legendWidth, legendY + legendHeight + 15);
    ctx.fillText(((minVal + maxVal) / 2).toFixed(2), legendX + legendWidth / 2, legendY + legendHeight + 15);

  }, [data, width, height, cellSize, colorScheme, showLabels]);

  // Handle mouse move for hover tooltip
  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;

    const labelWidth = showLabels ? 100 : 0;
    const labelHeight = showLabels ? 50 : 0;
    const numCols = data.columns.length;
    const numRows = data.rows.length;
    const plotWidth = Math.min(numCols * cellSize, width - labelWidth);
    const plotHeight = Math.min(numRows * cellSize, height - labelHeight);
    const actualCellWidth = plotWidth / numCols;
    const actualCellHeight = plotHeight / numRows;

    // Check if mouse is over a cell
    if (x >= labelWidth && x < labelWidth + plotWidth &&
        y >= labelHeight && y < labelHeight + plotHeight) {
      const col = Math.floor((x - labelWidth) / actualCellWidth);
      const row = Math.floor((y - labelHeight) / actualCellHeight);
      
      if (row >= 0 && row < numRows && col >= 0 && col < numCols) {
        setHoveredCell({
          row,
          col,
          value: data.values[row][col],
        });
        return;
      }
    }

    setHoveredCell(null);
  };

  const handleMouseLeave = () => {
    setHoveredCell(null);
  };

  const handleClick = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (!hoveredCell || !onCellClick) return;
    onCellClick(hoveredCell.row, hoveredCell.col, hoveredCell.value);
  };

  return (
    <div className="relative inline-block">
      <canvas
        ref={canvasRef}
        width={width}
        height={height}
        onMouseMove={handleMouseMove}
        onMouseLeave={handleMouseLeave}
        onClick={handleClick}
        className={onCellClick ? 'cursor-pointer' : ''}
        style={{ border: '1px solid #ddd' }}
      />
      
      {/* Hover tooltip */}
      {hoveredCell && (
        <div className="absolute bg-gray-900 text-white text-sm px-3 py-2 rounded shadow-lg pointer-events-none"
             style={{
               left: `${(hoveredCell.col + 1) * cellSize + (showLabels ? 100 : 0)}px`,
               top: `${(hoveredCell.row + 1) * cellSize + (showLabels ? 50 : 0)}px`,
             }}>
          <div><strong>{data.rows[hoveredCell.row]}</strong></div>
          <div>{data.columns[hoveredCell.col]}</div>
          <div>Value: {hoveredCell.value.toFixed(4)}</div>
        </div>
      )}
    </div>
  );
}
