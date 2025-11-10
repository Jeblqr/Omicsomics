import React, { useRef, useEffect } from 'react';

export interface ScatterPoint {
  x: number;
  y: number;
  label?: string;
  color?: string;
  size?: number;
}

export interface ScatterPlotProps {
  data: ScatterPoint[];
  width?: number;
  height?: number;
  xLabel?: string;
  yLabel?: string;
  title?: string;
  showGrid?: boolean;
  logScaleX?: boolean;
  logScaleY?: boolean;
  onPointClick?: (point: ScatterPoint, index: number) => void;
}

/**
 * Scatter plot component for displaying 2D point data.
 * 
 * Common use cases:
 * - Principal Component Analysis (PCA)
 * - t-SNE visualization
 * - Quality control plots
 * - Correlation plots
 */
export function ScatterPlot({
  data,
  width = 600,
  height = 500,
  xLabel = 'X',
  yLabel = 'Y',
  title,
  showGrid = true,
  logScaleX = false,
  logScaleY = false,
  onPointClick,
}: ScatterPlotProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [hoveredPoint, setHoveredPoint] = React.useState<{
    point: ScatterPoint;
    index: number;
  } | null>(null);

  const padding = {
    top: title ? 50 : 30,
    right: 30,
    bottom: 60,
    left: 70,
  };

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !data.length) return;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    // Clear canvas
    ctx.clearRect(0, 0, width, height);

    // Calculate data ranges
    let xMin = Infinity, xMax = -Infinity;
    let yMin = Infinity, yMax = -Infinity;

    data.forEach(point => {
      const x = logScaleX ? Math.log10(point.x) : point.x;
      const y = logScaleY ? Math.log10(point.y) : point.y;
      
      if (x < xMin) xMin = x;
      if (x > xMax) xMax = x;
      if (y < yMin) yMin = y;
      if (y > yMax) yMax = y;
    });

    // Add padding to ranges
    const xRange = xMax - xMin;
    const yRange = yMax - yMin;
    xMin -= xRange * 0.05;
    xMax += xRange * 0.05;
    yMin -= yRange * 0.05;
    yMax += yRange * 0.05;

    // Calculate plot area
    const plotWidth = width - padding.left - padding.right;
    const plotHeight = height - padding.top - padding.bottom;

    // Scale functions
    const scaleX = (x: number) => {
      const value = logScaleX ? Math.log10(x) : x;
      return padding.left + ((value - xMin) / (xMax - xMin)) * plotWidth;
    };

    const scaleY = (y: number) => {
      const value = logScaleY ? Math.log10(y) : y;
      return height - padding.bottom - ((value - yMin) / (yMax - yMin)) * plotHeight;
    };

    // Draw grid
    if (showGrid) {
      ctx.strokeStyle = '#e0e0e0';
      ctx.lineWidth = 0.5;

      // Vertical grid lines
      for (let i = 0; i <= 10; i++) {
        const x = padding.left + (i / 10) * plotWidth;
        ctx.beginPath();
        ctx.moveTo(x, padding.top);
        ctx.lineTo(x, height - padding.bottom);
        ctx.stroke();
      }

      // Horizontal grid lines
      for (let i = 0; i <= 10; i++) {
        const y = padding.top + (i / 10) * plotHeight;
        ctx.beginPath();
        ctx.moveTo(padding.left, y);
        ctx.lineTo(width - padding.right, y);
        ctx.stroke();
      }
    }

    // Draw axes
    ctx.strokeStyle = '#333';
    ctx.lineWidth = 2;
    
    // X-axis
    ctx.beginPath();
    ctx.moveTo(padding.left, height - padding.bottom);
    ctx.lineTo(width - padding.right, height - padding.bottom);
    ctx.stroke();
    
    // Y-axis
    ctx.beginPath();
    ctx.moveTo(padding.left, padding.top);
    ctx.lineTo(padding.left, height - padding.bottom);
    ctx.stroke();

    // Draw axis labels
    ctx.fillStyle = '#333';
    ctx.font = '14px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    
    // X-axis label
    ctx.fillText(xLabel, width / 2, height - 20);
    
    // Y-axis label (rotated)
    ctx.save();
    ctx.translate(15, height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(yLabel, 0, 0);
    ctx.restore();

    // Draw title
    if (title) {
      ctx.font = 'bold 16px sans-serif';
      ctx.fillText(title, width / 2, 25);
    }

    // Draw tick marks and labels
    ctx.font = '11px sans-serif';
    ctx.fillStyle = '#374151';
    
    // X-axis ticks
    for (let i = 0; i <= 5; i++) {
      const value = xMin + (i / 5) * (xMax - xMin);
      const x = padding.left + (i / 5) * plotWidth;
      
      ctx.beginPath();
      ctx.moveTo(x, height - padding.bottom);
      ctx.lineTo(x, height - padding.bottom + 5);
      ctx.stroke();
      
      const label = logScaleX ? Math.pow(10, value).toExponential(1) : value.toFixed(1);
      ctx.fillText(label, x, height - padding.bottom + 15);
    }

    // Y-axis ticks
    for (let i = 0; i <= 5; i++) {
      const value = yMin + (i / 5) * (yRange);
      const y = height - padding.bottom - (i / 5) * plotHeight;
      
      ctx.beginPath();
      ctx.moveTo(padding.left - 5, y);
      ctx.lineTo(padding.left, y);
      ctx.stroke();
      
      const label = logScaleY ? Math.pow(10, value).toExponential(1) : value.toFixed(1);
      ctx.textAlign = 'right';
      ctx.fillText(label, padding.left - 10, y);
    }

    // Draw points
    data.forEach((point, idx) => {
      const x = scaleX(point.x);
      const y = scaleY(point.y);
      const size = point.size || 4;
      const color = point.color || '#4299e1';

      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(x, y, size, 0, 2 * Math.PI);
      ctx.fill();

      // Highlight hovered point
      if (hoveredPoint && hoveredPoint.index === idx) {
        ctx.strokeStyle = '#000';
        ctx.lineWidth = 2;
        ctx.stroke();
      }
    });

  }, [data, width, height, xLabel, yLabel, title, showGrid, logScaleX, logScaleY, hoveredPoint]);

  // Handle mouse interactions
  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const rect = canvas.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    // Calculate plot area
    const plotWidth = width - padding.left - padding.right;
    const plotHeight = height - padding.top - padding.bottom;

    // Find data ranges
    let xMin = Infinity, xMax = -Infinity;
    let yMin = Infinity, yMax = -Infinity;

    data.forEach(point => {
      const x = logScaleX ? Math.log10(point.x) : point.x;
      const y = logScaleY ? Math.log10(point.y) : point.y;
      
      if (x < xMin) xMin = x;
      if (x > xMax) xMax = x;
      if (y < yMin) yMin = y;
      if (y > yMax) yMax = y;
    });

    const xRange = xMax - xMin;
    const yRange = yMax - yMin;
    xMin -= xRange * 0.05;
    xMax += xRange * 0.05;
    yMin -= yRange * 0.05;
    yMax += yRange * 0.05;

    // Check if mouse is near any point
    let closestPoint: { point: ScatterPoint; index: number; distance: number } | undefined;

    data.forEach((point, idx) => {
      const xValue = logScaleX ? Math.log10(point.x) : point.x;
      const yValue = logScaleY ? Math.log10(point.y) : point.y;
      
      const x = padding.left + ((xValue - xMin) / (xMax - xMin)) * plotWidth;
      const y = height - padding.bottom - ((yValue - yMin) / (yMax - yMin)) * plotHeight;
      
      const distance = Math.sqrt(Math.pow(mouseX - x, 2) + Math.pow(mouseY - y, 2));
      
      if (distance < 10 && (!closestPoint || distance < closestPoint.distance)) {
        closestPoint = { point, index: idx, distance };
      }
    });

    if (closestPoint) {
      setHoveredPoint({ point: closestPoint.point, index: closestPoint.index });
    } else {
      setHoveredPoint(null);
    }
  };

  const handleMouseLeave = () => {
    setHoveredPoint(null);
  };

  const handleClick = () => {
    if (hoveredPoint && onPointClick) {
      onPointClick(hoveredPoint.point, hoveredPoint.index);
    }
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
        className={onPointClick ? 'cursor-pointer' : ''}
        style={{ border: '1px solid #ddd' }}
      />
      
      {/* Hover tooltip */}
      {hoveredPoint && (
        <div className="absolute bg-gray-900 text-white text-xs px-2 py-1 rounded shadow-lg pointer-events-none"
             style={{
               left: `${(hoveredPoint.point.x / data[0].x) * width}px`,
               top: `${height - (hoveredPoint.point.y / data[0].y) * height}px`,
               transform: 'translate(-50%, -120%)',
             }}>
          {hoveredPoint.point.label && <div className="font-bold">{hoveredPoint.point.label}</div>}
          <div>X: {hoveredPoint.point.x.toFixed(3)}</div>
          <div>Y: {hoveredPoint.point.y.toFixed(3)}</div>
        </div>
      )}
    </div>
  );
}
