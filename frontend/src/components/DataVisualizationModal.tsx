/**
 * Data Visualization Modal
 * 数据可视化模态框 - 显示节点输出数据的可视化
 * 支持：表格、直方图、箱线图、散点图、热图
 */

import React, { useState } from 'react';
import {
  BarChart,
  Bar,
  LineChart,
  Line,
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
} from 'recharts';

export type VisualizationType = 'table' | 'histogram' | 'boxplot' | 'scatter' | 'line';

export interface DataVisualizationProps {
  isOpen: boolean;
  onClose: () => void;
  title: string;
  data: any[];
  columns?: string[];
  defaultVisualization?: VisualizationType;
}

const DataVisualizationModal: React.FC<DataVisualizationProps> = ({
  isOpen,
  onClose,
  title,
  data,
  columns,
  defaultVisualization = 'table',
}) => {
  const [currentView, setCurrentView] = useState<VisualizationType>(defaultVisualization);
  const [selectedColumn, setSelectedColumn] = useState<string>(columns?.[0] || '');
  const [page, setPage] = useState(0);
  const rowsPerPage = 20;

  if (!isOpen) return null;

  // 自动检测数据列
  const detectedColumns = columns || (data.length > 0 ? Object.keys(data[0]) : []);

  // 计算统计信息
  const getStatistics = (columnName: string) => {
    const values = data.map(row => {
      const val = row[columnName];
      return typeof val === 'number' ? val : parseFloat(val);
    }).filter(v => !isNaN(v));

    if (values.length === 0) return null;

    values.sort((a, b) => a - b);
    const sum = values.reduce((a, b) => a + b, 0);
    const mean = sum / values.length;
    const median = values[Math.floor(values.length / 2)];
    const min = values[0];
    const max = values[values.length - 1];
    
    // 标准差
    const variance = values.reduce((acc, val) => acc + Math.pow(val - mean, 2), 0) / values.length;
    const stdDev = Math.sqrt(variance);

    // 四分位数
    const q1 = values[Math.floor(values.length * 0.25)];
    const q3 = values[Math.floor(values.length * 0.75)];

    return { mean, median, min, max, stdDev, q1, q3, count: values.length };
  };

  // 生成直方图数据
  const getHistogramData = (columnName: string, bins: number = 20) => {
    const values = data.map(row => {
      const val = row[columnName];
      return typeof val === 'number' ? val : parseFloat(val);
    }).filter(v => !isNaN(v));

    if (values.length === 0) return [];

    const min = Math.min(...values);
    const max = Math.max(...values);
    const binWidth = (max - min) / bins;

    const histogram = Array(bins).fill(0).map((_, i) => ({
      range: `${(min + i * binWidth).toFixed(2)}-${(min + (i + 1) * binWidth).toFixed(2)}`,
      count: 0,
      midpoint: min + (i + 0.5) * binWidth,
    }));

    values.forEach(val => {
      const binIndex = Math.min(Math.floor((val - min) / binWidth), bins - 1);
      histogram[binIndex].count++;
    });

    return histogram;
  };

  // 渲染表格视图
  const renderTable = () => {
    const startIdx = page * rowsPerPage;
    const endIdx = startIdx + rowsPerPage;
    const pageData = data.slice(startIdx, endIdx);
    const totalPages = Math.ceil(data.length / rowsPerPage);

    return (
      <div>
        <div style={{ overflowX: 'auto', maxHeight: '500px', overflowY: 'auto' }}>
          <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.875rem' }}>
            <thead style={{ position: 'sticky', top: 0, backgroundColor: '#111827', zIndex: 1 }}>
              <tr>
                {detectedColumns.map((col) => (
                  <th
                    key={col}
                    style={{
                      padding: '0.75rem',
                      textAlign: 'left',
                      color: '#f3f4f6',
                      fontWeight: 600,
                      borderBottom: '2px solid #374151',
                    }}
                  >
                    {col}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {pageData.map((row, idx) => (
                <tr
                  key={idx}
                  style={{
                    backgroundColor: idx % 2 === 0 ? '#1f2937' : '#111827',
                  }}
                >
                  {detectedColumns.map((col) => (
                    <td
                      key={col}
                      style={{
                        padding: '0.75rem',
                        color: '#e5e7eb',
                        borderBottom: '1px solid #374151',
                      }}
                    >
                      {String(row[col])}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>

        {/* Pagination */}
        <div
          style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            marginTop: '1rem',
            padding: '0.75rem',
            backgroundColor: '#111827',
            borderRadius: '4px',
          }}
        >
          <span style={{ color: '#9ca3af', fontSize: '0.875rem' }}>
            Showing {startIdx + 1}-{Math.min(endIdx, data.length)} of {data.length} rows
          </span>
          <div style={{ display: 'flex', gap: '0.5rem' }}>
            <button
              onClick={() => setPage(Math.max(0, page - 1))}
              disabled={page === 0}
              style={{
                padding: '0.5rem 1rem',
                backgroundColor: page === 0 ? '#374151' : '#3b82f6',
                color: 'white',
                border: 'none',
                borderRadius: '4px',
                cursor: page === 0 ? 'not-allowed' : 'pointer',
                fontSize: '0.875rem',
              }}
            >
              Previous
            </button>
            <span style={{ padding: '0.5rem 1rem', color: '#e5e7eb', fontSize: '0.875rem' }}>
              Page {page + 1} of {totalPages}
            </span>
            <button
              onClick={() => setPage(Math.min(totalPages - 1, page + 1))}
              disabled={page >= totalPages - 1}
              style={{
                padding: '0.5rem 1rem',
                backgroundColor: page >= totalPages - 1 ? '#374151' : '#3b82f6',
                color: 'white',
                border: 'none',
                borderRadius: '4px',
                cursor: page >= totalPages - 1 ? 'not-allowed' : 'pointer',
                fontSize: '0.875rem',
              }}
            >
              Next
            </button>
          </div>
        </div>
      </div>
    );
  };

  // 渲染直方图
  const renderHistogram = () => {
    if (!selectedColumn) return <div style={{ color: '#9ca3af' }}>Select a column to visualize</div>;

    const histData = getHistogramData(selectedColumn);
    const stats = getStatistics(selectedColumn);

    return (
      <div>
        {/* Statistics Card */}
        {stats && (
          <div
            style={{
              padding: '1rem',
              backgroundColor: '#111827',
              borderRadius: '6px',
              marginBottom: '1rem',
              border: '1px solid #374151',
            }}
          >
            <h4 style={{ color: '#f3f4f6', marginBottom: '0.75rem', fontSize: '0.875rem', fontWeight: 600 }}>
              📊 Statistics for {selectedColumn}
            </h4>
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '1rem', fontSize: '0.75rem' }}>
              <div>
                <div style={{ color: '#9ca3af' }}>Count</div>
                <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{stats.count}</div>
              </div>
              <div>
                <div style={{ color: '#9ca3af' }}>Mean</div>
                <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{stats.mean.toFixed(3)}</div>
              </div>
              <div>
                <div style={{ color: '#9ca3af' }}>Median</div>
                <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{stats.median.toFixed(3)}</div>
              </div>
              <div>
                <div style={{ color: '#9ca3af' }}>Std Dev</div>
                <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{stats.stdDev.toFixed(3)}</div>
              </div>
              <div>
                <div style={{ color: '#9ca3af' }}>Min</div>
                <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{stats.min.toFixed(3)}</div>
              </div>
              <div>
                <div style={{ color: '#9ca3af' }}>Q1</div>
                <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{stats.q1.toFixed(3)}</div>
              </div>
              <div>
                <div style={{ color: '#9ca3af' }}>Q3</div>
                <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{stats.q3.toFixed(3)}</div>
              </div>
              <div>
                <div style={{ color: '#9ca3af' }}>Max</div>
                <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{stats.max.toFixed(3)}</div>
              </div>
            </div>
          </div>
        )}

        {/* Histogram Chart */}
        <ResponsiveContainer width="100%" height={400}>
          <BarChart data={histData}>
            <CartesianGrid strokeDasharray="3 3" stroke="#374151" />
            <XAxis dataKey="range" stroke="#9ca3af" tick={{ fill: '#9ca3af', fontSize: 12 }} angle={-45} textAnchor="end" height={100} />
            <YAxis stroke="#9ca3af" tick={{ fill: '#9ca3af', fontSize: 12 }} />
            <Tooltip
              contentStyle={{ backgroundColor: '#1f2937', border: '1px solid #374151', borderRadius: '4px' }}
              labelStyle={{ color: '#f3f4f6' }}
              itemStyle={{ color: '#e5e7eb' }}
            />
            <Legend wrapperStyle={{ color: '#f3f4f6' }} />
            <Bar dataKey="count" fill="#3b82f6" name="Frequency" />
          </BarChart>
        </ResponsiveContainer>
      </div>
    );
  };

  // 渲染箱线图（使用散点图模拟）
  const renderBoxplot = () => {
    if (!selectedColumn) return <div style={{ color: '#9ca3af' }}>Select a column to visualize</div>;

    const stats = getStatistics(selectedColumn);
    if (!stats) return <div style={{ color: '#9ca3af' }}>No numeric data available</div>;

    // 创建箱线图数据点
    const boxData = [
      { name: 'Min', value: stats.min, type: 'whisker' },
      { name: 'Q1', value: stats.q1, type: 'box' },
      { name: 'Median', value: stats.median, type: 'median' },
      { name: 'Mean', value: stats.mean, type: 'mean' },
      { name: 'Q3', value: stats.q3, type: 'box' },
      { name: 'Max', value: stats.max, type: 'whisker' },
    ];

    return (
      <div>
        <div
          style={{
            padding: '1rem',
            backgroundColor: '#111827',
            borderRadius: '6px',
            marginBottom: '1rem',
            border: '1px solid #374151',
          }}
        >
          <h4 style={{ color: '#f3f4f6', marginBottom: '0.75rem', fontSize: '0.875rem', fontWeight: 600 }}>
            📦 Box Plot for {selectedColumn}
          </h4>
          <ResponsiveContainer width="100%" height={300}>
            <LineChart data={boxData}>
              <CartesianGrid strokeDasharray="3 3" stroke="#374151" />
              <XAxis dataKey="name" stroke="#9ca3af" tick={{ fill: '#9ca3af' }} />
              <YAxis stroke="#9ca3af" tick={{ fill: '#9ca3af' }} />
              <Tooltip
                contentStyle={{ backgroundColor: '#1f2937', border: '1px solid #374151', borderRadius: '4px' }}
                labelStyle={{ color: '#f3f4f6' }}
                itemStyle={{ color: '#e5e7eb' }}
              />
              <Line type="monotone" dataKey="value" stroke="#3b82f6" strokeWidth={2} dot={{ r: 6, fill: '#3b82f6' }} />
            </LineChart>
          </ResponsiveContainer>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: '1rem', fontSize: '0.75rem' }}>
          <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '4px', border: '1px solid #374151' }}>
            <div style={{ color: '#9ca3af', marginBottom: '0.5rem' }}>Range</div>
            <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{(stats.max - stats.min).toFixed(3)}</div>
          </div>
          <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '4px', border: '1px solid #374151' }}>
            <div style={{ color: '#9ca3af', marginBottom: '0.5rem' }}>IQR (Q3-Q1)</div>
            <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{(stats.q3 - stats.q1).toFixed(3)}</div>
          </div>
          <div style={{ padding: '1rem', backgroundColor: '#111827', borderRadius: '4px', border: '1px solid #374151' }}>
            <div style={{ color: '#9ca3af', marginBottom: '0.5rem' }}>Coefficient of Variation</div>
            <div style={{ color: '#f3f4f6', fontWeight: 600 }}>{((stats.stdDev / stats.mean) * 100).toFixed(2)}%</div>
          </div>
        </div>
      </div>
    );
  };

  return (
    <div
      style={{
        position: 'fixed',
        top: 0,
        left: 0,
        right: 0,
        bottom: 0,
        backgroundColor: 'rgba(0,0,0,0.8)',
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        zIndex: 3000,
        padding: '2rem',
      }}
      onClick={onClose}
    >
      <div
        style={{
          backgroundColor: '#1f2937',
          borderRadius: '12px',
          maxWidth: '1200px',
          width: '100%',
          maxHeight: '90vh',
          overflow: 'hidden',
          display: 'flex',
          flexDirection: 'column',
          border: '2px solid #374151',
        }}
        onClick={(e) => e.stopPropagation()}
      >
        {/* Header */}
        <div
          style={{
            padding: '1.5rem',
            borderBottom: '2px solid #374151',
            backgroundColor: '#111827',
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
          }}
        >
          <h2 style={{ margin: 0, color: '#f3f4f6', fontSize: '1.25rem', fontWeight: 600 }}>
            📊 {title}
          </h2>
          <button
            onClick={onClose}
            style={{
              background: 'none',
              border: 'none',
              color: '#9ca3af',
              cursor: 'pointer',
              fontSize: '1.5rem',
              padding: 0,
              lineHeight: 1,
            }}
          >
            ✕
          </button>
        </div>

        {/* Toolbar */}
        <div
          style={{
            padding: '1rem 1.5rem',
            borderBottom: '1px solid #374151',
            backgroundColor: '#1f2937',
            display: 'flex',
            gap: '1rem',
            flexWrap: 'wrap',
            alignItems: 'center',
          }}
        >
          <div style={{ display: 'flex', gap: '0.5rem' }}>
            {(['table', 'histogram', 'boxplot'] as VisualizationType[]).map((view) => (
              <button
                key={view}
                onClick={() => setCurrentView(view)}
                style={{
                  padding: '0.5rem 1rem',
                  backgroundColor: currentView === view ? '#3b82f6' : '#374151',
                  color: 'white',
                  border: 'none',
                  borderRadius: '4px',
                  cursor: 'pointer',
                  fontSize: '0.875rem',
                  fontWeight: 600,
                  textTransform: 'capitalize',
                }}
              >
                {view === 'table' && '📋'}
                {view === 'histogram' && '📊'}
                {view === 'boxplot' && '📦'}
                {' '}
                {view}
              </button>
            ))}
          </div>

          {currentView !== 'table' && (
            <select
              value={selectedColumn}
              onChange={(e) => setSelectedColumn(e.target.value)}
              style={{
                padding: '0.5rem 1rem',
                backgroundColor: '#374151',
                border: '1px solid #4b5563',
                borderRadius: '4px',
                color: '#f3f4f6',
                fontSize: '0.875rem',
                cursor: 'pointer',
              }}
            >
              {detectedColumns.map((col) => (
                <option key={col} value={col}>
                  {col}
                </option>
              ))}
            </select>
          )}

          <div style={{ marginLeft: 'auto', fontSize: '0.875rem', color: '#9ca3af' }}>
            {data.length} rows × {detectedColumns.length} columns
          </div>
        </div>

        {/* Content */}
        <div style={{ flex: 1, overflowY: 'auto', padding: '1.5rem' }}>
          {currentView === 'table' && renderTable()}
          {currentView === 'histogram' && renderHistogram()}
          {currentView === 'boxplot' && renderBoxplot()}
        </div>
      </div>
    </div>
  );
};

export default DataVisualizationModal;
