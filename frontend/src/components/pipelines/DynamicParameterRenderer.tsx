/**
 * Dynamic Parameter Renderer
 * 动态参数渲染器 - 根据Tool JSON Schema自动生成参数配置UI
 * 
 * 支持所有参数类型的自动渲染：
 * - string, number, integer, boolean
 * - select, multiselect
 * - file, directory
 * - color, date, range
 * - textarea, json, array
 */

import React, { useState, useEffect } from 'react';
import { ToolParameterSchema } from '../../schemas/ToolSchema';

interface DynamicParameterRendererProps {
  parameter: ToolParameterSchema;
  value: any;
  onChange: (value: any) => void;
  disabled?: boolean;
}

const DynamicParameterRenderer: React.FC<DynamicParameterRendererProps> = ({
  parameter,
  value,
  onChange,
  disabled = false,
}) => {
  const [localValue, setLocalValue] = useState(value ?? parameter.default);

  useEffect(() => {
    setLocalValue(value ?? parameter.default);
  }, [value, parameter.default]);

  const handleChange = (newValue: any) => {
    setLocalValue(newValue);
    onChange(newValue);
  };

  // ===== 渲染标签和描述 =====
  const renderLabel = () => (
    <div style={{ marginBottom: '0.5rem' }}>
      <label
        style={{
          fontSize: '0.875rem',
          fontWeight: 600,
          color: '#e5e7eb',
          display: 'flex',
          alignItems: 'center',
          gap: '0.5rem',
        }}
      >
        {parameter.label}
        {parameter.required && <span style={{ color: '#ef4444' }}>*</span>}
        {parameter.unit && (
          <span style={{ fontSize: '0.75rem', color: '#9ca3af', fontWeight: 400 }}>
            ({parameter.unit})
          </span>
        )}
      </label>
      {parameter.description && (
        <div style={{ fontSize: '0.75rem', color: '#9ca3af', marginTop: '0.25rem' }}>
          {parameter.description}
        </div>
      )}
    </div>
  );

  // ===== 渲染帮助文本 =====
  const renderHelp = () => {
    if (!parameter.helpText) return null;
    return (
      <div
        style={{
          fontSize: '0.75rem',
          color: '#6b7280',
          marginTop: '0.25rem',
          fontStyle: 'italic',
        }}
      >
        💡 {parameter.helpText}
      </div>
    );
  };

  // ===== 渲染验证错误 =====
  const renderValidation = () => {
    const errors: string[] = [];

    if (parameter.required && (localValue === undefined || localValue === null || localValue === '')) {
      errors.push('This field is required');
    }

    if (parameter.type === 'number' || parameter.type === 'integer') {
      const numValue = Number(localValue);
      if (parameter.min !== undefined && numValue < parameter.min) {
        errors.push(`Must be >= ${parameter.min}`);
      }
      if (parameter.max !== undefined && numValue > parameter.max) {
        errors.push(`Must be <= ${parameter.max}`);
      }
    }

    if (parameter.type === 'string' && typeof localValue === 'string') {
      if (parameter.minLength && localValue.length < parameter.minLength) {
        errors.push(`Must be at least ${parameter.minLength} characters`);
      }
      if (parameter.maxLength && localValue.length > parameter.maxLength) {
        errors.push(`Must be at most ${parameter.maxLength} characters`);
      }
      if (parameter.pattern && !new RegExp(parameter.pattern).test(localValue)) {
        errors.push('Invalid format');
      }
    }

    if (errors.length === 0) return null;

    return (
      <div style={{ marginTop: '0.25rem' }}>
        {errors.map((error, idx) => (
          <div key={idx} style={{ fontSize: '0.75rem', color: '#ef4444' }}>
            ⚠️ {error}
          </div>
        ))}
      </div>
    );
  };

  // ===== 根据类型渲染输入控件 =====
  const renderInput = () => {
    const baseStyle: React.CSSProperties = {
      width: '100%',
      padding: '0.5rem',
      backgroundColor: '#111827',
      border: '1px solid #374151',
      borderRadius: '4px',
      color: '#f3f4f6',
      fontSize: '0.875rem',
    };

    switch (parameter.type) {
      case 'string':
        return (
          <input
            type="text"
            value={localValue || ''}
            onChange={(e) => handleChange(e.target.value)}
            placeholder={parameter.placeholder}
            disabled={disabled}
            style={baseStyle}
          />
        );

      case 'number':
      case 'integer':
        return (
          <input
            type="number"
            value={localValue ?? ''}
            onChange={(e) => {
              const val = parameter.type === 'integer' ? parseInt(e.target.value) : parseFloat(e.target.value);
              handleChange(isNaN(val) ? undefined : val);
            }}
            min={parameter.min}
            max={parameter.max}
            step={parameter.type === 'integer' ? 1 : 'any'}
            placeholder={parameter.placeholder}
            disabled={disabled}
            style={baseStyle}
          />
        );

      case 'boolean':
        return (
          <label style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', cursor: 'pointer' }}>
            <input
              type="checkbox"
              checked={localValue || false}
              onChange={(e) => handleChange(e.target.checked)}
              disabled={disabled}
              style={{ width: '1rem', height: '1rem', cursor: 'pointer' }}
            />
            <span style={{ fontSize: '0.875rem', color: '#e5e7eb' }}>
              {localValue ? 'Enabled' : 'Disabled'}
            </span>
          </label>
        );

      case 'select':
        return (
          <select
            value={localValue || ''}
            onChange={(e) => handleChange(e.target.value)}
            disabled={disabled}
            style={baseStyle}
          >
            <option value="">-- Select {parameter.label} --</option>
            {parameter.enum?.map((option: string) => (
              <option key={option} value={option}>
                {option}
              </option>
            ))}
          </select>
        );

      case 'multiselect':
        return (
          <select
            multiple
            value={localValue || []}
            onChange={(e) => {
              const selected = Array.from(e.target.selectedOptions, (option) => option.value);
              handleChange(selected);
            }}
            disabled={disabled}
            style={{ ...baseStyle, minHeight: '100px' }}
          >
            {parameter.enum?.map((option: string) => (
              <option key={option} value={option}>
                {option}
              </option>
            ))}
          </select>
        );

      case 'file':
      case 'directory':
        return (
          <input
            type="text"
            value={localValue || ''}
            onChange={(e) => handleChange(e.target.value)}
            placeholder={parameter.placeholder || `Path to ${parameter.type}`}
            disabled={disabled}
            style={baseStyle}
          />
        );

      case 'color':
        return (
          <div style={{ display: 'flex', gap: '0.5rem', alignItems: 'center' }}>
            <input
              type="color"
              value={localValue || '#000000'}
              onChange={(e) => handleChange(e.target.value)}
              disabled={disabled}
              style={{ width: '3rem', height: '2.5rem', cursor: 'pointer', border: '1px solid #374151' }}
            />
            <input
              type="text"
              value={localValue || ''}
              onChange={(e) => handleChange(e.target.value)}
              placeholder="#000000"
              disabled={disabled}
              style={{ ...baseStyle, flex: 1 }}
            />
          </div>
        );

      case 'date':
        return (
          <input
            type="date"
            value={localValue || ''}
            onChange={(e) => handleChange(e.target.value)}
            disabled={disabled}
            style={baseStyle}
          />
        );

      case 'range':
        return (
          <div>
            <input
              type="range"
              value={localValue ?? parameter.default ?? parameter.min ?? 0}
              onChange={(e) => handleChange(parseFloat(e.target.value))}
              min={parameter.min}
              max={parameter.max}
              step={(parameter.max && parameter.min) ? (parameter.max - parameter.min) / 100 : 1}
              disabled={disabled}
              style={{ width: '100%', cursor: 'pointer' }}
            />
            <div style={{ textAlign: 'center', fontSize: '0.875rem', color: '#9ca3af', marginTop: '0.25rem' }}>
              {localValue ?? parameter.default ?? parameter.min ?? 0}
            </div>
          </div>
        );

      case 'textarea':
        return (
          <textarea
            value={localValue || ''}
            onChange={(e) => handleChange(e.target.value)}
            placeholder={parameter.placeholder}
            disabled={disabled}
            rows={5}
            style={{ ...baseStyle, resize: 'vertical', fontFamily: 'monospace' }}
          />
        );

      case 'json':
        return (
          <textarea
            value={typeof localValue === 'string' ? localValue : JSON.stringify(localValue, null, 2)}
            onChange={(e) => {
              try {
                const parsed = JSON.parse(e.target.value);
                handleChange(parsed);
              } catch {
                handleChange(e.target.value);
              }
            }}
            placeholder={parameter.placeholder || '{}'}
            disabled={disabled}
            rows={8}
            style={{ ...baseStyle, resize: 'vertical', fontFamily: 'monospace' }}
          />
        );

      case 'array':
        const arrayValue = Array.isArray(localValue) ? localValue : [];
        return (
          <div>
            {arrayValue.map((item: any, idx: number) => (
              <div key={idx} style={{ display: 'flex', gap: '0.5rem', marginBottom: '0.5rem' }}>
                <input
                  type="text"
                  value={item}
                  onChange={(e) => {
                    const newArray = [...arrayValue];
                    newArray[idx] = e.target.value;
                    handleChange(newArray);
                  }}
                  disabled={disabled}
                  style={{ ...baseStyle, flex: 1 }}
                />
                <button
                  onClick={() => {
                    const newArray = arrayValue.filter((_, i) => i !== idx);
                    handleChange(newArray);
                  }}
                  disabled={disabled}
                  style={{
                    padding: '0.5rem',
                    backgroundColor: '#dc2626',
                    color: 'white',
                    border: 'none',
                    borderRadius: '4px',
                    cursor: 'pointer',
                  }}
                >
                  ✕
                </button>
              </div>
            ))}
            <button
              onClick={() => handleChange([...arrayValue, ''])}
              disabled={disabled}
              style={{
                padding: '0.5rem 1rem',
                backgroundColor: '#3b82f6',
                color: 'white',
                border: 'none',
                borderRadius: '4px',
                cursor: 'pointer',
                fontSize: '0.875rem',
              }}
            >
              + Add Item
            </button>
          </div>
        );

      default:
        return (
          <div style={{ padding: '0.5rem', backgroundColor: '#7f1d1d', borderRadius: '4px', color: '#fecaca' }}>
            Unsupported parameter type: {parameter.type}
          </div>
        );
    }
  };

  return (
    <div style={{ marginBottom: '1rem' }}>
      {renderLabel()}
      {renderInput()}
      {renderValidation()}
      {renderHelp()}
    </div>
  );
};

export default DynamicParameterRenderer;
