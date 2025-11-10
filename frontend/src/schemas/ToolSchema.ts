/**
 * Tool JSON Schema Definition
 * 标准化的工具定义格式，支持用户自定义工具上传和动态UI生成
 * 
 * 这是一个完全可扩展的JSON格式，用户可以：
 * 1. 上传自己的工具JSON文件
 * 2. 在线编辑工具定义
 * 3. 从社区导入工具
 * 4. 导出工具配置为JSON
 */

// ===== JSON Schema Version =====
export const TOOL_SCHEMA_VERSION = '1.0.0';

// ===== Parameter Types =====
export type ParameterType = 
  | 'string'      // 文本输入
  | 'number'      // 数值输入
  | 'integer'     // 整数输入
  | 'boolean'     // 布尔开关
  | 'select'      // 下拉选择
  | 'multiselect' // 多选下拉
  | 'file'        // 文件路径
  | 'directory'   // 目录路径
  | 'color'       // 颜色选择器
  | 'date'        // 日期选择
  | 'range'       // 滑块范围
  | 'textarea'    // 多行文本
  | 'json'        // JSON编辑器
  | 'array';      // 数组（可添加多个值）

// ===== Parameter Definition =====
export interface ToolParameterSchema {
  // 基础属性
  name: string;                    // 参数名称（用于后端传递）
  label: string;                   // 显示标签
  type: ParameterType;             // 参数类型
  description?: string;            // 参数描述
  required?: boolean;              // 是否必需
  
  // 默认值
  default?: any;                   // 默认值
  
  // 验证规则
  min?: number;                    // 最小值（number类型）
  max?: number;                    // 最大值（number类型）
  minLength?: number;              // 最小长度（string类型）
  maxLength?: number;              // 最大长度（string类型）
  pattern?: string;                // 正则表达式（string类型）
  enum?: string[];                 // 枚举值（select类型）
  
  // UI提示
  placeholder?: string;            // 占位符文本
  helpText?: string;               // 帮助文本
  unit?: string;                   // 单位（如 'bp', 'GB', 'threads'）
  
  // 高级配置
  dependsOn?: {                    // 条件显示（依赖其他参数）
    parameter: string;             // 依赖的参数名
    value: any;                    // 触发显示的值
  };
  group?: string;                  // 参数分组（用于UI折叠）
  advanced?: boolean;              // 是否为高级参数
  
  // 文件类型特定
  accept?: string[];               // 接受的文件扩展名（file类型）
  multiple?: boolean;              // 是否允许多文件（file类型）
}

// ===== Tool Definition Schema =====
export interface ToolDefinitionSchema {
  // 元数据
  $schema: string;                 // Schema版本（固定为TOOL_SCHEMA_VERSION）
  id: string;                      // 唯一标识符（建议: author.toolname）
  name: string;                    // 工具名称
  version: string;                 // 工具版本
  author?: string;                 // 作者
  license?: string;                // 许可证
  homepage?: string;               // 主页链接
  repository?: string;             // 代码仓库
  
  // 分类
  omicsType: string;               // 组学类型（genomics, transcriptomics, etc.）
  category: string;                // 操作类别（qc, alignment, etc.）
  tags?: string[];                 // 标签（用于搜索）
  
  // 描述
  description: string;             // 简短描述
  longDescription?: string;        // 详细描述（支持Markdown）
  documentation?: string;          // 文档链接
  citation?: string;               // 引用信息
  
  // 参数定义
  parameters: ToolParameterSchema[]; // 参数列表
  
  // 输入输出
  inputs?: {
    type: string;                  // 输入类型（fastq, bam, vcf, etc.）
    label: string;                 // 输入标签
    required: boolean;             // 是否必需
    multiple?: boolean;            // 是否支持多个输入
  }[];
  
  outputs?: {
    type: string;                  // 输出类型
    label: string;                 // 输出标签
    description?: string;          // 输出描述
  }[];
  
  // 执行配置
  execution?: {
    container?: string;            // Docker镜像
    command?: string;              // 命令模板（支持变量替换）
    environment?: Record<string, string>; // 环境变量
    resources?: {
      cpu?: number;                // CPU核心数
      memory?: string;             // 内存限制（如 "4GB"）
      gpu?: boolean;               // 是否需要GPU
    };
  };
  
  // UI配置
  ui?: {
    icon?: string;                 // 图标名称或URL
    color?: string;                // 主题颜色
    layout?: 'compact' | 'detailed'; // UI布局样式
    parameterGroups?: {            // 参数分组配置
      name: string;
      label: string;
      collapsed?: boolean;
    }[];
  };
  
  // 扩展字段（用户可以添加自定义元数据）
  [key: string]: any;
}

// ===== Tool Collection Schema (用于批量导入) =====
export interface ToolCollectionSchema {
  $schema: string;
  name: string;
  description: string;
  version: string;
  author?: string;
  tools: ToolDefinitionSchema[];
}

// ===== 验证函数 =====
export function validateToolSchema(tool: any): { valid: boolean; errors: string[] } {
  const errors: string[] = [];
  
  // 必需字段检查
  if (!tool.$schema) errors.push('Missing $schema field');
  if (!tool.id) errors.push('Missing id field');
  if (!tool.name) errors.push('Missing name field');
  if (!tool.version) errors.push('Missing version field');
  if (!tool.omicsType) errors.push('Missing omicsType field');
  if (!tool.category) errors.push('Missing category field');
  if (!tool.description) errors.push('Missing description field');
  if (!tool.parameters || !Array.isArray(tool.parameters)) {
    errors.push('Missing or invalid parameters field');
  }
  
  // 参数验证
  if (tool.parameters && Array.isArray(tool.parameters)) {
    tool.parameters.forEach((param: any, idx: number) => {
      if (!param.name) errors.push(`Parameter ${idx}: missing name`);
      if (!param.label) errors.push(`Parameter ${idx}: missing label`);
      if (!param.type) errors.push(`Parameter ${idx}: missing type`);
      
      // 类型特定验证
      if (param.type === 'select' && !param.enum) {
        errors.push(`Parameter ${param.name}: select type requires enum field`);
      }
    });
  }
  
  // ID格式检查（建议格式: author.toolname）
  if (tool.id && !/^[a-zA-Z0-9._-]+$/.test(tool.id)) {
    errors.push('Tool ID should only contain alphanumeric characters, dots, hyphens, and underscores');
  }
  
  return {
    valid: errors.length === 0,
    errors,
  };
}

// ===== 从旧格式转换 =====
export function convertLegacyToolToSchema(legacyTool: any): ToolDefinitionSchema {
  return {
    $schema: TOOL_SCHEMA_VERSION,
    id: legacyTool.id || `legacy.${legacyTool.name.toLowerCase().replace(/\s+/g, '-')}`,
    name: legacyTool.name,
    version: legacyTool.version || '1.0.0',
    omicsType: legacyTool.omicsType || 'general',
    category: legacyTool.operationCategory || 'general',
    description: legacyTool.description || '',
    parameters: (legacyTool.parameterTemplate || []).map((param: any) => ({
      name: param.name || param.label?.toLowerCase().replace(/\s+/g, '_'),
      label: param.label,
      type: param.type,
      description: param.description,
      required: param.required,
      default: param.default,
      min: param.min,
      max: param.max,
      enum: param.options,
      unit: param.unit,
    })),
    documentation: legacyTool.documentation,
  };
}

// ===== JSON Schema示例 =====
export const TOOL_SCHEMA_EXAMPLE: ToolDefinitionSchema = {
  $schema: TOOL_SCHEMA_VERSION,
  id: 'community.custom_trimmer',
  name: 'Custom Trimmer',
  version: '1.0.0',
  author: 'John Doe',
  license: 'MIT',
  omicsType: 'genomics',
  category: 'qc',
  description: 'A custom quality trimming tool',
  longDescription: `
# Custom Trimmer

This tool performs quality trimming on FASTQ files with customizable parameters.

## Features
- Quality score filtering
- Adapter trimming
- Length filtering
`,
  parameters: [
    {
      name: 'quality_threshold',
      label: 'Quality Threshold',
      type: 'integer',
      description: 'Minimum quality score',
      required: true,
      default: 20,
      min: 0,
      max: 40,
      unit: 'Q',
      group: 'Quality Control',
    },
    {
      name: 'min_length',
      label: 'Minimum Length',
      type: 'integer',
      description: 'Minimum read length after trimming',
      required: true,
      default: 50,
      min: 10,
      max: 1000,
      unit: 'bp',
      group: 'Length Filtering',
    },
    {
      name: 'adapter_sequence',
      label: 'Adapter Sequence',
      type: 'string',
      description: 'Adapter sequence to remove',
      required: false,
      placeholder: 'AGATCGGAAGAGC',
      group: 'Adapter Trimming',
    },
    {
      name: 'output_format',
      label: 'Output Format',
      type: 'select',
      description: 'Output file format',
      required: true,
      default: 'fastq',
      enum: ['fastq', 'fasta', 'fastq.gz'],
      group: 'Output',
    },
  ],
  inputs: [
    {
      type: 'fastq',
      label: 'Input FASTQ',
      required: true,
      multiple: false,
    },
  ],
  outputs: [
    {
      type: 'fastq',
      label: 'Trimmed FASTQ',
      description: 'Quality trimmed reads',
    },
  ],
  execution: {
    container: 'biocontainers/custom-trimmer:1.0.0',
    command: 'custom-trimmer -q {{quality_threshold}} -l {{min_length}} -a {{adapter_sequence}} -o {{output}}',
    resources: {
      cpu: 4,
      memory: '4GB',
    },
  },
  ui: {
    icon: 'scissors',
    color: '#3b82f6',
    layout: 'detailed',
    parameterGroups: [
      {
        name: 'quality_control',
        label: 'Quality Control',
        collapsed: false,
      },
      {
        name: 'adapter_trimming',
        label: 'Adapter Trimming',
        collapsed: true,
      },
    ],
  },
};
