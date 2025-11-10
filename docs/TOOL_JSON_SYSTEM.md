# 通用工具 JSON 配置系统

## 概述

该系统实现了一个完全可扩展的、基于 JSON Schema 的工具配置框架，允许用户上传自定义工具定义并自动生成图形化参数配置界面。

## 核心特性

### 1. 标准化 JSON Schema (`ToolSchema.ts`)

定义了统一的工具定义格式，包含：

**元数据**

- `$schema`: Schema 版本号
- `id`: 唯一标识符（建议格式：`author.toolname`）
- `name`: 工具名称
- `version`: 版本号
- `author`, `license`, `homepage`, `repository`: 可选元信息

**分类**

- `omicsType`: 组学类型（genomics, transcriptomics 等）
- `category`: 操作类别（qc, alignment, peak_calling 等）
- `tags`: 标签数组（用于搜索）

**参数定义** (`parameters`)

- 支持 15 种参数类型：
  - 基础类型：`string`, `number`, `integer`, `boolean`
  - 选择类型：`select`, `multiselect`
  - 文件类型：`file`, `directory`
  - 特殊类型：`color`, `date`, `range`, `textarea`, `json`, `array`
- 验证规则：`min`, `max`, `minLength`, `maxLength`, `pattern`, `enum`
- UI 提示：`placeholder`, `helpText`, `unit`, `group`
- 条件显示：`dependsOn`

**输入输出**

- `inputs`: 定义输入数据类型和要求
- `outputs`: 定义输出数据类型

**执行配置**

- `execution.container`: Docker 镜像
- `execution.command`: 命令模板
- `execution.environment`: 环境变量
- `execution.resources`: 资源限制（CPU、内存、GPU）

**UI 配置**

- `ui.icon`: 图标
- `ui.color`: 主题颜色
- `ui.layout`: 布局样式
- `ui.parameterGroups`: 参数分组配置

### 2. 工具注册中心 (`toolRegistry.ts`)

核心功能：

- **加载内置工具**：自动从旧的`ToolDefinitions.ts`转换
- **管理自定义工具**：localStorage 持久化
- **启用/禁用工具**：灵活控制工具可见性
- **导入/导出**：
  - 导入单个工具 JSON
  - 导入工具集合
  - 导出为 JSON 文件
  - 批量导出所有自定义工具
- **搜索和过滤**：
  - 按组学类型
  - 按操作类别
  - 关键字搜索（名称、描述、ID、标签）
- **验证**：自动验证工具定义的完整性和正确性

API 示例：

```typescript
import { toolRegistry } from "../utils/toolRegistry";

// 获取所有工具
const tools = toolRegistry.getAllTools();

// 按类型获取
const genomicsTools = toolRegistry.getToolsByOmicsType("genomics");

// 搜索
const results = toolRegistry.searchTools("trimmer");

// 添加自定义工具
const result = toolRegistry.addCustomTool(toolDefinition);

// 导入JSON文件
const importResult = await toolRegistry.importToolFromJSON(file);

// 导出工具
toolRegistry.downloadToolJSON("custom.mytool");
```

### 3. 动态参数渲染器 (`DynamicParameterRenderer.tsx`)

根据 JSON Schema 自动生成 UI 组件：

**支持的参数类型渲染**：

- `string` → 文本输入框
- `number`/`integer` → 数字输入框（带 min/max 验证）
- `boolean` → 复选框
- `select` → 下拉选择框
- `multiselect` → 多选下拉框
- `file`/`directory` → 文件路径输入
- `color` → 颜色选择器 + 文本输入
- `date` → 日期选择器
- `range` → 滑块（带实时值显示）
- `textarea` → 多行文本编辑器
- `json` → JSON 编辑器（带验证）
- `array` → 动态列表（可添加/删除项）

**自动验证**：

- 必需字段检查
- 数值范围验证
- 字符串长度验证
- 正则表达式匹配

**UI 特性**：

- 显示参数标签、描述、单位
- 必需参数标记（\*）
- 实时验证错误提示
- 帮助文本提示（💡）

### 4. 工具管理页面 (`ToolManagerPage.tsx`)

完整的工具库管理界面：

**统计面板**

- 总工具数
- 内置工具数
- 自定义工具数
- 已启用工具数

**操作功能**

- 📤 上传工具 JSON
- ➕ 创建新工具（使用模板）
- 💾 导出所有自定义工具

**筛选和搜索**

- 🔍 关键字搜索
- 按组学类型筛选
- 按操作类别筛选
- 仅显示自定义工具

**工具卡片**

- 显示工具名称、版本、分类
- 参数数量统计
- 操作按钮：
  - View: 查看工具 JSON 定义
  - Export: 导出单个工具
  - Enable/Disable: 启用/禁用工具
  - Delete: 删除自定义工具（内置工具不可删除）

**访问路径**: `/tools`

### 5. 通用配置面板 (`ConfigPanelV3.tsx`)

升级版参数配置面板：

**特性**

- 从`toolRegistry`动态加载工具定义
- 使用`DynamicParameterRenderer`渲染所有参数
- 支持参数分组（可折叠）
- 参数预设管理：
  - 加载预设
  - 保存当前配置为预设
  - 删除预设
  - 重置为默认值
- 工具信息卡片：
  - 名称、版本、组学类型、类别
  - 描述和文档链接

**参数分组**

- 未分组参数显示在顶部
- 按`group`字段自动分组
- 可折叠的分组面板

### 6. 工具 JSON 编辑器 (`ToolJSONEditor.tsx`)

在线编辑工具定义：

**功能**

- ✓ Valid / ⚠️ Validate: 实时验证
- 🎨 Format: 格式化 JSON
- 📄 Load Template: 加载示例模板
- 💾 Save Tool: 保存到工具库

**编辑器特性**

- Monaco 风格等宽字体
- 语法高亮（浏览器原生）
- 验证错误提示面板
- 内置 Schema 帮助文档

**用例**

- 创建新工具
- 编辑现有自定义工具
- 学习 JSON Schema 格式

## 使用流程

### 创建自定义工具

#### 方式 1：上传 JSON 文件

1. 访问 `/tools`
2. 点击 "📤 Upload Tool JSON"
3. 选择你的工具 JSON 文件
4. 系统自动验证并导入

#### 方式 2：在线创建

1. 访问 `/tools`
2. 点击 "➕ Create New Tool"
3. 在 JSON 编辑器中修改模板
4. 点击 "Validate" 验证
5. 点击 "Save Tool" 保存

#### 方式 3：从现有工具复制

1. 在工具列表中找到相似工具
2. 点击 "Export" 导出 JSON
3. 编辑 JSON 文件
4. 重新上传

### 使用自定义工具

1. 在 `/custom-pipelines` 创建管道
2. 从工具箱中拖拽你的自定义工具
3. 点击节点配置参数
4. 系统自动生成配置 UI（基于 JSON Schema）

### 分享工具

**导出单个工具**：

1. 在 `/tools` 找到工具
2. 点击 "Export"
3. 分享 JSON 文件给其他用户

**导出工具集合**：

1. 点击 "💾 Export All Custom Tools"
2. 获得包含所有自定义工具的 JSON 文件
3. 其他用户上传即可导入所有工具

## JSON Schema 示例

```json
{
  "$schema": "1.0.0",
  "id": "community.custom_trimmer",
  "name": "Custom Trimmer",
  "version": "1.0.0",
  "author": "John Doe",
  "license": "MIT",
  "omicsType": "genomics",
  "category": "qc",
  "description": "A custom quality trimming tool",
  "parameters": [
    {
      "name": "quality_threshold",
      "label": "Quality Threshold",
      "type": "integer",
      "description": "Minimum quality score",
      "required": true,
      "default": 20,
      "min": 0,
      "max": 40,
      "unit": "Q",
      "group": "Quality Control"
    },
    {
      "name": "output_format",
      "label": "Output Format",
      "type": "select",
      "required": true,
      "default": "fastq",
      "enum": ["fastq", "fasta", "fastq.gz"],
      "group": "Output"
    }
  ],
  "inputs": [
    {
      "type": "fastq",
      "label": "Input FASTQ",
      "required": true
    }
  ],
  "outputs": [
    {
      "type": "fastq",
      "label": "Trimmed FASTQ"
    }
  ]
}
```

## 存储

- **自定义工具**: `localStorage` → `omicsomics_custom_tools`
- **禁用工具列表**: `localStorage` → `omicsomics_disabled_tools`
- **参数预设**: `localStorage` → `pipeline_parameter_presets`

数据在浏览器本地持久化，可通过导出功能备份。

## 扩展性

### 添加新参数类型

1. 在`ToolSchema.ts`中添加到`ParameterType`枚举
2. 在`DynamicParameterRenderer.tsx`中添加渲染逻辑
3. 更新 JSON Schema 文档

### 添加验证规则

1. 在`ToolParameterSchema`接口添加新字段
2. 在`validateToolSchema()`函数添加验证逻辑
3. 在`DynamicParameterRenderer`添加 UI 错误提示

### 自定义 UI 组件

继承`DynamicParameterRenderer`或创建自定义渲染器：

```typescript
const CustomRenderer = ({ parameter, value, onChange }) => {
  // 自定义渲染逻辑
};
```

## 技术栈

- **TypeScript**: 类型安全
- **React**: UI 组件
- **localStorage**: 客户端持久化
- **JSON Schema**: 标准化定义格式
- **React Flow**: 管道编辑器集成

## 优势

1. **完全可扩展**：用户可添加任意工具，无需修改代码
2. **自动 UI 生成**：基于 JSON 自动生成配置界面
3. **标准化格式**：统一的 JSON Schema 便于分享和集成
4. **验证保障**：自动验证确保工具定义正确性
5. **社区友好**：支持导入导出，方便社区共享工具
6. **向后兼容**：自动转换旧的工具定义格式

## 未来扩展

- [ ] 工具市场：从云端浏览和安装社区工具
- [ ] 版本管理：工具定义版本控制和更新
- [ ] 依赖管理：工具间依赖关系自动解析
- [ ] 在线文档：自动生成工具使用文档
- [ ] 性能优化：虚拟列表、懒加载
- [ ] 高级编辑器：Monaco Editor 集成，语法高亮
- [ ] 工具测试：参数验证和 dry-run 测试
- [ ] 国际化：多语言工具定义支持

## 文件清单

```
frontend/src/
├── schemas/
│   └── ToolSchema.ts                    # JSON Schema定义和验证
├── utils/
│   └── toolRegistry.ts                  # 工具注册中心
├── components/
│   ├── ToolJSONEditor.tsx               # JSON编辑器
│   └── pipelines/
│       ├── DynamicParameterRenderer.tsx # 动态参数渲染器
│       └── ConfigPanelV3.tsx            # 通用配置面板
└── pages/
    └── ToolManagerPage.tsx              # 工具管理页面
```

## 总结

这个系统彻底改变了工具配置的方式，从硬编码到完全数据驱动。用户现在可以：

1. **上传自己的工具**：无需编程即可扩展系统
2. **自动生成 UI**：系统根据 JSON 自动创建配置界面
3. **分享工具定义**：JSON 格式便于分享和协作
4. **统一管理**：所有工具在一个地方管理
5. **灵活配置**：支持复杂的参数类型和验证规则

这是一个真正的**可扩展生物信息学工具平台**！ 🚀
