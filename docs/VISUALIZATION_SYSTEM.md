# 可视化系统实现总结

## 概述

本文档总结了 Omicsomics 平台中完整的可视化系统实现，包括内置数据可视化工具、工具输出预览、Run 执行监控和结果详情页面。

## 1. 内置数据可视化工具节点（TODO #9）

### 新增可视化工具

在 `ToolDefinitions.ts` 中添加了 8 个专业的数据可视化工具节点：

#### 1.1 已有工具（扩展前）

- **Heatmap（热图）**: 层次聚类、颜色方案选择
- **Volcano Plot（火山图）**: 差异分析可视化
- **PCA Plot（主成分分析图）**: 降维可视化
- **Box Plot（箱线图）**: 数据分布展示

#### 1.2 新增工具（本次扩展）

- **Violin Plot（小提琴图）**

  - 参数：显示内部箱线图、核密度带宽
  - 用途：结合密度和分布的可视化

- **MA Plot（MA 图）**

  - 参数：P 值阈值、倍数变化阈值
  - 用途：差异表达分析（M vs A）

- **Scatter Plot（散点图）**

  - 参数：X/Y 轴列选择、回归线、颜色分组
  - 用途：相关性分析

- **Histogram（直方图）**

  - 参数：分箱数、核密度估计
  - 用途：数据频率分布

- **Correlation Matrix（相关性矩阵）**

  - 参数：相关性方法（Pearson/Spearman/Kendall）、显示数值
  - 用途：多变量相关性热图

- **Dendrogram（树状图）**

  - 参数：连接方法、距离度量
  - 用途：层次聚类可视化

- **Bar Plot（条形图）**
  - 参数：方向、误差条
  - 用途：分类数据比较

### 工具特性

- **类型**: `type: 'visualization'`
- **图标**: 📊（在 ToolboxV2 中显示为 📈）
- **分类**: `operationCategory: 'Visualization'`
- **通用性**: `omicsType: 'General'`（适用于所有组学类型）
- **可配置**: 每个工具都有详细的参数模板

## 2. 工具输出预览系统（TODO #10）

### 组件：ToolOutputPreview.tsx

#### 功能特性

- **多文件支持**: 可切换查看多个输出文件
- **文件类型支持**:
  - **Text**: 代码/日志文件预览（语法高亮）
  - **Image**: 图片预览（支持缩放）
  - **Table**: 表格数据展示（带分页）
  - **JSON**: JSON 数据美化显示

#### UI 特性

- **文件信息面板**:
  - 文件名、类型、大小、修改时间
  - 文件选择器（多文件时）
- **内容渲染**:

  - 文本：monospace 字体，深色背景
  - 表格：斑马纹行、粘性表头
  - 图片：居中显示、自适应大小
  - JSON：格式化缩进显示

- **操作按钮**:
  - 💾 Download：下载文件
  - 📊 Export to Data：导出到数据集

#### 使用场景

- 在节点配置面板中点击"输出预览"
- 查看工具执行后的结果文件
- 快速验证数据处理效果

## 3. Run 执行状态可视化（TODO #11）

### 组件：RunExecutionVisualizer.tsx

#### 核心功能

##### 3.1 Pipeline 拓扑可视化

- **ReactFlow 集成**: 交互式流程图
- **实时状态高亮**:

  - ⏳ Waiting（灰色 #374151）
  - ▶️ Running（蓝色 #3b82f6，带脉冲动画）
  - ✅ Completed（绿色 #10b981）
  - ❌ Failed（红色 #dc2626）

- **边缘动画**: 当前执行路径带动画效果
- **节点点击**: 点击节点查看详细信息

##### 3.2 总体统计面板

- **进度条**: 实时显示整体完成百分比
- **状态计数器**:
  - Total（总节点数）
  - Waiting（等待中）
  - Running（执行中）
  - Completed（已完成）
  - Failed（失败）

##### 3.3 工具执行详情列表

- **每个工具显示**:
  - 工具名称和状态徽章
  - 开始时间和执行时长
  - CPU 使用率和内存占用
  - 执行进度百分比（运行中时）
  - 错误信息（失败时）

#### 性能监控

- **资源追踪**:

  - CPU 使用率（%）
  - 内存使用量（MB）
  - 执行时长（秒）

- **实时更新**: 支持 WebSocket 实时更新状态

## 4. Run 结果详情页面（TODO #12）

### 组件：RunDetailsPage.tsx

#### 页面结构

##### 4.1 Header 区域

- **Run 信息**:
  - Pipeline 名称
  - Run ID
  - 状态徽章（Running/Completed/Failed）
- **总体指标**:
  - 开始/结束时间
  - 总执行时长
  - 平均 CPU 使用率
  - 总内存使用量

##### 4.2 四个视图标签

###### Overview（概览）

- 嵌入 `RunExecutionVisualizer` 组件
- Pipeline 拓扑图 + 实时状态
- 工具执行详情列表

###### Timeline（时间线）

- **垂直时间轴设计**:
  - 左侧时间线（垂直线 + 状态点）
  - 每个工具为一个时间点
  - 状态点颜色对应执行状态
- **时间点信息**:
  - 工具名称和类型
  - 状态徽章
  - 开始时间和执行时长
  - 输入/输出文件数量
  - 点击展开详细日志

###### Logs（日志）

- **全局日志视图**: 所有工具的日志列表
- **工具筛选**: 点击工具名查看单个工具日志
- **日志查看器**:
  - Monospace 字体
  - 深色主题
  - 滚动查看
  - 日志项数统计

###### Files（文件）

- **输入文件列表**: 每个工具的输入文件
- **输出文件列表**: 每个工具的输出文件
- **文件信息**:
  - 文件名
  - 文件类型
  - 文件大小
- **点击预览**: 点击输出文件打开 `ToolOutputPreview`

#### 交互功能

- **工具选择**: 点击时间线或概览中的工具
- **文件预览**: 点击输出文件查看内容
- **日志过滤**: 切换显示全部/单个工具日志
- **下载功能**: 下载输出文件（TODO）

## 5. 集成点

### 5.1 Pipeline Editor 集成

- 可视化工具节点可从工具库拖拽到画布
- 支持所有标准节点功能（配置、连接、删除）
- 参数配置面板显示可视化特定选项

### 5.2 RUNS 页面集成

需要在现有 RUNS 页面中添加：

```typescript
// 示例：集成到 RunsPage.tsx
import { useNavigate } from 'react-router-dom';

const navigate = useNavigate();

// 在 run 列表中添加"查看详情"按钮
<button onClick={() => navigate(`/runs/${run.id}`)}>
  View Details
</button>

// 路由配置
// App.tsx 或 routes.tsx
<Route path="/runs/:runId" element={<RunDetailsPage />} />
```

### 5.3 ConfigPanel 集成

在 `ConfigPanelV2.tsx` 中添加"输出预览"标签：

```typescript
// 添加预览按钮
{
  selectedNode.type === "visualization" && (
    <button onClick={() => setShowOutputPreview(true)}>
      📊 Preview Output
    </button>
  );
}

<ToolOutputPreview
  isOpen={showOutputPreview}
  onClose={() => setShowOutputPreview(false)}
  toolName={selectedNode.data.label}
  outputFiles={mockOutputFiles}
/>;
```

## 6. 数据流

### 6.1 Pipeline 执行流程

```
Pipeline Submit
    ↓
Backend 创建 Run
    ↓
WebSocket 连接建立
    ↓
工具逐个执行
    ↓
实时状态更新 → RunExecutionVisualizer
    ↓
输出文件生成 → ToolOutputPreview
    ↓
Run 完成 → RunDetailsPage
```

### 6.2 状态更新机制

```typescript
// 后端推送更新（WebSocket）
{
  "type": "tool_status_update",
  "runId": "run-123",
  "toolId": "tool-456",
  "status": "running",
  "progress": 45,
  "cpuUsage": 78.5,
  "memoryUsage": 1024
}

// 前端处理
useEffect(() => {
  const ws = new WebSocket(`ws://api/runs/${runId}/stream`);
  ws.onmessage = (event) => {
    const update = JSON.parse(event.data);
    updateToolExecution(update);
  };
}, [runId]);
```

## 7. 样式系统

### 颜色方案（深色主题）

```typescript
const theme = {
  background: {
    primary: "#0f172a", // 主背景
    secondary: "#1f2937", // 卡片背景
    tertiary: "#111827", // 深层背景
  },
  text: {
    primary: "#f3f4f6", // 主文本
    secondary: "#e5e7eb", // 次要文本
    tertiary: "#9ca3af", // 三级文本
  },
  border: "#374151",
  status: {
    waiting: "#6b7280",
    running: "#3b82f6",
    completed: "#10b981",
    failed: "#dc2626",
  },
};
```

### 响应式设计

- 使用 `display: grid` 和 `minmax()` 实现自适应布局
- 移动端优化（待完善）

## 8. 性能优化

### 8.1 大数据处理

- 表格分页显示（默认 20 行/页）
- 日志懒加载
- 图片延迟加载

### 8.2 实时更新优化

- WebSocket 连接复用
- 状态更新节流（避免频繁重渲染）
- 仅更新变化的节点

## 9. 待完善功能

### 9.1 输出预览

- [ ] 实际文件下载功能
- [ ] 更多文件类型支持（PDF、HDF5）
- [ ] 图片缩放和平移
- [ ] 大文件流式加载

### 9.2 Run 监控

- [ ] WebSocket 实时连接
- [ ] 性能图表（时间序列）
- [ ] 资源使用趋势图
- [ ] 多 Run 对比

### 9.3 可视化工具

- [ ] 交互式参数调整
- [ ] 实时预览（小样本）
- [ ] 导出高分辨率图片
- [ ] 自定义配色方案

## 10. 测试建议

### 10.1 单元测试

```typescript
// ToolOutputPreview.test.tsx
describe("ToolOutputPreview", () => {
  it("renders text files correctly", () => {
    // 测试文本文件渲染
  });

  it("switches between multiple files", () => {
    // 测试多文件切换
  });
});
```

### 10.2 集成测试

- 完整 Pipeline 执行流程
- 状态更新正确性
- 文件预览功能

### 10.3 性能测试

- 大 Pipeline（100+ 节点）渲染
- 1000+ 行日志显示
- 10+ MB 文件预览

## 11. 文档和示例

### 11.1 用户文档

需要创建用户指南：

- 如何使用可视化工具节点
- 如何查看 Run 执行状态
- 如何分析输出文件

### 11.2 开发者文档

- 如何添加新的可视化工具
- 如何扩展文件类型支持
- WebSocket 协议规范

## 12. 总结

### 已完成

✅ 8 个新的可视化工具节点  
✅ ToolOutputPreview 组件（支持 4 种文件类型）  
✅ RunExecutionVisualizer 组件（实时状态 + 性能监控）  
✅ RunDetailsPage 组件（4 个视图标签）  
✅ 深色主题 UI 设计  
✅ 响应式布局

### 进行中

🔄 后端 API 集成  
🔄 WebSocket 实时更新  
🔄 实际文件处理

### 待开始

⏳ 性能优化  
⏳ 移动端适配  
⏳ 用户文档  
⏳ 测试覆盖

---

**创建时间**: 2025-01-10  
**最后更新**: 2025-01-10  
**版本**: 1.0
