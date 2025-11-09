# Frontend Data Processing UI - Manual Test Guide

## 测试日期

2025 年 1 月 9 日

## 前提条件

- ✅ 后端服务运行中 (http://localhost:8001)
- ✅ 前端服务运行中 (http://localhost:5173)
- ✅ 已有测试账号和项目

## 新功能清单

### 1. 上传时处理选项

**位置**: Data Catalog 页面 → Upload File 表单

**测试步骤**:

1. 点击 "+ Upload File" 按钮
2. 观察新增的复选框: "Process file immediately"
3. 验证默认勾选状态
4. 验证说明文字显示: "Convert to unified format for analysis..."

**预期结果**:

- ✅ 复选框默认勾选
- ✅ 可以取消勾选(仅上传原始文件)
- ✅ 说明文字清晰易懂

### 2. 文件处理状态显示

**位置**: Data Catalog 页面 → 文件列表

**测试步骤**:

1. 上传一个文件(勾选"Process file immediately")
2. 等待上传完成,刷新文件列表
3. 观察文件状态列

**预期结果**:

- ✅ 已处理文件显示绿色 "✓ Processed"
- ✅ 未处理文件显示灰色 "Raw only"
- ✅ 显示 omics type 标签(transcriptomics/genomics/proteomics 等)

### 3. 查看处理后数据

**位置**: Data Catalog 页面 → 已处理文件的 "👁️ View Data" 按钮

**测试步骤**:

1. 点击已处理文件的 "👁️ View Data" 按钮
2. 观察弹出的模态窗口
3. 检查显示的数据内容

**预期结果**:

- ✅ 模态窗口正确弹出
- ✅ 显示元数据(Omics Type, Source, Total Records 等)
- ✅ 显示前 10 条数据记录(表格形式)
- ✅ 显示处理信息(Processing Info JSON)
- ✅ 可以滚动查看完整内容
- ✅ "✕ Close" 按钮可关闭窗口

### 4. 多文件选择功能

**位置**: Data Catalog 页面 → 文件列表复选框

**测试步骤**:

1. 上传多个文件(都勾选处理)
2. 勾选文件列表中的复选框
3. 观察顶部的选择提示

**预期结果**:

- ✅ 只有已处理的文件可以勾选
- ✅ 未处理的文件复选框禁用
- ✅ 全选复选框可一次选中所有已处理文件
- ✅ 选中文件后显示蓝色提示条
- ✅ 提示条显示选中数量

### 5. 创建分析运行按钮

**位置**: Data Catalog 页面 → 选中文件后显示的操作栏

**测试步骤**:

1. 选中 1 个或多个已处理文件
2. 观察蓝色操作栏
3. 点击 "🚀 Create Analysis Run" 按钮

**预期结果**:

- ✅ 操作栏只在有文件被选中时显示
- ✅ 显示选中文件数量
- ✅ 点击按钮显示提示(目前是占位功能)
- ✅ 提示说明将整合到 Runs 页面

## 完整测试流程

### 流程 A: 转录组数据处理

1. 登录系统
2. 选择一个项目
3. 进入 Data Catalog 页面
4. 上传 `test_data/sample_transcriptomics.csv`
   - ✅ 勾选 "Process file immediately"
5. 等待上传完成
6. 验证文件状态显示 "✓ Processed" 和 "transcriptomics" 标签
7. 点击 "👁️ View Data" 查看处理后数据
8. 验证显示:
   - Metadata: omics_type = transcriptomics
   - 5 条基因表达记录
   - 包含 gene_id, gene_name, sample1-3 列
9. 勾选该文件
10. 验证 "🚀 Create Analysis Run" 按钮出现

### 流程 B: 蛋白质组数据处理

1. 上传 `test_data/sample_proteomics.csv`
2. 验证状态显示 "proteomics" 标签
3. 查看处理后数据
4. 验证显示 4 条蛋白质记录
5. 包含 protein_id, mz, intensity, charge 列

### 流程 C: 多文件选择

1. 上传多个不同类型的文件(transcriptomics + proteomics)
2. 都勾选处理
3. 使用全选复选框选中所有已处理文件
4. 验证可以创建多组学分析运行

### 流程 D: 原始文件上传(不处理)

1. 上传文件时取消勾选 "Process file immediately"
2. 验证文件状态显示 "Raw only"
3. 验证该文件的复选框被禁用
4. 验证没有 "👁️ View Data" 按钮

## UI/UX 检查点

### 视觉设计

- [ ] 处理状态用颜色区分(绿色=已处理,灰色=未处理)
- [ ] Omics type 标签有明显的背景色
- [ ] 复选框对齐且易于点击
- [ ] 按钮间距合理,不拥挤
- [ ] 模态窗口居中且有半透明背景遮罩

### 交互体验

- [ ] 复选框点击响应快速
- [ ] "View Data" 按钮点击后数据加载有提示
- [ ] 模态窗口可以滚动查看长内容
- [ ] 关闭按钮位置明显且易于点击
- [ ] 表格数据可以水平滚动(如果列太多)

### 错误处理

- [ ] 加载处理数据失败时显示错误提示
- [ ] 网络错误时有友好的错误消息
- [ ] 无数据时显示空状态提示

## API 调用验证

### 上传文件

```http
POST /data/upload
Content-Type: multipart/form-data

- file: [文件]
- project_id: [项目ID]
- process_file: true/false
```

**验证点**:

- ✅ process_file=true 时返回 processing_info
- ✅ processing_info 包含 omics_type, processed, converted
- ✅ metadata\_ 包含 processed_file_id

### 获取处理数据

```http
GET /data/{raw_file_id}/processed
```

**验证点**:

- ✅ 返回 unified_data 对象
- ✅ unified_data.metadata 包含 omics_type, source 等
- ✅ unified_data.records 是数组
- ✅ processing_info 包含处理详情

## 浏览器测试

### Chrome/Edge

- [ ] 所有功能正常
- [ ] 表格渲染正确
- [ ] 模态窗口正常

### Firefox

- [ ] 所有功能正常
- [ ] 复选框样式一致

### Safari

- [ ] 所有功能正常
- [ ] 文件上传正常

## 性能测试

### 大文件处理

- [ ] 上传 100MB 文件时 UI 不卡顿
- [ ] 处理进度有反馈(目前立即返回)
- [ ] 查看大数据集时只显示前 10 条(性能良好)

### 多文件操作

- [ ] 选择 100 个文件时复选框响应正常
- [ ] 文件列表滚动流畅

## 已知限制与后续改进

### 当前实现

✅ 上传时可选处理
✅ 状态显示
✅ 查看统一格式数据
✅ 多文件选择
⚠️ 创建分析运行(占位功能)

### 待开发功能

- [ ] 与 Runs 页面集成(选中文件 → 创建 Run)
- [ ] 处理进度实时反馈(异步处理+进度条)
- [ ] 数据导出功能(CSV/TSV/Excel)
- [ ] 数据可视化预览(图表)
- [ ] 过滤和搜索处理后的文件
- [ ] 批量操作(批量处理/删除)

## 集成测试命令

```bash
# 启动所有服务
./docker-start.sh

# 访问前端
open http://localhost:5173

# 检查后端健康
curl http://localhost:8001/healthz

# 运行集成测试
pytest tests/integration/test_end_to_end.py -v
```

## 问题记录

### 问题 1: [示例]

**描述**: ...
**重现步骤**: ...
**预期**: ...
**实际**: ...
**状态**: 待修复/已修复

---

**测试人员**: **********\_\_\_**********
**测试日期**: **********\_\_\_**********
**通过**: ☐ 是 ☐ 否
**备注**: ************\_\_\_************
