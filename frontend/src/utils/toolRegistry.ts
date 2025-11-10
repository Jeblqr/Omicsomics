/**
 * Tool Registry
 * 工具注册中心 - 管理内置工具和用户自定义工具
 * 
 * 功能：
 * 1. 加载内置工具（从ToolDefinitions.ts转换）
 * 2. 用户上传自定义工具JSON
 * 3. 导入/导出工具集合
 * 4. 工具搜索和过滤
 * 5. 持久化到localStorage
 */

import {
  ToolDefinitionSchema,
  ToolCollectionSchema,
  validateToolSchema,
  convertLegacyToolToSchema,
  TOOL_SCHEMA_VERSION,
} from '../schemas/ToolSchema';
import { allTools, getToolById } from '../components/pipelines/ToolDefinitions';

// ===== Storage Keys =====
const STORAGE_KEY_CUSTOM_TOOLS = 'omicsomics_custom_tools';
const STORAGE_KEY_DISABLED_TOOLS = 'omicsomics_disabled_tools';

// ===== Tool Registry =====
class ToolRegistry {
  private customTools: Map<string, ToolDefinitionSchema> = new Map();
  private disabledTools: Set<string> = new Set();
  private builtInTools: Map<string, ToolDefinitionSchema> = new Map();

  constructor() {
    this.loadBuiltInTools();
    this.loadCustomTools();
    this.loadDisabledTools();
  }

  // ===== 加载内置工具 =====
  private loadBuiltInTools() {
    allTools.forEach((tool) => {
      const schema = convertLegacyToolToSchema(tool);
      this.builtInTools.set(schema.id, schema);
    });
    console.log(`Loaded ${this.builtInTools.size} built-in tools`);
  }

  // ===== 从localStorage加载自定义工具 =====
  private loadCustomTools() {
    try {
      const stored = localStorage.getItem(STORAGE_KEY_CUSTOM_TOOLS);
      if (stored) {
        const tools: ToolDefinitionSchema[] = JSON.parse(stored);
        tools.forEach((tool) => {
          const validation = validateToolSchema(tool);
          if (validation.valid) {
            this.customTools.set(tool.id, tool);
          } else {
            console.warn(`Invalid custom tool ${tool.id}:`, validation.errors);
          }
        });
        console.log(`Loaded ${this.customTools.size} custom tools`);
      }
    } catch (error) {
      console.error('Failed to load custom tools:', error);
    }
  }

  // ===== 保存自定义工具到localStorage =====
  private saveCustomTools() {
    try {
      const tools = Array.from(this.customTools.values());
      localStorage.setItem(STORAGE_KEY_CUSTOM_TOOLS, JSON.stringify(tools, null, 2));
    } catch (error) {
      console.error('Failed to save custom tools:', error);
    }
  }

  // ===== 加载禁用工具列表 =====
  private loadDisabledTools() {
    try {
      const stored = localStorage.getItem(STORAGE_KEY_DISABLED_TOOLS);
      if (stored) {
        const disabled: string[] = JSON.parse(stored);
        this.disabledTools = new Set(disabled);
      }
    } catch (error) {
      console.error('Failed to load disabled tools:', error);
    }
  }

  // ===== 保存禁用工具列表 =====
  private saveDisabledTools() {
    try {
      const disabled = Array.from(this.disabledTools);
      localStorage.setItem(STORAGE_KEY_DISABLED_TOOLS, JSON.stringify(disabled));
    } catch (error) {
      console.error('Failed to save disabled tools:', error);
    }
  }

  // ===== 获取所有工具 =====
  getAllTools(includeDisabled = false): ToolDefinitionSchema[] {
    const allTools = [
      ...Array.from(this.builtInTools.values()),
      ...Array.from(this.customTools.values()),
    ];

    if (includeDisabled) {
      return allTools;
    }

    return allTools.filter((tool) => !this.disabledTools.has(tool.id));
  }

  // ===== 按ID获取工具 =====
  getTool(id: string): ToolDefinitionSchema | undefined {
    return this.customTools.get(id) || this.builtInTools.get(id);
  }

  // ===== 按组学类型获取工具 =====
  getToolsByOmicsType(omicsType: string): ToolDefinitionSchema[] {
    return this.getAllTools().filter((tool) => tool.omicsType === omicsType);
  }

  // ===== 按分类获取工具 =====
  getToolsByCategory(category: string): ToolDefinitionSchema[] {
    return this.getAllTools().filter((tool) => tool.category === category);
  }

  // ===== 搜索工具 =====
  searchTools(query: string): ToolDefinitionSchema[] {
    const lowerQuery = query.toLowerCase();
    return this.getAllTools().filter((tool) => {
      return (
        tool.name.toLowerCase().includes(lowerQuery) ||
        tool.description.toLowerCase().includes(lowerQuery) ||
        tool.id.toLowerCase().includes(lowerQuery) ||
        tool.tags?.some((tag) => tag.toLowerCase().includes(lowerQuery))
      );
    });
  }

  // ===== 添加自定义工具 =====
  addCustomTool(tool: ToolDefinitionSchema): { success: boolean; error?: string } {
    // 验证工具定义
    const validation = validateToolSchema(tool);
    if (!validation.valid) {
      return {
        success: false,
        error: `Tool validation failed: ${validation.errors.join(', ')}`,
      };
    }

    // 检查ID冲突
    if (this.builtInTools.has(tool.id)) {
      return {
        success: false,
        error: `Tool ID ${tool.id} conflicts with a built-in tool. Please use a different ID.`,
      };
    }

    if (this.customTools.has(tool.id)) {
      return {
        success: false,
        error: `Tool ID ${tool.id} already exists. Use updateCustomTool() to modify existing tools.`,
      };
    }

    // 添加工具
    this.customTools.set(tool.id, tool);
    this.saveCustomTools();

    return { success: true };
  }

  // ===== 更新自定义工具 =====
  updateCustomTool(id: string, tool: ToolDefinitionSchema): { success: boolean; error?: string } {
    if (!this.customTools.has(id)) {
      return {
        success: false,
        error: `Custom tool ${id} not found`,
      };
    }

    // 验证工具定义
    const validation = validateToolSchema(tool);
    if (!validation.valid) {
      return {
        success: false,
        error: `Tool validation failed: ${validation.errors.join(', ')}`,
      };
    }

    // 如果ID改变了，检查新ID是否冲突
    if (tool.id !== id) {
      if (this.builtInTools.has(tool.id) || this.customTools.has(tool.id)) {
        return {
          success: false,
          error: `New tool ID ${tool.id} conflicts with an existing tool`,
        };
      }
      // 删除旧ID
      this.customTools.delete(id);
    }

    // 更新工具
    this.customTools.set(tool.id, tool);
    this.saveCustomTools();

    return { success: true };
  }

  // ===== 删除自定义工具 =====
  deleteCustomTool(id: string): { success: boolean; error?: string } {
    if (!this.customTools.has(id)) {
      return {
        success: false,
        error: `Custom tool ${id} not found`,
      };
    }

    this.customTools.delete(id);
    this.saveCustomTools();

    return { success: true };
  }

  // ===== 启用/禁用工具 =====
  setToolEnabled(id: string, enabled: boolean) {
    if (enabled) {
      this.disabledTools.delete(id);
    } else {
      this.disabledTools.add(id);
    }
    this.saveDisabledTools();
  }

  isToolEnabled(id: string): boolean {
    return !this.disabledTools.has(id);
  }

  // ===== 导入工具JSON =====
  async importToolFromJSON(file: File): Promise<{ success: boolean; error?: string; tool?: ToolDefinitionSchema }> {
    try {
      const text = await file.text();
      const json = JSON.parse(text);

      // 检查是否为工具集合
      if (json.$schema && json.tools && Array.isArray(json.tools)) {
        return this.importToolCollection(json);
      }

      // 单个工具
      const validation = validateToolSchema(json);
      if (!validation.valid) {
        return {
          success: false,
          error: `Invalid tool JSON: ${validation.errors.join(', ')}`,
        };
      }

      const result = this.addCustomTool(json);
      if (!result.success) {
        return result;
      }

      return {
        success: true,
        tool: json,
      };
    } catch (error) {
      return {
        success: false,
        error: `Failed to parse JSON: ${error}`,
      };
    }
  }

  // ===== 导入工具集合 =====
  importToolCollection(collection: ToolCollectionSchema): { success: boolean; error?: string; imported?: number } {
    let imported = 0;
    const errors: string[] = [];

    collection.tools.forEach((tool) => {
      const result = this.addCustomTool(tool);
      if (result.success) {
        imported++;
      } else {
        errors.push(`${tool.id}: ${result.error}`);
      }
    });

    if (errors.length > 0) {
      return {
        success: false,
        error: `Imported ${imported} tools with ${errors.length} errors: ${errors.join('; ')}`,
        imported,
      };
    }

    return {
      success: true,
      imported,
    };
  }

  // ===== 导出自定义工具为JSON =====
  exportCustomTools(): ToolCollectionSchema {
    return {
      $schema: TOOL_SCHEMA_VERSION,
      name: 'Custom Tools',
      description: 'User-defined custom tools',
      version: '1.0.0',
      tools: Array.from(this.customTools.values()),
    };
  }

  // ===== 导出单个工具 =====
  exportTool(id: string): ToolDefinitionSchema | undefined {
    return this.getTool(id);
  }

  // ===== 下载工具JSON =====
  downloadToolJSON(id: string, filename?: string) {
    const tool = this.exportTool(id);
    if (!tool) {
      console.error(`Tool ${id} not found`);
      return;
    }

    const json = JSON.stringify(tool, null, 2);
    const blob = new Blob([json], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename || `${tool.id}.json`;
    a.click();
    URL.revokeObjectURL(url);
  }

  // ===== 下载所有自定义工具 =====
  downloadAllCustomTools(filename?: string) {
    const collection = this.exportCustomTools();
    const json = JSON.stringify(collection, null, 2);
    const blob = new Blob([json], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename || 'custom-tools.json';
    a.click();
    URL.revokeObjectURL(url);
  }

  // ===== 获取统计信息 =====
  getStats() {
    return {
      builtIn: this.builtInTools.size,
      custom: this.customTools.size,
      disabled: this.disabledTools.size,
      total: this.builtInTools.size + this.customTools.size,
      enabled: this.builtInTools.size + this.customTools.size - this.disabledTools.size,
    };
  }

  // ===== 清除所有自定义工具 =====
  clearAllCustomTools() {
    this.customTools.clear();
    this.saveCustomTools();
  }

  // ===== 重置禁用列表 =====
  resetDisabledTools() {
    this.disabledTools.clear();
    this.saveDisabledTools();
  }
}

// ===== 导出单例实例 =====
export const toolRegistry = new ToolRegistry();

// ===== 导出类型和工具函数 =====
export { ToolRegistry };
export type { ToolDefinitionSchema };
