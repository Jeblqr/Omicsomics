/**
 * Pipeline Validation and Management Utilities
 * Validates pipeline integrity, detects issues, and handles import/export
 */

import { Node, Edge } from 'reactflow';
import { PipelineNodeData } from '../components/pipelines/PipelineNode';
import { getToolById } from '../components/pipelines/ToolDefinitions';

export interface ValidationIssue {
  type: 'error' | 'warning' | 'info';
  category: 'connection' | 'parameter' | 'dataflow' | 'dependency';
  nodeId?: string;
  message: string;
  suggestion?: string;
}

export interface ValidationResult {
  valid: boolean;
  issues: ValidationIssue[];
  errors: ValidationIssue[];
  warnings: ValidationIssue[];
}

export interface PipelineExport {
  version: string;
  metadata: {
    name: string;
    description?: string;
    author?: string;
    createdAt: string;
    exportedAt: string;
  };
  nodes: Node<PipelineNodeData>[];
  edges: Edge[];
}

/**
 * Validate pipeline for errors and warnings
 */
export const validatePipeline = (nodes: Node<PipelineNodeData>[], edges: Edge[]): ValidationResult => {
  const issues: ValidationIssue[] = [];

  // 1. Check for isolated nodes (no connections)
  const connectedNodeIds = new Set<string>();
  edges.forEach((edge) => {
    connectedNodeIds.add(edge.source);
    connectedNodeIds.add(edge.target);
  });

  nodes.forEach((node) => {
    if (!connectedNodeIds.has(node.id) && nodes.length > 1) {
      issues.push({
        type: 'warning',
        category: 'connection',
        nodeId: node.id,
        message: `Node "${node.data.label}" is isolated (not connected to other nodes)`,
        suggestion: 'Connect this node to other nodes or remove it',
      });
    }
  });

  // 2. Check for missing input nodes
  const inputNodes = nodes.filter((n) => n.data.nodeType === 'input');
  if (inputNodes.length === 0) {
    issues.push({
      type: 'error',
      category: 'dataflow',
      message: 'Pipeline has no input nodes',
      suggestion: 'Add at least one input node to provide data to the pipeline',
    });
  }

  // 3. Check for missing output nodes
  const outputNodes = nodes.filter((n) => n.data.nodeType === 'output');
  if (outputNodes.length === 0) {
    issues.push({
      type: 'warning',
      category: 'dataflow',
      message: 'Pipeline has no output nodes',
      suggestion: 'Add an output node to save results',
    });
  }

  // 4. Check for cycles (circular dependencies)
  const cycles = detectCycles(nodes, edges);
  if (cycles.length > 0) {
    cycles.forEach((cycle) => {
      issues.push({
        type: 'error',
        category: 'dependency',
        message: `Circular dependency detected: ${cycle.join(' → ')}`,
        suggestion: 'Remove one of the connections to break the cycle',
      });
    });
  }

  // 5. Check for invalid data flow (e.g., output → input)
  edges.forEach((edge) => {
    const sourceNode = nodes.find((n) => n.id === edge.source);
    const targetNode = nodes.find((n) => n.id === edge.target);

    if (sourceNode && targetNode) {
      // Output nodes shouldn't have outgoing connections
      if (sourceNode.data.nodeType === 'output') {
        issues.push({
          type: 'error',
          category: 'dataflow',
          nodeId: sourceNode.id,
          message: `Output node "${sourceNode.data.label}" cannot have outgoing connections`,
          suggestion: 'Remove connections from output nodes',
        });
      }

      // Input nodes shouldn't have incoming connections
      if (targetNode.data.nodeType === 'input') {
        issues.push({
          type: 'error',
          category: 'dataflow',
          nodeId: targetNode.id,
          message: `Input node "${targetNode.data.label}" cannot have incoming connections`,
          suggestion: 'Remove connections to input nodes',
        });
      }

      // Check logical flow order
      const flowOrder = ['input', 'process', 'filter', 'transform', 'analysis', 'output'];
      const sourceIndex = flowOrder.indexOf(sourceNode.data.nodeType);
      const targetIndex = flowOrder.indexOf(targetNode.data.nodeType);

      if (sourceIndex > targetIndex && targetNode.data.nodeType !== 'process') {
        issues.push({
          type: 'warning',
          category: 'dataflow',
          message: `Unusual data flow: ${sourceNode.data.nodeType} → ${targetNode.data.nodeType}`,
          suggestion: 'Verify this connection makes sense for your analysis',
        });
      }
    }
  });

  // 6. Check for required parameters
  nodes.forEach((node) => {
    const toolDef = getToolById(node.id.replace(/^node_\d+_/, ''));
    if (toolDef) {
      const missingParams: string[] = [];
      toolDef.parameterTemplate.forEach((param) => {
        if (param.required && !(param.name in node.data.parameters)) {
          missingParams.push(param.label);
        }
      });

      if (missingParams.length > 0) {
        issues.push({
          type: 'error',
          category: 'parameter',
          nodeId: node.id,
          message: `Node "${node.data.label}" is missing required parameters: ${missingParams.join(', ')}`,
          suggestion: 'Configure the node and fill in all required parameters',
        });
      }
    }
  });

  // 7. Check for nodes with no outgoing connections (except output nodes)
  const nodesWithOutgoing = new Set(edges.map((e) => e.source));
  nodes.forEach((node) => {
    if (!nodesWithOutgoing.has(node.id) && node.data.nodeType !== 'output') {
      issues.push({
        type: 'warning',
        category: 'connection',
        nodeId: node.id,
        message: `Node "${node.data.label}" has no outgoing connections`,
        suggestion: 'Connect this node to downstream processing steps',
      });
    }
  });

  // 8. Check for reasonable pipeline length
  if (nodes.length > 50) {
    issues.push({
      type: 'warning',
      category: 'dataflow',
      message: `Pipeline is very large (${nodes.length} nodes)`,
      suggestion: 'Consider breaking into smaller sub-pipelines for better maintainability',
    });
  }

  const errors = issues.filter((i) => i.type === 'error');
  const warnings = issues.filter((i) => i.type === 'warning');

  return {
    valid: errors.length === 0,
    issues,
    errors,
    warnings,
  };
};

/**
 * Detect cycles in the pipeline graph
 */
function detectCycles(nodes: Node<PipelineNodeData>[], edges: Edge[]): string[][] {
  const adjacencyList = new Map<string, string[]>();
  const cycles: string[][] = [];

  // Build adjacency list
  nodes.forEach((node) => adjacencyList.set(node.id, []));
  edges.forEach((edge) => {
    const neighbors = adjacencyList.get(edge.source) || [];
    neighbors.push(edge.target);
    adjacencyList.set(edge.source, neighbors);
  });

  // DFS to detect cycles
  const visited = new Set<string>();
  const recursionStack = new Set<string>();
  const path: string[] = [];

  function dfs(nodeId: string): boolean {
    visited.add(nodeId);
    recursionStack.add(nodeId);
    path.push(nodeId);

    const neighbors = adjacencyList.get(nodeId) || [];
    for (const neighbor of neighbors) {
      if (!visited.has(neighbor)) {
        if (dfs(neighbor)) return true;
      } else if (recursionStack.has(neighbor)) {
        // Found a cycle
        const cycleStart = path.indexOf(neighbor);
        const cycle = path.slice(cycleStart).map((id) => {
          const node = nodes.find((n) => n.id === id);
          return node?.data.label || id;
        });
        cycles.push([...cycle, cycle[0]]); // Add first node again to show cycle
        return true;
      }
    }

    recursionStack.delete(nodeId);
    path.pop();
    return false;
  }

  nodes.forEach((node) => {
    if (!visited.has(node.id)) {
      dfs(node.id);
    }
  });

  return cycles;
}

/**
 * Export pipeline to JSON
 */
export const exportPipeline = (
  nodes: Node<PipelineNodeData>[],
  edges: Edge[],
  metadata: {
    name: string;
    description?: string;
    author?: string;
    createdAt?: string;
  }
): PipelineExport => {
  return {
    version: '1.0',
    metadata: {
      ...metadata,
      createdAt: metadata.createdAt || new Date().toISOString(),
      exportedAt: new Date().toISOString(),
    },
    nodes,
    edges,
  };
};

/**
 * Export pipeline as downloadable JSON file
 */
export const downloadPipelineJSON = (
  nodes: Node<PipelineNodeData>[],
  edges: Edge[],
  metadata: {
    name: string;
    description?: string;
    author?: string;
  }
): void => {
  const pipeline = exportPipeline(nodes, edges, metadata);
  const data = JSON.stringify(pipeline, null, 2);
  const blob = new Blob([data], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = `pipeline-${metadata.name.toLowerCase().replace(/\s+/g, '-')}-${Date.now()}.json`;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};

/**
 * Import pipeline from JSON
 */
export const importPipeline = (file: File): Promise<PipelineExport> => {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();

    reader.onload = (e) => {
      try {
        const imported = JSON.parse(e.target?.result as string) as PipelineExport;

        // Validate structure
        if (!imported.nodes || !imported.edges || !imported.metadata) {
          throw new Error('Invalid pipeline file format');
        }

        // Validate version compatibility
        if (imported.version !== '1.0') {
          console.warn(`Pipeline version ${imported.version} may not be fully compatible`);
        }

        resolve(imported);
      } catch (error) {
        reject(error);
      }
    };

    reader.onerror = () => reject(reader.error);
    reader.readAsText(file);
  });
};

/**
 * Merge multiple pipelines (for reference, not auto-execution)
 */
export const mergePipelines = (
  pipelines: PipelineExport[],
  layout: 'horizontal' | 'vertical' = 'horizontal'
): { nodes: Node<PipelineNodeData>[]; edges: Edge[] } => {
  const allNodes: Node<PipelineNodeData>[] = [];
  const allEdges: Edge[] = [];
  let offsetX = 0;
  let offsetY = 0;

  pipelines.forEach((pipeline, index) => {
    // Offset nodes to avoid overlap
    const offset = layout === 'horizontal' ? { x: offsetX, y: 0 } : { x: 0, y: offsetY };

    // Add prefix to node IDs to avoid collisions
    const prefix = `p${index}_`;
    const offsetNodes = pipeline.nodes.map((node) => ({
      ...node,
      id: prefix + node.id,
      position: {
        x: node.position.x + offset.x,
        y: node.position.y + offset.y,
      },
    }));

    const offsetEdges = pipeline.edges.map((edge) => ({
      ...edge,
      id: prefix + edge.id,
      source: prefix + edge.source,
      target: prefix + edge.target,
    }));

    allNodes.push(...offsetNodes);
    allEdges.push(...offsetEdges);

    // Calculate next offset
    if (layout === 'horizontal') {
      const maxX = Math.max(...offsetNodes.map((n) => n.position.x));
      offsetX = maxX + 400; // Add spacing
    } else {
      const maxY = Math.max(...offsetNodes.map((n) => n.position.y));
      offsetY = maxY + 300; // Add spacing
    }
  });

  return { nodes: allNodes, edges: allEdges };
};

/**
 * Get predefined pipeline templates
 */
export const getPipelineTemplates = (): PipelineExport[] => {
  return [
    {
      version: '1.0',
      metadata: {
        name: 'Basic RNA-seq Analysis',
        description: 'Standard RNA-seq pipeline: QC → Alignment → Quantification → DE Analysis',
        author: 'Omicsomics',
        createdAt: '2025-01-10T00:00:00Z',
        exportedAt: '2025-01-10T00:00:00Z',
      },
      nodes: [
        {
          id: 'template_rnaseq_1',
          type: 'pipelineNode',
          position: { x: 100, y: 100 },
          data: {
            label: 'FASTQ Input',
            tool: 'file-input',
            version: '1.0',
            parameters: { paired_end: true },
            nodeType: 'input',
          },
        },
        {
          id: 'template_rnaseq_2',
          type: 'pipelineNode',
          position: { x: 100, y: 200 },
          data: {
            label: 'FastQC',
            tool: 'fastqc',
            version: '0.12.1',
            parameters: { threads: 4, kmers: 7 },
            nodeType: 'process',
          },
        },
        {
          id: 'template_rnaseq_3',
          type: 'pipelineNode',
          position: { x: 100, y: 300 },
          data: {
            label: 'Trim Galore',
            tool: 'trim_galore',
            version: '0.6.10',
            parameters: { quality: 20, length: 20 },
            nodeType: 'filter',
          },
        },
        {
          id: 'template_rnaseq_4',
          type: 'pipelineNode',
          position: { x: 350, y: 300 },
          data: {
            label: 'STAR Alignment',
            tool: 'star',
            version: '2.7.11',
            parameters: { threads: 8, quant_mode: 'GeneCounts' },
            nodeType: 'process',
          },
        },
        {
          id: 'template_rnaseq_5',
          type: 'pipelineNode',
          position: { x: 350, y: 400 },
          data: {
            label: 'featureCounts',
            tool: 'featureCounts',
            version: '2.0.6',
            parameters: { threads: 4, is_paired_end: true },
            nodeType: 'process',
          },
        },
        {
          id: 'template_rnaseq_6',
          type: 'pipelineNode',
          position: { x: 600, y: 400 },
          data: {
            label: 'DESeq2',
            tool: 'deseq2',
            version: '1.40.0',
            parameters: { alpha: 0.05, lfc_threshold: 0 },
            nodeType: 'analysis',
          },
        },
        {
          id: 'template_rnaseq_7',
          type: 'pipelineNode',
          position: { x: 600, y: 500 },
          data: {
            label: 'CSV Output',
            tool: 'file-output',
            version: '1.0',
            parameters: { format: 'csv' },
            nodeType: 'output',
          },
        },
      ],
      edges: [
        { id: 'e1-2', source: 'template_rnaseq_1', target: 'template_rnaseq_2', animated: true },
        { id: 'e2-3', source: 'template_rnaseq_2', target: 'template_rnaseq_3', animated: true },
        { id: 'e3-4', source: 'template_rnaseq_3', target: 'template_rnaseq_4', animated: true },
        { id: 'e4-5', source: 'template_rnaseq_4', target: 'template_rnaseq_5', animated: true },
        { id: 'e5-6', source: 'template_rnaseq_5', target: 'template_rnaseq_6', animated: true },
        { id: 'e6-7', source: 'template_rnaseq_6', target: 'template_rnaseq_7', animated: true },
      ],
    },
    {
      version: '1.0',
      metadata: {
        name: 'Variant Calling Pipeline',
        description: 'DNA-seq variant calling: QC → Alignment → Variant Calling → VCF Output',
        author: 'Omicsomics',
        createdAt: '2025-01-10T00:00:00Z',
        exportedAt: '2025-01-10T00:00:00Z',
      },
      nodes: [
        {
          id: 'template_variant_1',
          type: 'pipelineNode',
          position: { x: 100, y: 100 },
          data: {
            label: 'FASTQ Input',
            tool: 'file-input',
            version: '1.0',
            parameters: { paired_end: true },
            nodeType: 'input',
          },
        },
        {
          id: 'template_variant_2',
          type: 'pipelineNode',
          position: { x: 100, y: 200 },
          data: {
            label: 'FastQC',
            tool: 'fastqc',
            version: '0.12.1',
            parameters: { threads: 4 },
            nodeType: 'process',
          },
        },
        {
          id: 'template_variant_3',
          type: 'pipelineNode',
          position: { x: 350, y: 200 },
          data: {
            label: 'BWA-MEM',
            tool: 'bwa',
            version: '0.7.17',
            parameters: { threads: 8 },
            nodeType: 'process',
          },
        },
        {
          id: 'template_variant_4',
          type: 'pipelineNode',
          position: { x: 600, y: 200 },
          data: {
            label: 'GATK HaplotypeCaller',
            tool: 'gatk_haplotypecaller',
            version: '4.4.0',
            parameters: { ploidy: 2 },
            nodeType: 'analysis',
          },
        },
        {
          id: 'template_variant_5',
          type: 'pipelineNode',
          position: { x: 600, y: 300 },
          data: {
            label: 'VCF Output',
            tool: 'file-output',
            version: '1.0',
            parameters: { format: 'vcf', compress: true },
            nodeType: 'output',
          },
        },
      ],
      edges: [
        { id: 'e1-2', source: 'template_variant_1', target: 'template_variant_2', animated: true },
        { id: 'e2-3', source: 'template_variant_2', target: 'template_variant_3', animated: true },
        { id: 'e3-4', source: 'template_variant_3', target: 'template_variant_4', animated: true },
        { id: 'e4-5', source: 'template_variant_4', target: 'template_variant_5', animated: true },
      ],
    },
  ];
};
