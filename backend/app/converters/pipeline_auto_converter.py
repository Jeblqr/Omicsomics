"""
Auto-conversion utilities for pipeline execution.

Detects format mismatches between connected tools and auto-inserts
conversion nodes into the pipeline DAG.
"""

from typing import List, Dict, Optional, Tuple
from app.converters.format_converter import get_format_converter


class PipelineAutoConverter:
    """
    Detect and handle format conversions in pipeline execution.
    """

    def __init__(self):
        self.converter = get_format_converter()

    def analyze_pipeline(self, pipeline: Dict) -> Dict:
        """
        Analyze pipeline for format mismatches.

        Args:
            pipeline: Pipeline definition with nodes and edges

        Returns:
            Analysis result with conversion requirements
        """
        nodes = pipeline.get("nodes", [])
        edges = pipeline.get("edges", [])

        conversions_needed = []
        warnings = []

        # Build node lookup
        node_map = {node["id"]: node for node in nodes}

        # Analyze each connection
        for edge in edges:
            source_node_id = edge["source"]
            target_node_id = edge["target"]
            source_port = edge.get("sourceHandle", "output")
            target_port = edge.get("targetHandle", "input")

            source_node = node_map.get(source_node_id)
            target_node = node_map.get(target_node_id)

            if not source_node or not target_node:
                continue

            # Get output format of source and input format of target
            source_format = self._get_output_format(source_node, source_port)
            target_format = self._get_input_format(target_node, target_port)

            if not source_format or not target_format:
                continue

            # Check if formats match
            if source_format != target_format:
                try:
                    # Get conversion path
                    conversion_path = self.converter.get_conversion_path(
                        source_format, target_format
                    )

                    # Estimate time
                    # Note: Can't estimate without file, use default
                    estimated_time = self._estimate_default_time(
                        source_format, target_format
                    )

                    conversions_needed.append(
                        {
                            "edge_id": edge.get("id"),
                            "source_node": source_node_id,
                            "target_node": target_node_id,
                            "source_format": source_format,
                            "target_format": target_format,
                            "conversion_path": conversion_path,
                            "estimated_seconds": estimated_time,
                            "num_steps": len(conversion_path) - 1,
                        }
                    )

                    if len(conversion_path) > 2:
                        warnings.append(
                            f"Multi-step conversion required between "
                            f"{source_node.get('data', {}).get('label', source_node_id)} "
                            f"and {target_node.get('data', {}).get('label', target_node_id)}: "
                            f"{' → '.join(conversion_path)}"
                        )

                except ValueError as e:
                    warnings.append(
                        f"Cannot convert between "
                        f"{source_node.get('data', {}).get('label', source_node_id)} "
                        f"({source_format}) and "
                        f"{target_node.get('data', {}).get('label', target_node_id)} "
                        f"({target_format}): {str(e)}"
                    )

        return {
            "conversions_needed": conversions_needed,
            "num_conversions": len(conversions_needed),
            "warnings": warnings,
            "total_estimated_time": sum(
                c["estimated_seconds"] for c in conversions_needed
            ),
        }

    def insert_conversion_nodes(self, pipeline: Dict) -> Dict:
        """
        Insert conversion nodes into pipeline where needed.

        Args:
            pipeline: Original pipeline definition

        Returns:
            Modified pipeline with conversion nodes inserted
        """
        analysis = self.analyze_pipeline(pipeline)

        if analysis["num_conversions"] == 0:
            return pipeline

        nodes = list(pipeline.get("nodes", []))
        edges = list(pipeline.get("edges", []))

        # Track new nodes and edges
        new_nodes = []
        new_edges = []
        edges_to_remove = []

        # Generate conversion nodes
        for idx, conversion in enumerate(analysis["conversions_needed"]):
            # Create conversion node
            conversion_node = {
                "id": f"conversion_{idx}",
                "type": "tool",
                "data": {
                    "label": f'Convert {conversion["source_format"]} → {conversion["target_format"]}',
                    "tool_id": "format_converter",
                    "tool_type": "converter",
                    "is_auto_inserted": True,
                    "conversion_info": conversion,
                },
                "position": self._calculate_middle_position(
                    nodes, conversion["source_node"], conversion["target_node"]
                ),
            }
            new_nodes.append(conversion_node)

            # Find the edge to replace
            edge_to_replace = None
            for edge in edges:
                if (
                    edge["source"] == conversion["source_node"]
                    and edge["target"] == conversion["target_node"]
                ):
                    edge_to_replace = edge
                    break

            if edge_to_replace:
                edges_to_remove.append(edge_to_replace)

                # Create new edges: source → converter → target
                new_edges.append(
                    {
                        "id": f"edge_to_conversion_{idx}",
                        "source": conversion["source_node"],
                        "target": conversion_node["id"],
                        "sourceHandle": edge_to_replace.get("sourceHandle", "output"),
                        "targetHandle": "input",
                        "type": "smoothstep",
                        "animated": True,
                        "style": {"stroke": "#fbbf24"},  # Amber color for conversion
                    }
                )

                new_edges.append(
                    {
                        "id": f"edge_from_conversion_{idx}",
                        "source": conversion_node["id"],
                        "target": conversion["target_node"],
                        "sourceHandle": "output",
                        "targetHandle": edge_to_replace.get("targetHandle", "input"),
                        "type": "smoothstep",
                        "animated": True,
                        "style": {"stroke": "#fbbf24"},
                    }
                )

        # Remove old edges
        for edge in edges_to_remove:
            edges.remove(edge)

        # Add new nodes and edges
        nodes.extend(new_nodes)
        edges.extend(new_edges)

        return {
            **pipeline,
            "nodes": nodes,
            "edges": edges,
            "auto_conversion_applied": True,
            "conversion_analysis": analysis,
        }

    def _get_output_format(self, node: Dict, port: str) -> Optional[str]:
        """Get output format for a node's port."""
        tool_data = node.get("data", {})
        outputs = tool_data.get("outputs", [])

        # Find matching output
        for output in outputs:
            if output.get("id") == port or output.get("name") == port:
                return output.get("format")

        # Default output format
        if outputs and len(outputs) > 0:
            return outputs[0].get("format")

        return None

    def _get_input_format(self, node: Dict, port: str) -> Optional[str]:
        """Get input format for a node's port."""
        tool_data = node.get("data", {})
        inputs = tool_data.get("inputs", [])

        # Find matching input
        for input_spec in inputs:
            if input_spec.get("id") == port or input_spec.get("name") == port:
                return input_spec.get("format")

        # Default input format
        if inputs and len(inputs) > 0:
            return inputs[0].get("format")

        return None

    def _estimate_default_time(self, from_format: str, to_format: str) -> float:
        """Estimate conversion time without file (assumes 1GB)."""
        # Use 1GB as default file size
        default_size_gb = 1.0

        conversion_path = self.converter.get_conversion_path(from_format, to_format)
        total_time = 0

        for i in range(len(conversion_path) - 1):
            pair = (conversion_path[i], conversion_path[i + 1])
            time_per_gb = self.converter.CONVERSION_TIME_ESTIMATES.get(pair, 10)

            if pair not in self.converter.CONVERSION_TIME_ESTIMATES:
                reverse_pair = (pair[1], pair[0])
                time_per_gb = self.converter.CONVERSION_TIME_ESTIMATES.get(
                    reverse_pair, 10
                )

            total_time += time_per_gb * default_size_gb

        return max(total_time, 1.0)

    def _calculate_middle_position(
        self, nodes: List[Dict], source_id: str, target_id: str
    ) -> Dict:
        """Calculate position for conversion node between source and target."""
        source_node = next((n for n in nodes if n["id"] == source_id), None)
        target_node = next((n for n in nodes if n["id"] == target_id), None)

        if not source_node or not target_node:
            return {"x": 0, "y": 0}

        source_pos = source_node.get("position", {"x": 0, "y": 0})
        target_pos = target_node.get("position", {"x": 0, "y": 0})

        return {
            "x": (source_pos["x"] + target_pos["x"]) / 2,
            "y": (source_pos["y"] + target_pos["y"]) / 2,
        }


# Singleton instance
_auto_converter = None


def get_pipeline_auto_converter() -> PipelineAutoConverter:
    """Get singleton auto-converter instance."""
    global _auto_converter
    if _auto_converter is None:
        _auto_converter = PipelineAutoConverter()
    return _auto_converter
