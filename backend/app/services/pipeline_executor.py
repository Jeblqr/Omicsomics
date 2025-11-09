"""
Pipeline Execution Engine
Handles the actual execution of analysis pipelines
"""

import asyncio
from datetime import datetime, timezone
from typing import Any, Dict, List
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select
import logging

from app.models.run import Run
from app.models.datafile import DataFile

logger = logging.getLogger(__name__)


class PipelineExecutor:
    """Executes analysis pipelines and updates progress"""

    def __init__(self, db: AsyncSession, run: Run):
        self.db = db
        self.run = run
        self.steps_completed = 0
        self.total_steps = 0

    async def execute(self):
        """Execute the pipeline and update progress"""
        try:
            # Update status to running
            self.run.status = "running"
            self.run.started_at = datetime.now(timezone.utc)
            self.run.logs = "Pipeline execution started\n"
            self.run.progress = 0.0
            await self.db.commit()

            # Get pipeline configuration
            pipeline_config = self.run.pipeline_config or {}
            pipeline_type = self.run.pipeline_type

            # Log start
            await self._log(f"Starting {pipeline_type} pipeline: {self.run.name}")
            await self._log(f"Pipeline template: {self.run.pipeline_template_id}")
            await self._log(f"Input files: {len(self.run.input_files or [])}")
            await self._log(f"Parameters: {self.run.parameters}")

            # Execute based on pipeline type
            if pipeline_type == "template":
                await self._execute_template_pipeline()
            elif pipeline_type == "custom":
                await self._execute_custom_pipeline()
            elif pipeline_type == "merged":
                await self._execute_merged_pipeline()
            else:
                await self._log(
                    f"Unknown pipeline type: {pipeline_type}", level="ERROR"
                )
                raise ValueError(f"Unknown pipeline type: {pipeline_type}")

            # Mark as completed
            self.run.status = "completed"
            self.run.progress = 100.0
            self.run.finished_at = datetime.now(timezone.utc)
            await self._log("Pipeline execution completed successfully")
            await self.db.commit()

        except Exception as e:
            logger.exception(f"Pipeline execution failed for run {self.run.id}")
            self.run.status = "failed"
            self.run.error_message = str(e)
            self.run.finished_at = datetime.now(timezone.utc)
            await self._log(f"Pipeline execution failed: {str(e)}", level="ERROR")
            await self.db.commit()
            raise

    async def _execute_template_pipeline(self):
        """Execute a template-based pipeline"""
        template_id = self.run.pipeline_template_id

        # Define steps based on template
        steps = self._get_template_steps(template_id)
        self.total_steps = len(steps)

        await self._log(f"Executing template pipeline with {self.total_steps} steps")

        for i, step in enumerate(steps):
            await self._execute_step(step, i + 1)

    async def _execute_custom_pipeline(self):
        """Execute a custom pipeline"""
        pipeline_config = self.run.pipeline_config or {}
        definition = pipeline_config.get("definition", {})

        if not definition:
            # Try to get from custom_pipeline_id
            if self.run.custom_pipeline_id:
                await self._log(
                    f"Loading custom pipeline {self.run.custom_pipeline_id}"
                )
                # TODO: Load custom pipeline definition from database
                nodes = []
                edges = []
            else:
                nodes = []
                edges = []
        else:
            nodes = definition.get("nodes", [])
            edges = definition.get("edges", [])

        self.total_steps = len(nodes)
        await self._log(f"Executing custom pipeline with {self.total_steps} nodes")

        # Execute nodes in topological order based on edges
        executed_nodes = set()
        node_results = {}

        for i, node in enumerate(nodes):
            await self._execute_custom_node(node, i + 1, node_results)
            executed_nodes.add(node.get("id"))

    async def _execute_merged_pipeline(self):
        """Execute a merged workflow pipeline"""
        await self._log("Executing merged pipeline")

        # Execute template pipeline first if specified
        if self.run.pipeline_template_id:
            await self._log(f"Executing template part: {self.run.pipeline_template_id}")
            template_steps = self._get_template_steps(self.run.pipeline_template_id)
            self.total_steps += len(template_steps)

            for i, step in enumerate(template_steps):
                await self._execute_step(step, i + 1)

        # Then execute custom pipeline if specified
        if self.run.custom_pipeline_id:
            await self._log(f"Executing custom part: {self.run.custom_pipeline_id}")
            # TODO: Load and execute custom pipeline
            pass

    def _get_template_steps(self, template_id: str) -> List[Dict[str, Any]]:
        """Get steps for a template pipeline"""
        # Define steps for each template type
        template_steps = {
            "rna-seq-basic": [
                {"name": "Quality Control", "tool": "fastqc", "duration": 2},
                {"name": "Adapter Trimming", "tool": "trimmomatic", "duration": 3},
                {"name": "Alignment", "tool": "star", "duration": 5},
                {"name": "Quantification", "tool": "featurecounts", "duration": 2},
                {"name": "Differential Expression", "tool": "deseq2", "duration": 3},
            ],
            "variant-calling": [
                {"name": "Quality Control", "tool": "fastqc", "duration": 2},
                {"name": "Read Mapping", "tool": "bwa", "duration": 4},
                {"name": "Sort and Index", "tool": "samtools", "duration": 2},
                {"name": "Variant Calling", "tool": "gatk", "duration": 5},
                {"name": "Variant Filtering", "tool": "gatk", "duration": 2},
                {"name": "Annotation", "tool": "snpeff", "duration": 3},
            ],
            "chip-seq": [
                {"name": "Quality Control", "tool": "fastqc", "duration": 2},
                {"name": "Alignment", "tool": "bowtie2", "duration": 4},
                {"name": "Remove Duplicates", "tool": "picard", "duration": 2},
                {"name": "Peak Calling", "tool": "macs2", "duration": 3},
                {"name": "Motif Analysis", "tool": "homer", "duration": 4},
            ],
            "proteomics-label-free": [
                {"name": "Raw Data Processing", "tool": "maxquant", "duration": 3},
                {"name": "Protein Identification", "tool": "mascot", "duration": 4},
                {"name": "Quantification", "tool": "maxquant", "duration": 3},
                {"name": "Statistical Analysis", "tool": "perseus", "duration": 2},
                {"name": "Pathway Analysis", "tool": "reactome", "duration": 2},
            ],
            "metabolomics-untargeted": [
                {"name": "Peak Detection", "tool": "xcms", "duration": 3},
                {"name": "Peak Alignment", "tool": "xcms", "duration": 2},
                {"name": "Normalization", "tool": "metaboanalyst", "duration": 2},
                {
                    "name": "Statistical Analysis",
                    "tool": "metaboanalyst",
                    "duration": 2,
                },
                {"name": "Metabolite Identification", "tool": "hmdb", "duration": 3},
            ],
            "single-cell-rna": [
                {"name": "Quality Control", "tool": "cellranger", "duration": 3},
                {"name": "Normalization", "tool": "seurat", "duration": 2},
                {"name": "Dimensionality Reduction", "tool": "seurat", "duration": 3},
                {"name": "Clustering", "tool": "seurat", "duration": 2},
                {"name": "Marker Identification", "tool": "seurat", "duration": 3},
            ],
            "gwas": [
                {"name": "Quality Control", "tool": "plink", "duration": 2},
                {"name": "Imputation", "tool": "impute2", "duration": 4},
                {"name": "Association Testing", "tool": "plink", "duration": 5},
                {"name": "Multiple Testing Correction", "tool": "plink", "duration": 2},
                {"name": "Annotation", "tool": "annovar", "duration": 3},
            ],
            "metagenomics": [
                {"name": "Quality Control", "tool": "fastqc", "duration": 2},
                {"name": "Host Removal", "tool": "bowtie2", "duration": 3},
                {"name": "Taxonomic Classification", "tool": "kraken2", "duration": 4},
                {"name": "Abundance Estimation", "tool": "metaphlan", "duration": 3},
                {"name": "Functional Profiling", "tool": "humann", "duration": 4},
            ],
        }

        return template_steps.get(
            template_id,
            [
                {"name": "Step 1", "tool": "generic", "duration": 2},
                {"name": "Step 2", "tool": "generic", "duration": 3},
                {"name": "Step 3", "tool": "generic", "duration": 2},
            ],
        )

    async def _execute_step(self, step: Dict[str, Any], step_number: int):
        """Execute a single pipeline step"""
        step_name = step.get("name", f"Step {step_number}")
        tool = step.get("tool", "unknown")
        duration = step.get("duration", 2)  # seconds

        await self._log(
            f"[{step_number}/{self.total_steps}] Starting: {step_name} (tool: {tool})"
        )

        # Simulate tool execution
        await asyncio.sleep(duration)

        self.steps_completed += 1
        progress = (self.steps_completed / self.total_steps) * 100
        self.run.progress = round(progress, 2)

        await self._log(
            f"[{step_number}/{self.total_steps}] Completed: {step_name} (Progress: {self.run.progress}%)"
        )
        await self.db.commit()

    async def _execute_custom_node(
        self, node: Dict[str, Any], node_number: int, results: Dict[str, Any]
    ):
        """Execute a custom pipeline node"""
        node_id = node.get("id", f"node_{node_number}")
        node_type = node.get("type", "process")
        node_label = node.get("label", f"Node {node_number}")
        node_data = node.get("data", {})

        await self._log(
            f"[{node_number}/{self.total_steps}] Starting node: {node_label} (type: {node_type})"
        )

        # Execute based on node type
        if node_type == "input":
            results[node_id] = await self._handle_input_node(node_data)
        elif node_type == "process":
            results[node_id] = await self._handle_process_node(node_data)
        elif node_type == "output":
            results[node_id] = await self._handle_output_node(node_data, results)
        else:
            results[node_id] = {
                "status": "skipped",
                "reason": f"Unknown node type: {node_type}",
            }

        # Simulate execution time
        await asyncio.sleep(2)

        self.steps_completed += 1
        progress = (self.steps_completed / self.total_steps) * 100
        self.run.progress = round(progress, 2)

        await self._log(
            f"[{node_number}/{self.total_steps}] Completed node: {node_label} (Progress: {self.run.progress}%)"
        )
        await self.db.commit()

    async def _handle_input_node(self, node_data: Dict[str, Any]) -> Dict[str, Any]:
        """Handle input node - load data files"""
        await self._log(f"Loading input data: {node_data}")
        return {"status": "success", "type": "input"}

    async def _handle_process_node(self, node_data: Dict[str, Any]) -> Dict[str, Any]:
        """Handle process node - run tool"""
        tool = node_data.get("tool", "unknown")
        params = node_data.get("parameters", {})
        await self._log(f"Running tool: {tool} with parameters: {params}")
        return {"status": "success", "type": "process", "tool": tool}

    async def _handle_output_node(
        self, node_data: Dict[str, Any], results: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Handle output node - save results"""
        await self._log(f"Saving output: {node_data}")
        return {"status": "success", "type": "output"}

    async def _log(self, message: str, level: str = "INFO"):
        """Add log message to run logs"""
        timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
        log_entry = f"[{timestamp}] [{level}] {message}\n"

        if self.run.logs is None:
            self.run.logs = ""
        self.run.logs += log_entry

        # Also log to application logger
        if level == "ERROR":
            logger.error(f"Run {self.run.id}: {message}")
        elif level == "WARNING":
            logger.warning(f"Run {self.run.id}: {message}")
        else:
            logger.info(f"Run {self.run.id}: {message}")


async def execute_run_async(run_id: int):
    """Execute a run asynchronously in a new database session"""
    from app.database import AsyncSessionMaker
    from app.services.runs import get_run

    # Create a new database session for the background task
    async with AsyncSessionMaker() as db:
        run = await get_run(db, run_id)
        if not run:
            logger.error(f"Run {run_id} not found")
            return

        executor = PipelineExecutor(db, run)
        await executor.execute()
