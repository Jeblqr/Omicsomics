"""Single-cell analysis pipeline for scRNA-seq and scATAC-seq."""

import asyncio
import json
import logging
from pathlib import Path
from typing import Any

from app.models.workflow import WorkflowStatus
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)


class SingleCellAnalyzer:
    """Single-cell analysis pipeline orchestrator."""

    def __init__(self, work_dir: Path | None = None):
        """Initialize single-cell analyzer."""
        self.work_dir = work_dir or Path("/tmp/omicsomics_singlecell")
        self.work_dir.mkdir(parents=True, exist_ok=True)

    async def run_cellranger_count(
        self,
        workflow_id: int,
        fastq_dir: str,
        sample_id: str,
        transcriptome: str,
        output_dir: str,
        chemistry: str = "auto",
        threads: int = 16,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run Cell Ranger count for 10x Genomics data.

        Args:
            workflow_id: Workflow database ID
            fastq_dir: Directory containing FASTQ files
            sample_id: Sample identifier
            transcriptome: Path to Cell Ranger transcriptome reference
            output_dir: Output directory
            chemistry: Chemistry version (auto, SC3Pv3, SC5P-PE, etc.)
            threads: Number of threads
            db: Database session

        Returns:
            Result with Cell Ranger outputs
        """
        cmd = [
            "cellranger",
            "count",
            "--id",
            sample_id,
            "--fastqs",
            fastq_dir,
            "--sample",
            sample_id,
            "--transcriptome",
            transcriptome,
            "--chemistry",
            chemistry,
            "--localcores",
            str(threads),
        ]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running Cell Ranger: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=output_dir,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_path = Path(output_dir) / sample_id / "outs"
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={
                                "filtered_h5": str(
                                    output_path / "filtered_feature_bc_matrix.h5"
                                ),
                                "raw_h5": str(output_path / "raw_feature_bc_matrix.h5"),
                                "web_summary": str(output_path / "web_summary.html"),
                                "metrics": str(output_path / "metrics_summary.csv"),
                            },
                        ),
                    )
                return {
                    "status": "success",
                    "output_path": str(output_path),
                }
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"Cell Ranger failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"Cell Ranger failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_scanpy_preprocessing(
        self,
        workflow_id: int,
        input_h5ad: str,
        output_h5ad: str,
        params: dict[str, Any],
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run Scanpy preprocessing pipeline.

        Args:
            workflow_id: Workflow database ID
            input_h5ad: Input AnnData file (h5ad)
            output_h5ad: Output processed AnnData file
            params: Processing parameters (min_genes, min_cells, etc.)
            db: Database session

        Returns:
            Result with processed AnnData file
        """
        # Python script for Scanpy preprocessing
        python_script = f"""
import scanpy as sc
import numpy as np

# Load data
adata = sc.read_h5ad("{input_h5ad}")

# QC metrics
sc.pp.calculate_qc_metrics(adata, inplace=True)

# Filter cells and genes
min_genes = {params.get('min_genes', 200)}
min_cells = {params.get('min_cells', 3)}
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

# Filter by mitochondrial percentage
if {params.get('filter_mito', True)}:
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    max_pct_mt = {params.get('max_pct_mt', 20)}
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt, :]

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes={params.get('n_top_genes', 2000)})

# Save raw counts
adata.raw = adata

# Scale data
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver='arpack')

# Neighbors and UMAP
sc.pp.neighbors(adata, n_neighbors={params.get('n_neighbors', 10)})
sc.tl.umap(adata)

# Leiden clustering
sc.tl.leiden(adata, resolution={params.get('leiden_resolution', 0.5)})

# Save
adata.write_h5ad("{output_h5ad}")

print("Preprocessing complete!")
print(f"Cells: {{adata.n_obs}}, Genes: {{adata.n_vars}}")
print(f"Clusters: {{len(adata.obs['leiden'].unique())}}")
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_scanpy.py"
        script_path.write_text(python_script)

        cmd = ["python", str(script_path)]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running Scanpy preprocessing: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={"processed_h5ad": output_h5ad},
                        ),
                    )
                return {"status": "success", "processed_h5ad": output_h5ad}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"Scanpy failed: exit code {process.returncode}",
                        ),
                    )
                return {
                    "status": "failed",
                    "error": f"Exit code {process.returncode}",
                    "logs": logs,
                }

        except Exception as e:
            logger.error(f"Scanpy preprocessing failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_seurat_integration(
        self,
        workflow_id: int,
        input_h5_files: list[str],
        output_rds: str,
        batch_key: str = "batch",
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run Seurat integration for batch correction.

        Args:
            workflow_id: Workflow database ID
            input_h5_files: List of input H5 files (one per batch)
            output_rds: Output Seurat object RDS file
            batch_key: Batch identifier column name
            db: Database session

        Returns:
            Result with integrated Seurat object
        """
        # R script for Seurat integration
        r_script = f"""
library(Seurat)
library(SeuratDisk)

# Load data
h5_files <- c({', '.join([f'"{f}"' for f in input_h5_files])})
seurat_list <- list()

for (i in seq_along(h5_files)) {{
  h5_data <- Read10X_h5(h5_files[i])
  seurat_obj <- CreateSeuratObject(counts = h5_data, project = paste0("batch_", i))
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_list[[i]] <- seurat_obj
}}

# Integration
features <- SelectIntegrationFeatures(object.list = seurat_list)
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)

# Standard workflow on integrated data
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)

# Save
saveRDS(integrated, "{output_rds}")

cat("Integration complete!\\n")
cat(paste0("Total cells: ", ncol(integrated), "\\n"))
cat(paste0("Clusters: ", length(unique(integrated$seurat_clusters)), "\\n"))
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_seurat.R"
        script_path.write_text(r_script)

        cmd = ["Rscript", str(script_path)]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running Seurat integration: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={"integrated_rds": output_rds},
                        ),
                    )
                return {"status": "success", "integrated_rds": output_rds}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"Seurat integration failed: exit code {process.returncode}",
                        ),
                    )
                return {
                    "status": "failed",
                    "error": f"Exit code {process.returncode}",
                    "logs": logs,
                }

        except Exception as e:
            logger.error(f"Seurat integration failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_cell_annotation(
        self,
        workflow_id: int,
        input_h5ad: str,
        output_h5ad: str,
        marker_genes: dict[str, list[str]],
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Annotate cell types based on marker genes.

        Args:
            workflow_id: Workflow database ID
            input_h5ad: Input AnnData file
            output_h5ad: Output annotated AnnData file
            marker_genes: Dictionary mapping cell types to marker gene lists
            db: Database session

        Returns:
            Result with annotated AnnData file
        """
        # Python script for cell annotation
        marker_json = json.dumps(marker_genes)
        python_script = f"""
import scanpy as sc
import numpy as np
import json

# Load data
adata = sc.read_h5ad("{input_h5ad}")

# Load marker genes
marker_genes = json.loads('''{marker_json}''')

# Calculate mean expression of marker genes for each cluster
cluster_annotations = {{}}

for cluster in adata.obs['leiden'].unique():
    cluster_cells = adata[adata.obs['leiden'] == cluster]
    cluster_scores = {{}}
    
    for cell_type, markers in marker_genes.items():
        available_markers = [g for g in markers if g in adata.var_names]
        if available_markers:
            score = cluster_cells[:, available_markers].X.mean()
            cluster_scores[cell_type] = score
    
    # Assign cell type with highest score
    if cluster_scores:
        best_type = max(cluster_scores, key=cluster_scores.get)
        cluster_annotations[cluster] = best_type
    else:
        cluster_annotations[cluster] = "Unknown"

# Add annotations to AnnData
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)

# Save
adata.write_h5ad("{output_h5ad}")

print("Annotation complete!")
for ct, count in adata.obs['cell_type'].value_counts().items():
    print(f"{{ct}}: {{count}} cells")
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_annotation.py"
        script_path.write_text(python_script)

        cmd = ["python", str(script_path)]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running cell annotation: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={"annotated_h5ad": output_h5ad},
                        ),
                    )
                return {"status": "success", "annotated_h5ad": output_h5ad}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"Cell annotation failed: exit code {process.returncode}",
                        ),
                    )
                return {
                    "status": "failed",
                    "error": f"Exit code {process.returncode}",
                    "logs": logs,
                }

        except Exception as e:
            logger.error(f"Cell annotation failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}


# Global instance
singlecell_analyzer = SingleCellAnalyzer()
