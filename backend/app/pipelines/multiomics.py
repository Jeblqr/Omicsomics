"""Multi-omics integration pipeline for joint analysis."""

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


class MultiOmicsIntegrator:
    """Multi-omics data integration pipeline orchestrator."""

    def __init__(self, work_dir: Path | None = None):
        """Initialize multi-omics integrator."""
        self.work_dir = work_dir or Path("/tmp/omicsomics_multiomics")
        self.work_dir.mkdir(parents=True, exist_ok=True)
        (self.work_dir / "outputs").mkdir(parents=True, exist_ok=True)

    def _sandbox_path(self, path: Path) -> Path:
        """Construct a sandboxed path under the integrator work directory."""
        relative_parts = [part for part in path.parts if part not in ("", "/")]
        return (self.work_dir / "outputs").joinpath(*relative_parts)

    def _prepare_output_dir(self, output_dir: str) -> Path:
        """Ensure the output directory exists, falling back to sandbox if needed."""
        target = Path(output_dir)
        try:
            target.mkdir(parents=True, exist_ok=True)
            return target
        except OSError:
            sandboxed = self._sandbox_path(target)
            sandboxed.mkdir(parents=True, exist_ok=True)
            logger.warning(
                "Using sandboxed output directory %s for requested path %s",
                sandboxed,
                target,
            )
            return sandboxed

    def _prepare_output_file(self, output_file: str) -> Path:
        """Ensure the output file parent exists, falling back to sandbox if needed."""
        target = Path(output_file)
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            return target
        except OSError:
            sandboxed = self._sandbox_path(target)
            sandboxed.parent.mkdir(parents=True, exist_ok=True)
            logger.warning(
                "Using sandboxed output file %s for requested path %s",
                sandboxed,
                target,
            )
            return sandboxed

    async def run_mofa2_integration(
        self,
        workflow_id: int,
        data_matrices: dict[str, str],
        output_dir: str,
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run MOFA2 (Multi-Omics Factor Analysis v2) for data integration.

        Args:
            workflow_id: Workflow database ID
            data_matrices: Dictionary of omics type to data matrix file path
                          e.g., {"transcriptomics": "rna.csv", "proteomics": "protein.csv"}
            output_dir: Output directory
            params: MOFA2 parameters
            db: Database session

        Returns:
            Result with integrated factors
        """
        output_path = self._prepare_output_dir(output_dir)
        output_path_str = str(output_path)
        params = params or {}

        # R script for MOFA2
        r_script = self._generate_mofa2_script(data_matrices, output_path_str, params)
        script_path = self.work_dir / f"workflow_{workflow_id}_mofa2.R"
        script_path.write_text(r_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["Rscript", str(script_path)]
            logger.info(f"Running MOFA2: {' '.join(cmd)}")

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "model": f"{output_path_str}/mofa2_model.hdf5",
                    "factors": f"{output_path_str}/mofa2_factors.csv",
                    "weights": f"{output_path_str}/mofa2_weights.csv",
                    "variance": f"{output_path_str}/mofa2_variance.csv",
                    "plots": f"{output_path_str}/mofa2_plots.pdf",
                }

                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files=output_files,
                        ),
                    )
                return {"status": "success", **output_files}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"MOFA2 failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"MOFA2 integration failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    def _generate_mofa2_script(
        self, data_matrices: dict[str, str], output_dir: str, params: dict[str, Any]
    ) -> str:
        """Generate R script for MOFA2."""
        n_factors = params.get("n_factors", 10)
        convergence_mode = params.get("convergence_mode", "fast")

        # Build data loading code
        data_load_lines = []
        for view_name, file_path in data_matrices.items():
            data_load_lines.append(
                f'{view_name}_data <- as.matrix(read.csv("{file_path}", row.names=1))'
            )

        view_names = list(data_matrices.keys())
        data_list = ", ".join([f"{v}_data = {v}_data" for v in view_names])

        r_script = f"""
library(MOFA2)
library(ggplot2)

# Load omics data matrices
{chr(10).join(data_load_lines)}

# Create MOFA object
MOFAobject <- create_mofa(list({data_list}))

# Define data options
data_opts <- get_default_data_options(MOFAobject)

# Define model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- {n_factors}

# Define training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "{convergence_mode}"
train_opts$seed <- 42

# Prepare MOFA object
MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
)

# Train model
MOFAobject <- run_mofa(MOFAobject, outfile = "{output_dir}/mofa2_model.hdf5")

# Extract factors
factors <- get_factors(MOFAobject, as.data.frame = TRUE)
write.csv(factors, "{output_dir}/mofa2_factors.csv", row.names = FALSE)

# Extract weights
weights <- get_weights(MOFAobject, as.data.frame = TRUE)
write.csv(weights, "{output_dir}/mofa2_weights.csv", row.names = FALSE)

# Calculate variance explained
variance <- calculate_variance_explained(MOFAobject)
variance_df <- as.data.frame(variance$r2_per_factor)
write.csv(variance_df, "{output_dir}/mofa2_variance.csv", row.names = TRUE)

# Generate plots
pdf("{output_dir}/mofa2_plots.pdf", width = 10, height = 8)

# Plot variance explained
plot_variance_explained(MOFAobject, plot_total = TRUE)

# Plot factor values
plot_factor(MOFAobject, factors = 1:min(5, {n_factors}))

# Plot top weights
for (view in views_names(MOFAobject)) {{
    print(plot_top_weights(
        MOFAobject,
        view = view,
        factor = 1,
        nfeatures = 10
    ))
}}

# Correlation heatmap
plot_factor_cor(MOFAobject)

dev.off()

cat("MOFA2 integration completed successfully\\n")
cat("Factors:", model_opts$num_factors, "\\n")
cat("Views:", length(views_names(MOFAobject)), "\\n")
"""
        return r_script

    async def run_diablo_integration(
        self,
        workflow_id: int,
        data_matrices: dict[str, str],
        phenotype_file: str,
        output_dir: str,
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run DIABLO (Data Integration Analysis for Biomarker discovery using Latent cOmponents).

        Args:
            workflow_id: Workflow database ID
            data_matrices: Dictionary of omics type to data matrix file path
            phenotype_file: CSV file with sample phenotype/class labels
            output_dir: Output directory
            params: DIABLO parameters
            db: Database session

        Returns:
            Result with integrated components and biomarkers
        """
        output_path = self._prepare_output_dir(output_dir)
        output_path_str = str(output_path)
        params = params or {}

        # R script for DIABLO
        r_script = self._generate_diablo_script(
            data_matrices, phenotype_file, output_path_str, params
        )
        script_path = self.work_dir / f"workflow_{workflow_id}_diablo.R"
        script_path.write_text(r_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["Rscript", str(script_path)]
            logger.info(f"Running DIABLO: {' '.join(cmd)}")

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "model": f"{output_path_str}/diablo_model.rds",
                    "components": f"{output_path_str}/diablo_components.csv",
                    "loadings": f"{output_path_str}/diablo_loadings.csv",
                    "performance": f"{output_path_str}/diablo_performance.csv",
                    "selected_features": f"{output_path_str}/diablo_selected_features.csv",
                    "plots": f"{output_path_str}/diablo_plots.pdf",
                }

                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files=output_files,
                        ),
                    )
                return {"status": "success", **output_files}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"DIABLO failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"DIABLO integration failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    def _generate_diablo_script(
        self,
        data_matrices: dict[str, str],
        phenotype_file: str,
        output_dir: str,
        params: dict[str, Any],
    ) -> str:
        """Generate R script for DIABLO."""
        n_comp = params.get("n_components", 2)
        design_correlation = params.get("design_correlation", 0.1)

        # Build data loading code
        data_load_lines = []
        for view_name, file_path in data_matrices.items():
            data_load_lines.append(
                f'{view_name}_data <- as.matrix(read.csv("{file_path}", row.names=1))'
            )

        view_names = list(data_matrices.keys())

        r_script = f"""
library(mixOmics)
library(ggplot2)

# Load omics data matrices
{chr(10).join(data_load_lines)}

# Load phenotype data
phenotype <- read.csv("{phenotype_file}", row.names=1)
Y <- as.factor(phenotype$class)  # Assuming column named 'class'

# Create list of data matrices
X <- list(
{', '.join([f'  {v} = {v}_data' for v in view_names])}
)

# Define design matrix (correlation between datasets)
design <- matrix({design_correlation}, ncol = length(X), nrow = length(X))
diag(design) <- 0
rownames(design) <- colnames(design) <- names(X)

# Initial DIABLO model
diablo_model <- block.plsda(
    X = X,
    Y = Y,
    ncomp = {n_comp},
    design = design
)

# Parameter tuning (feature selection)
test_keepX <- list()
for (block in names(X)) {{
    test_keepX[[block]] <- c(5, 10, 15, 20, 30, 50)
}}

tune_result <- tune.block.splsda(
    X = X,
    Y = Y,
    ncomp = {n_comp},
    test.keepX = test_keepX,
    design = design,
    validation = "Mfold",
    folds = 5,
    nrepeat = 10,
    dist = "centroids.dist",
    cpus = 4
)

# Final model with optimal parameters
optimal_keepX <- tune_result$choice.keepX
final_model <- block.splsda(
    X = X,
    Y = Y,
    ncomp = {n_comp},
    keepX = optimal_keepX,
    design = design
)

# Save model
saveRDS(final_model, "{output_dir}/diablo_model.rds")

# Extract components
components <- final_model$variates
components_df <- do.call(cbind, lapply(names(components), function(block) {{
    df <- as.data.frame(components[[block]])
    colnames(df) <- paste0(block, "_comp", 1:{n_comp})
    df
}}))
components_df$sample <- rownames(components_df)
components_df$class <- Y
write.csv(components_df, "{output_dir}/diablo_components.csv", row.names = FALSE)

# Extract loadings
loadings <- final_model$loadings
loadings_df <- do.call(rbind, lapply(names(loadings), function(block) {{
    df <- as.data.frame(loadings[[block]])
    df$block <- block
    df$feature <- rownames(df)
    df
}}))
write.csv(loadings_df, "{output_dir}/diablo_loadings.csv", row.names = FALSE)

# Performance metrics
perf <- perf(final_model, validation = "Mfold", folds = 5, nrepeat = 10)
perf_df <- data.frame(
    component = 1:{n_comp},
    error_rate = perf$error.rate$overall["centroids.dist", ],
    auc = if (!is.null(perf$auc)) perf$auc$comp1 else NA
)
write.csv(perf_df, "{output_dir}/diablo_performance.csv", row.names = FALSE)

# Selected features
selected_features <- selectVar(final_model, comp = 1)
selected_df <- do.call(rbind, lapply(names(selected_features), function(block) {{
    if (nrow(selected_features[[block]]$value) > 0) {{
        df <- selected_features[[block]]$value
        df$block <- block
        df$feature <- rownames(df)
        df
    }}
}}))
write.csv(selected_df, "{output_dir}/diablo_selected_features.csv", row.names = FALSE)

# Generate plots
pdf("{output_dir}/diablo_plots.pdf", width = 10, height = 8)

# Sample plot
plotIndiv(final_model, ind.names = FALSE, legend = TRUE, title = "DIABLO Sample Plot")

# Variable plot (arrow plot)
plotVar(final_model, cutoff = 0.5, title = "DIABLO Variable Plot")

# Loadings plot
plotLoadings(final_model, comp = 1, contrib = "max")

# Circos plot
circosPlot(final_model, cutoff = 0.7, line = TRUE)

# Heatmap for each block
for (block in names(X)) {{
    cim(final_model, comp = 1:{n_comp}, margins = c(5, 10))
}}

dev.off()

cat("DIABLO integration completed successfully\\n")
cat("Components:", {n_comp}, "\\n")
cat("Blocks:", length(X), "\\n")
cat("Selected features per block:\\n")
print(optimal_keepX)
"""
        return r_script

    async def run_pathway_enrichment(
        self,
        workflow_id: int,
        feature_lists: dict[str, str],
        organism: str,
        output_dir: str,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run pathway enrichment analysis across multiple omics layers.

        Args:
            workflow_id: Workflow database ID
            feature_lists: Dictionary of omics type to feature list file
            organism: Organism (e.g., "hsapiens", "mmusculus")
            output_dir: Output directory
            db: Database session

        Returns:
            Result with enriched pathways
        """
        output_path = self._prepare_output_dir(output_dir)
        output_path_str = str(output_path)

        # Python script for pathway enrichment
        python_script = f"""
import pandas as pd
import json

# Note: This is a placeholder for pathway enrichment
# Full implementation would use gprofiler, DAVID, or enrichR APIs

feature_lists = {json.dumps(feature_lists)}
organism = "{organism}"

enrichment_results = {{}}

for omics_type, feature_file in feature_lists.items():
    print(f"Processing {{omics_type}} features...")
    
    # Load features
    features = pd.read_csv(feature_file)
    
    # Placeholder for enrichment results
    enrichment_results[omics_type] = {{
        "pathways": [],
        "n_features": len(features),
        "message": "Pathway enrichment requires API setup (gprofiler/DAVID/enrichR)"
    }}

# Save results
with open("{output_path_str}/pathway_enrichment.json", "w") as f:
    json.dump(enrichment_results, f, indent=2)

print("Pathway enrichment analysis completed")
print(f"Analyzed {{len(feature_lists)}} omics layers")
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_pathways.py"
        script_path.write_text(python_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["python", str(script_path)]
            logger.info(f"Running pathway enrichment: {' '.join(cmd)}")

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "pathways": f"{output_path_str}/pathway_enrichment.json",
                }

                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files=output_files,
                        ),
                    )
                return {"status": "success", **output_files}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"Pathway enrichment failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"Pathway enrichment failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def match_samples_across_omics(
        self,
        workflow_id: int,
        omics_tables: dict[str, str],
        output_file: str,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Match samples across different omics datasets.

        Args:
            workflow_id: Workflow database ID
            omics_tables: Dictionary of omics type to data table
            output_file: Output matched sample mapping
            db: Database session

        Returns:
            Result with sample mapping
        """
        # Python script for sample matching
        output_path = self._prepare_output_file(output_file)
        output_path_str = str(output_path)

        python_script = f"""
import pandas as pd
import json

omics_tables = {json.dumps(omics_tables)}

# Load all tables and extract sample names
sample_sets = {{}}
for omics_type, table_file in omics_tables.items():
    df = pd.read_csv(table_file, index_col=0)
    sample_sets[omics_type] = set(df.columns)
    print(f"{{omics_type}}: {{len(sample_sets[omics_type])}} samples")

# Find common samples
common_samples = set.intersection(*sample_sets.values())
print(f"\\nCommon samples across all omics: {{len(common_samples)}}")

# Create sample mapping
sample_mapping = {{
    "common_samples": sorted(list(common_samples)),
    "n_common": len(common_samples),
    "omics_specific": {{}}
}}

for omics_type, samples in sample_sets.items():
    specific = samples - common_samples
    sample_mapping["omics_specific"][omics_type] = {{
        "total": len(samples),
        "specific": sorted(list(specific)),
        "n_specific": len(specific)
    }}

# Save mapping
with open("{output_path_str}", "w") as f:
    json.dump(sample_mapping, f, indent=2)

print(f"\nSample mapping saved to {output_path_str}")
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_match_samples.py"
        script_path.write_text(python_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["python", str(script_path)]
            logger.info(f"Running sample matching: {' '.join(cmd)}")

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
                            output_files={"sample_mapping": output_path_str},
                        ),
                    )
                return {"status": "success", "sample_mapping": output_path_str}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"Sample matching failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"Sample matching failed: {e}")
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
multiomics_integrator = MultiOmicsIntegrator()
