"""Visualization generation and data preparation for frontend."""

import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class VisualizationGenerator:
    """Generate visualization data for various omics analyses."""

    def __init__(self, output_dir: Path | None = None):
        """Initialize visualization generator."""
        self.output_dir = output_dir or Path("/tmp/omicsomics_viz")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_volcano_plot_data(
        self,
        results_file: str,
        log2fc_col: str = "log2FoldChange",
        pvalue_col: str = "pvalue",
        padj_col: str = "padj",
        gene_col: str = "gene",
        fc_threshold: float = 1.0,
        pval_threshold: float = 0.05,
    ) -> dict[str, Any]:
        """
        Generate volcano plot data from differential expression results.

        Args:
            results_file: Path to DE results CSV/TSV
            log2fc_col: Column name for log2 fold change
            pvalue_col: Column name for p-value
            padj_col: Column name for adjusted p-value
            gene_col: Column name for gene identifiers
            fc_threshold: Log2 fold change threshold for significance
            pval_threshold: Adjusted p-value threshold

        Returns:
            Dictionary with plot data for frontend (Plotly format)
        """
        try:
            import pandas as pd
            import numpy as np

            # Load results
            df = pd.read_csv(results_file)
            df = df.dropna(subset=[log2fc_col, pvalue_col])

            # Calculate -log10(pvalue)
            df["neg_log10_pval"] = -np.log10(df[pvalue_col])

            # Determine significance
            df["significant"] = (np.abs(df[log2fc_col]) > fc_threshold) & (
                df[padj_col] < pval_threshold
            )
            df["direction"] = "not_sig"
            df.loc[
                (df[log2fc_col] > fc_threshold) & (df[padj_col] < pval_threshold),
                "direction",
            ] = "up"
            df.loc[
                (df[log2fc_col] < -fc_threshold) & (df[padj_col] < pval_threshold),
                "direction",
            ] = "down"

            # Prepare data for Plotly
            plot_data = {
                "x": df[log2fc_col].tolist(),
                "y": df["neg_log10_pval"].tolist(),
                "genes": (
                    df[gene_col].tolist()
                    if gene_col in df.columns
                    else df.index.tolist()
                ),
                "colors": df["direction"].tolist(),
                "layout": {
                    "title": "Volcano Plot",
                    "xaxis": {"title": "Log2 Fold Change"},
                    "yaxis": {"title": "-Log10(P-value)"},
                },
                "counts": {
                    "up": int((df["direction"] == "up").sum()),
                    "down": int((df["direction"] == "down").sum()),
                    "not_sig": int((df["direction"] == "not_sig").sum()),
                },
            }

            return {"status": "success", "data": plot_data}

        except Exception as e:
            logger.error(f"Error generating volcano plot data: {e}")
            return {"status": "error", "error": str(e)}

    def generate_pca_plot_data(
        self,
        h5ad_file: str,
        color_by: str = "leiden",
        n_components: int = 2,
    ) -> dict[str, Any]:
        """
        Generate PCA plot data from AnnData object.

        Args:
            h5ad_file: Path to AnnData H5AD file
            color_by: Observation column to color by
            n_components: Number of PCA components to include

        Returns:
            Dictionary with PCA coordinates and metadata
        """
        try:
            import scanpy as sc
            import pandas as pd

            # Load AnnData
            adata = sc.read_h5ad(h5ad_file)

            if "X_pca" not in adata.obsm:
                return {
                    "status": "error",
                    "error": "PCA not computed. Run preprocessing first.",
                }

            # Extract PCA coordinates
            pca_coords = adata.obsm["X_pca"][:, :n_components]

            # Prepare data
            plot_data = {
                "pca": pca_coords.tolist(),
                "metadata": (
                    adata.obs[[color_by]].to_dict(orient="list")
                    if color_by in adata.obs
                    else {}
                ),
                "cell_ids": adata.obs_names.tolist(),
                "layout": {
                    "title": "PCA Plot",
                    "xaxis": {"title": "PC1"},
                    "yaxis": {"title": "PC2"},
                },
            }

            if n_components >= 3:
                plot_data["layout"]["zaxis"] = {"title": "PC3"}

            return {"status": "success", "data": plot_data}

        except Exception as e:
            logger.error(f"Error generating PCA plot data: {e}")
            return {"status": "error", "error": str(e)}

    def generate_umap_plot_data(
        self,
        h5ad_file: str,
        color_by: str = "leiden",
    ) -> dict[str, Any]:
        """
        Generate UMAP plot data from AnnData object.

        Args:
            h5ad_file: Path to AnnData H5AD file
            color_by: Observation column to color by (leiden, cell_type, etc.)

        Returns:
            Dictionary with UMAP coordinates and metadata
        """
        try:
            import scanpy as sc

            # Load AnnData
            adata = sc.read_h5ad(h5ad_file)

            if "X_umap" not in adata.obsm:
                return {
                    "status": "error",
                    "error": "UMAP not computed. Run preprocessing first.",
                }

            # Extract UMAP coordinates
            umap_coords = adata.obsm["X_umap"]

            # Prepare data
            plot_data = {
                "umap": umap_coords.tolist(),
                "metadata": (
                    adata.obs[[color_by]].to_dict(orient="list")
                    if color_by in adata.obs
                    else {}
                ),
                "cell_ids": adata.obs_names.tolist(),
                "layout": {
                    "title": "UMAP Plot",
                    "xaxis": {"title": "UMAP 1"},
                    "yaxis": {"title": "UMAP 2"},
                },
                "cluster_counts": (
                    adata.obs[color_by].value_counts().to_dict()
                    if color_by in adata.obs
                    else {}
                ),
            }

            return {"status": "success", "data": plot_data}

        except Exception as e:
            logger.error(f"Error generating UMAP plot data: {e}")
            return {"status": "error", "error": str(e)}

    def generate_heatmap_data(
        self,
        h5ad_file: str,
        gene_list: list[str],
        groupby: str = "leiden",
        standard_scale: str = "var",
    ) -> dict[str, Any]:
        """
        Generate heatmap data for gene expression.

        Args:
            h5ad_file: Path to AnnData H5AD file
            gene_list: List of genes to include
            groupby: Observation column to group by
            standard_scale: 'var' or 'obs' for scaling

        Returns:
            Dictionary with heatmap matrix and annotations
        """
        try:
            import scanpy as sc
            import numpy as np

            # Load AnnData
            adata = sc.read_h5ad(h5ad_file)

            # Filter genes
            available_genes = [g for g in gene_list if g in adata.var_names]
            if not available_genes:
                return {"status": "error", "error": "No genes found in dataset"}

            # Extract expression matrix
            adata_subset = adata[:, available_genes]

            # Get mean expression per group
            if groupby in adata.obs:
                groups = adata.obs[groupby].unique()
                expr_matrix = []
                for group in groups:
                    group_cells = adata[adata.obs[groupby] == group]
                    group_expr = group_cells[:, available_genes].X.mean(axis=0)
                    if hasattr(group_expr, "A1"):  # Sparse matrix
                        group_expr = group_expr.A1
                    expr_matrix.append(group_expr.tolist())

                plot_data = {
                    "matrix": expr_matrix,
                    "genes": available_genes,
                    "groups": groups.tolist(),
                    "layout": {
                        "title": "Gene Expression Heatmap",
                        "xaxis": {"title": "Genes"},
                        "yaxis": {"title": groupby.capitalize()},
                    },
                }
            else:
                return {
                    "status": "error",
                    "error": f"Column '{groupby}' not found in data",
                }

            return {"status": "success", "data": plot_data}

        except Exception as e:
            logger.error(f"Error generating heatmap data: {e}")
            return {"status": "error", "error": str(e)}

    def generate_igv_data(
        self,
        bam_file: str,
        region: str,
        reference_genome: str = "hg38",
    ) -> dict[str, Any]:
        """
        Prepare data for IGV.js browser.

        Args:
            bam_file: Path to BAM file
            region: Genomic region (e.g., "chr1:1000-2000")
            reference_genome: Reference genome identifier

        Returns:
            Configuration for IGV.js
        """
        try:
            # IGV.js configuration
            config = {
                "genome": reference_genome,
                "locus": region,
                "tracks": [
                    {
                        "name": "Alignments",
                        "type": "alignment",
                        "format": "bam",
                        "url": bam_file,
                        "indexURL": f"{bam_file}.bai",
                    }
                ],
            }

            return {"status": "success", "config": config}

        except Exception as e:
            logger.error(f"Error generating IGV data: {e}")
            return {"status": "error", "error": str(e)}

    def generate_quality_metrics_plot(
        self,
        h5ad_file: str,
        metrics: list[str] = None,
    ) -> dict[str, Any]:
        """
        Generate quality metrics visualization data.

        Args:
            h5ad_file: Path to AnnData H5AD file
            metrics: List of metrics to plot (n_genes, n_counts, pct_counts_mt)

        Returns:
            Dictionary with quality metrics data
        """
        if metrics is None:
            metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

        try:
            import scanpy as sc

            # Load AnnData
            adata = sc.read_h5ad(h5ad_file)

            plot_data = {"metrics": {}}

            for metric in metrics:
                if metric in adata.obs:
                    plot_data["metrics"][metric] = {
                        "values": adata.obs[metric].tolist(),
                        "mean": float(adata.obs[metric].mean()),
                        "median": float(adata.obs[metric].median()),
                        "std": float(adata.obs[metric].std()),
                    }

            return {"status": "success", "data": plot_data}

        except Exception as e:
            logger.error(f"Error generating quality metrics: {e}")
            return {"status": "error", "error": str(e)}


# Global instance
viz_generator = VisualizationGenerator()
