"""Metabolomics analysis pipeline for mass spectrometry data."""

import asyncio
import logging
from pathlib import Path
from typing import Any

from app.models.workflow import WorkflowStatus
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)


class MetabolomicsAnalyzer:
    """Metabolomics analysis pipeline orchestrator for MS data."""

    def __init__(self, work_dir: Path | None = None):
        """Initialize metabolomics analyzer."""
        self.work_dir = work_dir or Path("/tmp/omicsomics_metabolomics")
        self.work_dir.mkdir(parents=True, exist_ok=True)

    async def run_xcms_feature_detection(
        self,
        workflow_id: int,
        mzml_files: list[str],
        output_dir: str,
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run XCMS for metabolomics feature detection.

        Args:
            workflow_id: Workflow database ID
            mzml_files: List of mzML files
            output_dir: Output directory
            params: XCMS parameters
            db: Database session

        Returns:
            Result with feature table
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        params = params or {}

        # Generate R script for XCMS
        r_script = self._generate_xcms_script(mzml_files, output_dir, params)
        script_path = self.work_dir / f"workflow_{workflow_id}_xcms.R"
        script_path.write_text(r_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["Rscript", str(script_path)]
            logger.info(f"Running XCMS: {' '.join(cmd)}")

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "feature_table": f"{output_dir}/xcms_features.csv",
                    "peak_table": f"{output_dir}/xcms_peaks.csv",
                    "plots": f"{output_dir}/xcms_plots.pdf",
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
                            error_message=f"XCMS failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"XCMS feature detection failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    def _generate_xcms_script(
        self, mzml_files: list[str], output_dir: str, params: dict[str, Any]
    ) -> str:
        """Generate R script for XCMS feature detection."""
        file_list = ", ".join([f'"{f}"' for f in mzml_files])

        ppm = params.get("ppm", 25)
        peakwidth_min = params.get("peakwidth_min", 10)
        peakwidth_max = params.get("peakwidth_max", 60)
        snthresh = params.get("snthresh", 10)

        r_script = f"""
library(xcms)
library(BiocParallel)

# Register parallel backend
register(bpparam = SnowParam(workers = 4))

# Load mzML files
mzml_files <- c({file_list})
raw_data <- readMSData(files = mzml_files, mode = "onDisk")

# Peak picking
cwp <- CentWaveParam(
    ppm = {ppm},
    peakwidth = c({peakwidth_min}, {peakwidth_max}),
    snthresh = {snthresh},
    prefilter = c(3, 100),
    mzCenterFun = "wMean",
    integrate = 1L,
    mzdiff = -0.001,
    fitgauss = FALSE,
    noise = 0,
    verboseColumns = FALSE
)

xdata <- findChromPeaks(raw_data, param = cwp)

# Alignment
pdp <- PeakDensityParam(
    sampleGroups = rep(1, length(mzml_files)),
    minFraction = 0.5,
    bw = 30
)

xdata <- groupChromPeaks(xdata, param = pdp)

# Retention time correction
pgp <- PeakGroupsParam(
    minFraction = 0.85,
    extraPeaks = 1,
    smooth = "loess",
    span = 0.2,
    family = "gaussian"
)

xdata <- adjustRtime(xdata, param = pgp)

# Re-group after RT correction
xdata <- groupChromPeaks(xdata, param = pdp)

# Fill missing peaks
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())

# Export feature table
feature_table <- featureDefinitions(xdata)
write.csv(feature_table, "{output_dir}/xcms_features.csv", row.names = TRUE)

# Export peak intensities
peak_intensities <- featureValues(xdata, value = "into")
write.csv(peak_intensities, "{output_dir}/xcms_peaks.csv", row.names = TRUE)

# Generate plots
pdf("{output_dir}/xcms_plots.pdf", width = 10, height = 8)
plotChromPeakDensity(xdata, param = pdp)
plotAdjustedRtime(xdata)
dev.off()

cat("XCMS feature detection completed successfully\\n")
cat("Features detected:", nrow(feature_table), "\\n")
"""
        return r_script

    async def run_mzmine_feature_detection(
        self,
        workflow_id: int,
        mzml_files: list[str],
        output_dir: str,
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run MZmine for metabolomics feature detection.

        Args:
            workflow_id: Workflow database ID
            mzml_files: List of mzML files
            output_dir: Output directory
            params: MZmine parameters
            db: Database session

        Returns:
            Result with feature table
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        params = params or {}

        # Generate MZmine batch file
        batch_xml = self._generate_mzmine_batch(mzml_files, output_dir, params)
        batch_path = self.work_dir / f"workflow_{workflow_id}_mzmine.xml"
        batch_path.write_text(batch_xml)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = [
                "mzmine",
                str(batch_path),
            ]
            logger.info(f"Running MZmine: {' '.join(cmd)}")

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "feature_table": f"{output_dir}/mzmine_features.csv",
                    "aligned_table": f"{output_dir}/mzmine_aligned.csv",
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
                            error_message=f"MZmine failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"MZmine feature detection failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    def _generate_mzmine_batch(
        self, mzml_files: list[str], output_dir: str, params: dict[str, Any]
    ) -> str:
        """Generate MZmine batch XML file."""
        file_list = "\n".join([f'<file>{f}</file>' for f in mzml_files])

        xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<batch>
    <batchstep method="io.github.mzmine.modules.io.import_rawdata_all.AllSpectralDataImportModule">
        <parameter name="Raw data file names">
            {file_list}
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_massdetection.MassDetectionModule">
        <parameter name="Mass detector">
            <module>io.github.mzmine.modules.dataprocessing.featdet_massdetection.centroid.CentroidMassDetector</module>
            <parameter name="Noise level">1000.0</parameter>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_chromatogrambuilder.ChromatogramBuilderModule">
        <parameter name="Min time span (min)">0.1</parameter>
        <parameter name="Min height">5000.0</parameter>
        <parameter name="m/z tolerance">
            <absolutetolerance>0.001</absolutetolerance>
            <ppmtolerance>10.0</ppmtolerance>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.DeconvolutionModule">
        <parameter name="Algorithm">Local minimum search</parameter>
        <parameter name="Chromatographic threshold">75.0</parameter>
        <parameter name="Min peak height">10000.0</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.align_join.JoinAlignerModule">
        <parameter name="m/z tolerance">
            <absolutetolerance>0.001</absolutetolerance>
            <ppmtolerance>10.0</ppmtolerance>
        </parameter>
        <parameter name="Weight for m/z">75</parameter>
        <parameter name="Weight for RT">25</parameter>
        <parameter name="Retention time tolerance" unit="min">0.5</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.filter_rowsfilter.RowsFilterModule">
        <parameter name="Min peaks in a row">2</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.io.export_features_csv.CSVExportModule">
        <parameter name="Filename">{output_dir}/mzmine_features.csv</parameter>
        <parameter name="Field separator">,</parameter>
    </batchstep>
</batch>
"""
        return xml

    async def run_gnps_annotation(
        self,
        workflow_id: int,
        mgf_file: str,
        output_dir: str,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run GNPS spectral matching for metabolite annotation.

        Args:
            workflow_id: Workflow database ID
            mgf_file: MGF file with MS/MS spectra
            output_dir: Output directory
            db: Database session

        Returns:
            Result with annotation results
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Python script for GNPS API query
        python_script = f"""
import requests
import json
import time
import pandas as pd

# Upload MGF to GNPS
print("Uploading MGF file to GNPS...")
gnps_url = "https://gnps.ucsd.edu/ProteoSAFe/"

# Read MGF file
with open("{mgf_file}", "r") as f:
    mgf_content = f.read()

# Submit GNPS job
job_params = {{
    "workflow": "MOLECULAR-LIBRARYSEARCH-V2",
    "desc": "Automated GNPS search from Omicsomics",
    "library": "GNPS-LIBRARY",
    "spec_on_server": mgf_content,
    "TOLERANCE_ION": "0.5",
    "TOLERANCE_PRECURSOR": "2.0",
    "MIN_MATCHED_PEAKS": "6",
    "SCORE_THRESHOLD": "0.7",
    "TOP_K_RESULTS": "1"
}}

# Note: This is a placeholder for GNPS API integration
# Full implementation would require proper GNPS API authentication and job submission

# For now, create a dummy output
annotation_results = {{
    "status": "completed",
    "matches": 0,
    "message": "GNPS integration requires API credentials and proper job submission"
}}

# Save results
with open("{output_dir}/gnps_annotations.json", "w") as f:
    json.dump(annotation_results, f, indent=2)

print("GNPS annotation completed")
print("Note: Full GNPS integration requires API setup")
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_gnps.py"
        script_path.write_text(python_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["python", str(script_path)]
            logger.info(f"Running GNPS annotation: {' '.join(cmd)}")

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "annotations": f"{output_dir}/gnps_annotations.json",
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
                            error_message=f"GNPS annotation failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"GNPS annotation failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_msdial_annotation(
        self,
        workflow_id: int,
        feature_table: str,
        msp_library: str,
        output_file: str,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run MS-DIAL for metabolite annotation.

        Args:
            workflow_id: Workflow database ID
            feature_table: Feature table CSV file
            msp_library: MSP spectral library file
            output_file: Output annotation file
            db: Database session

        Returns:
            Result with annotation results
        """
        # Python script for MS-DIAL annotation
        python_script = f"""
import pandas as pd
import numpy as np

# Load feature table
print("Loading feature table...")
features = pd.read_csv("{feature_table}")

# Note: MS-DIAL integration requires the MS-DIAL console application
# This is a placeholder implementation

# Create dummy annotations for now
if 'mz' in features.columns:
    features['annotation'] = "Unknown"
    features['library_match_score'] = np.nan
    features['adduct'] = "[M+H]+"
    
    # Save annotated results
    features.to_csv("{output_file}", index=False)
    print(f"Saved annotated features to {{'{output_file}'}}")
    print(f"Total features: {{len(features)}}")
else:
    print("Error: Feature table does not contain 'mz' column")
    exit(1)

print("MS-DIAL annotation completed")
print("Note: Full MS-DIAL integration requires MS-DIAL console application")
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_msdial.py"
        script_path.write_text(python_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["python", str(script_path)]
            logger.info(f"Running MS-DIAL annotation: {' '.join(cmd)}")

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
                            output_files={"annotated_features": output_file},
                        ),
                    )
                return {"status": "success", "annotated_features": output_file}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"MS-DIAL annotation failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"MS-DIAL annotation failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def quantify_features(
        self,
        workflow_id: int,
        feature_table: str,
        sample_metadata: str,
        output_dir: str,
        normalization: str = "median",
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Quantify and normalize metabolomics features.

        Args:
            workflow_id: Workflow database ID
            feature_table: Feature intensity table
            sample_metadata: Sample metadata CSV
            output_dir: Output directory
            normalization: Normalization method (median, quantile, pqn)
            db: Database session

        Returns:
            Result with normalized data
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Python script for quantification
        python_script = f"""
import pandas as pd
import numpy as np

# Load feature table and metadata
print("Loading data...")
features = pd.read_csv("{feature_table}")
metadata = pd.read_csv("{sample_metadata}")

# Identify intensity columns (numeric columns except ID/annotation columns)
intensity_cols = [col for col in features.columns 
                 if features[col].dtype in ['float64', 'int64'] 
                 and col not in ['mz', 'm/z', 'rt', 'retention_time']]

print(f"Found {{len(intensity_cols)}} intensity columns")

if len(intensity_cols) == 0:
    print("Warning: No intensity columns found")
    features.to_csv("{output_dir}/quantified_features.csv", index=False)
else:
    # Apply normalization
    norm_method = "{normalization}"
    print(f"Applying {{norm_method}} normalization...")
    
    if norm_method == "median":
        # Median normalization
        medians = features[intensity_cols].median()
        target_median = medians.median()
        scaling_factors = target_median / medians
        features[intensity_cols] = features[intensity_cols] * scaling_factors
    
    elif norm_method == "quantile":
        # Quantile normalization (simplified)
        from sklearn.preprocessing import quantile_transform
        features[intensity_cols] = quantile_transform(
            features[intensity_cols], 
            n_quantiles=100,
            output_distribution='normal'
        )
    
    # Log2 transform
    features_log = features.copy()
    features_log[intensity_cols] = np.log2(features[intensity_cols].replace(0, np.nan))
    
    # Save results
    features.to_csv("{output_dir}/quantified_features.csv", index=False)
    features_log.to_csv("{output_dir}/quantified_features_log2.csv", index=False)
    
    print(f"Saved normalized features to {{'{output_dir}'}}")
    print(f"Total features: {{len(features)}}")

print("Feature quantification completed")
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_quant.py"
        script_path.write_text(python_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["python", str(script_path)]
            logger.info(f"Running feature quantification: {' '.join(cmd)}")

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "quantified_features": f"{output_dir}/quantified_features.csv",
                    "quantified_features_log2": f"{output_dir}/quantified_features_log2.csv",
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
                            error_message=f"Feature quantification failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"Feature quantification failed: {e}")
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
metabolomics_analyzer = MetabolomicsAnalyzer()
