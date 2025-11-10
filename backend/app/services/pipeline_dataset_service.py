"""
Pipeline Dataset Integration Service

Handles integration between Pipeline Builder and Dataset Manager:
- Use datasets as pipeline input sources
- Auto-create output datasets from pipeline results
- Record lineage for pipeline operations
- Link datasets to pipeline runs
"""

import hashlib
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from sqlalchemy.orm import Session

from app.models.dataset import Dataset, DatasetFileEntry, DatasetLineage
from app.models.run import Run
from app.models.datafile import DataFile
from app.services.dataset_manager import DatasetManager

logger = logging.getLogger(__name__)


class PipelineDatasetService:
    """Service for integrating datasets with pipeline execution."""

    def __init__(self, db: Session):
        self.db = db
        self.dataset_manager = DatasetManager(db)

    def get_dataset_files_for_input(
        self, dataset_id: int, file_role: Optional[str] = None
    ) -> List[Dict]:
        """
        Get files from a dataset suitable for pipeline input.

        Args:
            dataset_id: Dataset ID
            file_role: Filter by file role (e.g., 'input', 'primary', 'reference')

        Returns:
            List of file information dicts with path, name, type, role
        """
        dataset = self.dataset_manager.get_dataset(dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found")

        files = []
        for entry in dataset.files:
            if file_role and entry.role != file_role:
                continue

            files.append(
                {
                    "id": entry.id,
                    "path": entry.file_path,
                    "name": entry.file_name,
                    "type": entry.file_type,
                    "role": entry.role,
                    "size": entry.file_size,
                }
            )

        return files

    def create_dataset_from_run(
        self,
        run: Run,
        dataset_name: Optional[str] = None,
        dataset_description: Optional[str] = None,
        data_type: Optional[str] = None,
        file_roles: Optional[Dict[str, str]] = None,
        tags: Optional[List[str]] = None,
    ) -> Dataset:
        """
        Create a dataset from pipeline run outputs.

        Args:
            run: Pipeline run object
            dataset_name: Name for the dataset (defaults to run name)
            dataset_description: Dataset description
            data_type: Type of data (genomics, transcriptomics, etc.)
            file_roles: Dict mapping file IDs to roles (primary, secondary, etc.)
            tags: List of tag names to apply

        Returns:
            Created Dataset object
        """
        # Get output files from run
        output_file_ids = run.output_files or []
        if not output_file_ids:
            raise ValueError(f"Run {run.id} has no output files")

        # Create dataset
        name = dataset_name or f"{run.name}_output"
        description = (
            dataset_description or f"Output dataset from pipeline run: {run.name}"
        )

        # Infer data type from pipeline type if not provided
        if not data_type:
            data_type = self._infer_data_type_from_run(run)

        dataset = self.dataset_manager.create_dataset(
            project_id=run.project_id,
            name=name,
            description=description,
            data_type=data_type,
            created_by=run.owner_id,
        )

        # Add output files to dataset
        file_roles = file_roles or {}
        for file_id in output_file_ids:
            datafile = self.db.query(DataFile).filter(DataFile.id == file_id).first()
            if not datafile:
                logger.warning(f"Output file {file_id} not found, skipping")
                continue

            # Determine file role
            role = file_roles.get(str(file_id), "output")

            # Add file to dataset with hash computation
            try:
                self.dataset_manager.add_file(
                    dataset_id=dataset.id,
                    file_path=datafile.file_path,
                    file_name=datafile.filename,
                    file_type=datafile.file_type,
                    role=role,
                    compute_hash=True,
                )
            except Exception as e:
                logger.error(f"Failed to add file {datafile.filename} to dataset: {e}")

        # Add tags
        if tags:
            for tag_name in tags:
                try:
                    tag = self.dataset_manager.get_or_create_tag(tag_name)
                    self.dataset_manager.add_tag_to_dataset(dataset.id, tag.id)
                except Exception as e:
                    logger.error(f"Failed to add tag {tag_name}: {e}")

        # Add pipeline run as metadata
        dataset.metadata = dataset.metadata or {}
        dataset.metadata["pipeline_run"] = {
            "run_id": run.id,
            "run_name": run.name,
            "pipeline_type": run.pipeline_type,
            "started_at": run.started_at.isoformat() if run.started_at else None,
            "finished_at": run.finished_at.isoformat() if run.finished_at else None,
        }
        self.db.commit()

        logger.info(f"Created dataset {dataset.id} from run {run.id}")
        return dataset

    def record_pipeline_lineage(
        self,
        dataset_id: int,
        run: Run,
        input_file_ids: Optional[List[int]] = None,
        output_file_ids: Optional[List[int]] = None,
    ) -> List[DatasetLineage]:
        """
        Record lineage for pipeline operations.

        Args:
            dataset_id: Output dataset ID
            run: Pipeline run object
            input_file_ids: List of input DatasetFileEntry IDs
            output_file_ids: List of output DatasetFileEntry IDs

        Returns:
            List of created lineage records
        """
        lineage_records = []

        # Get pipeline configuration
        pipeline_config = run.pipeline_config or {}
        definition = pipeline_config.get("definition", {})
        nodes = definition.get("nodes", [])

        # If no specific file mapping provided, record general pipeline lineage
        if not input_file_ids or not output_file_ids:
            lineage = self.dataset_manager.add_lineage(
                dataset_id=dataset_id,
                operation_type="pipeline",
                operation_id=run.id,
                operation_params={
                    "pipeline_type": run.pipeline_type,
                    "parameters": run.parameters,
                    "tools": [node.get("data", {}).get("label") for node in nodes],
                },
                tool_name=run.name,
                tool_version=pipeline_config.get("version", "1.0"),
            )
            lineage_records.append(lineage)
            return lineage_records

        # Create lineage records for each input->output pair
        output_file_ids = output_file_ids or []
        for output_file_id in output_file_ids:
            for input_file_id in input_file_ids or []:
                try:
                    lineage = self.dataset_manager.add_lineage(
                        dataset_id=dataset_id,
                        operation_type="pipeline",
                        operation_id=run.id,
                        source_file_id=input_file_id,
                        output_file_id=output_file_id,
                        operation_params={
                            "pipeline_type": run.pipeline_type,
                            "parameters": run.parameters,
                            "node_count": len(nodes),
                        },
                        tool_name=run.name,
                        tool_version=pipeline_config.get("version", "1.0"),
                    )
                    lineage_records.append(lineage)
                except Exception as e:
                    logger.error(f"Failed to create lineage record: {e}")

        logger.info(
            f"Created {len(lineage_records)} lineage records for dataset {dataset_id}"
        )
        return lineage_records

    def link_input_datasets_to_run(
        self, run_id: int, dataset_ids: List[int], input_mapping: Optional[Dict] = None
    ) -> Run:
        """
        Link input datasets to a pipeline run.

        Args:
            run_id: Pipeline run ID
            dataset_ids: List of dataset IDs used as inputs
            input_mapping: Optional mapping of node inputs to dataset files

        Returns:
            Updated Run object
        """
        run = self.db.query(Run).filter(Run.id == run_id).first()
        if not run:
            raise ValueError(f"Run {run_id} not found")

        # Store dataset references in run metadata
        run.parameters = run.parameters or {}
        run.parameters["input_datasets"] = dataset_ids

        if input_mapping:
            run.input_mapping = run.input_mapping or {}
            run.input_mapping["dataset_files"] = input_mapping

        self.db.commit()
        logger.info(f"Linked {len(dataset_ids)} datasets to run {run_id}")
        return run

    def get_run_output_datasets(self, run_id: int) -> List[Dataset]:
        """
        Get all datasets created from a pipeline run.

        Args:
            run_id: Pipeline run ID

        Returns:
            List of Dataset objects
        """
        # Query datasets with run_id in metadata
        datasets = (
            self.db.query(Dataset)
            .filter(Dataset.metadata["pipeline_run"]["run_id"].astext == str(run_id))
            .all()
        )
        return datasets

    def get_dataset_usage_in_runs(self, dataset_id: int) -> List[Dict]:
        """
        Get all pipeline runs that used this dataset as input.

        Args:
            dataset_id: Dataset ID

        Returns:
            List of run information dicts
        """
        # Query runs with dataset_id in parameters
        runs = (
            self.db.query(Run)
            .filter(Run.parameters["input_datasets"].astext.contains(str(dataset_id)))
            .all()
        )

        return [
            {
                "run_id": run.id,
                "run_name": run.name,
                "status": run.status,
                "started_at": run.started_at.isoformat() if run.started_at else None,
                "finished_at": run.finished_at.isoformat() if run.finished_at else None,
            }
            for run in runs
        ]

    def auto_create_output_dataset(
        self, run: Run, auto_tags: bool = True
    ) -> Optional[Dataset]:
        """
        Automatically create output dataset for completed pipeline run.

        Args:
            run: Completed pipeline run
            auto_tags: Whether to automatically add tags based on pipeline type

        Returns:
            Created Dataset or None if no outputs
        """
        if run.status != "completed":
            logger.warning(f"Run {run.id} not completed, skipping dataset creation")
            return None

        if not run.output_files:
            logger.warning(f"Run {run.id} has no output files")
            return None

        # Determine tags
        tags = []
        if auto_tags:
            tags.append("pipeline-output")
            if run.pipeline_type:
                tags.append(f"pipeline-{run.pipeline_type}")

        # Infer data type
        data_type = self._infer_data_type_from_run(run)

        # Create dataset
        try:
            dataset = self.create_dataset_from_run(
                run=run,
                data_type=data_type,
                tags=tags,
            )

            # Record lineage
            self.record_pipeline_lineage(dataset_id=dataset.id, run=run)

            return dataset
        except Exception as e:
            logger.error(f"Failed to auto-create dataset for run {run.id}: {e}")
            return None

    def _infer_data_type_from_run(self, run: Run) -> str:
        """Infer data type from run parameters and configuration."""
        # Check run parameters for hints
        params = run.parameters or {}
        if "data_type" in params:
            return params["data_type"]

        # Check pipeline template ID
        template_id = run.pipeline_template_id or ""
        if "rna" in template_id.lower() or "transcriptomics" in template_id.lower():
            return "transcriptomics"
        elif "dna" in template_id.lower() or "genomics" in template_id.lower():
            return "genomics"
        elif "protein" in template_id.lower() or "proteomics" in template_id.lower():
            return "proteomics"
        elif "metabol" in template_id.lower():
            return "metabolomics"
        elif "single" in template_id.lower() and "cell" in template_id.lower():
            return "single-cell"

        # Check output file types
        output_file_ids = run.output_files or []
        if output_file_ids:
            datafile = (
                self.db.query(DataFile)
                .filter(DataFile.id == output_file_ids[0])
                .first()
            )
            if datafile:
                file_type = datafile.file_type.lower()
                if file_type in ["vcf", "bam", "sam", "fastq"]:
                    return "genomics"
                elif file_type in ["h5ad", "loom"]:
                    return "single-cell"
                elif file_type in ["csv", "tsv", "txt"]:
                    return "general"

        return "general"

    def prepare_dataset_files_for_pipeline(
        self, dataset_id: int, node_id: str, input_name: str
    ) -> List[str]:
        """
        Prepare dataset files for use in pipeline node input.

        Args:
            dataset_id: Dataset ID
            node_id: Pipeline node ID
            input_name: Input parameter name

        Returns:
            List of file paths ready for pipeline execution
        """
        files = self.get_dataset_files_for_input(dataset_id)
        return [f["path"] for f in files]

    def validate_dataset_for_pipeline_input(
        self, dataset_id: int, required_file_types: Optional[List[str]] = None
    ) -> Tuple[bool, List[str]]:
        """
        Validate that a dataset is suitable for pipeline input.

        Args:
            dataset_id: Dataset ID
            required_file_types: List of required file types (optional)

        Returns:
            Tuple of (is_valid, list of error messages)
        """
        errors = []

        # Check dataset exists and is active
        dataset = self.dataset_manager.get_dataset(dataset_id)
        if not dataset:
            return False, ["Dataset not found"]

        if dataset.status != "active":
            errors.append(f"Dataset status is {dataset.status}, must be active")

        # Check has files
        if not dataset.files:
            errors.append("Dataset has no files")
            return False, errors

        # Check file types if specified
        if required_file_types:
            available_types = {f.file_type.lower() for f in dataset.files}
            required_types = {t.lower() for t in required_file_types}
            missing_types = required_types - available_types
            if missing_types:
                errors.append(
                    f"Missing required file types: {', '.join(missing_types)}"
                )

        # Check file paths exist
        for file_entry in dataset.files:
            file_path = Path(file_entry.file_path)
            if not file_path.exists():
                errors.append(f"File not found: {file_entry.file_name}")

        is_valid = len(errors) == 0
        return is_valid, errors


def get_pipeline_dataset_service(db: Session) -> PipelineDatasetService:
    """Get PipelineDatasetService instance."""
    return PipelineDatasetService(db)
