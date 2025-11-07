"""Workflow executor for running bioinformatics pipelines."""

import asyncio
import logging
import subprocess
from pathlib import Path
from typing import Any

from app.models.workflow import WorkflowStatus
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)


class WorkflowExecutor:
    """Execute bioinformatics workflows."""

    def __init__(self, work_dir: Path | None = None):
        """Initialize workflow executor."""
        self.work_dir = work_dir or Path("/tmp/omicsomics_workflows")
        self.work_dir.mkdir(parents=True, exist_ok=True)

    async def execute_nextflow(
        self,
        workflow_id: int,
        pipeline: str,
        params: dict[str, Any],
        db: AsyncSession,
    ) -> dict[str, Any]:
        """
        Execute a Nextflow pipeline.

        Args:
            workflow_id: Database workflow ID
            pipeline: Nextflow pipeline path or URL
            params: Pipeline parameters
            db: Database session

        Returns:
            Execution result with status and outputs
        """
        work_dir = self.work_dir / f"workflow_{workflow_id}"
        work_dir.mkdir(exist_ok=True)

        # Build Nextflow command
        cmd = ["nextflow", "run", pipeline]
        for key, value in params.items():
            cmd.extend([f"--{key}", str(value)])
        cmd.extend(["-work-dir", str(work_dir)])

        try:
            # Update status to running
            await workflow_service.update_workflow(
                db,
                workflow_id,
                workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
            )

            # Execute pipeline
            logger.info(f"Executing Nextflow: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=str(work_dir),
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                # Success
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.COMPLETED, logs=logs
                    ),
                )
                return {"status": "success", "logs": logs}
            else:
                # Failure
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED,
                        logs=logs,
                        error_message=f"Pipeline failed with exit code {process.returncode}",
                    ),
                )
                return {
                    "status": "failed",
                    "error": f"Exit code {process.returncode}",
                    "logs": logs,
                }

        except Exception as e:
            logger.error(f"Workflow {workflow_id} failed: {e}")
            await workflow_service.update_workflow(
                db,
                workflow_id,
                workflow_schema.WorkflowUpdate(
                    status=WorkflowStatus.FAILED, error_message=str(e)
                ),
            )
            return {"status": "error", "error": str(e)}

    async def execute_fastqc(
        self, workflow_id: int, input_files: list[str], output_dir: str, db: AsyncSession
    ) -> dict[str, Any]:
        """
        Execute FastQC quality control.

        Args:
            workflow_id: Database workflow ID
            input_files: List of FASTQ file paths
            output_dir: Output directory for QC results
            db: Database session

        Returns:
            Execution result
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        cmd = ["fastqc", "-o", output_dir, "-t", "4"] + input_files

        try:
            await workflow_service.update_workflow(
                db,
                workflow_id,
                workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
            )

            logger.info(f"Executing FastQC: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.COMPLETED,
                        logs=logs,
                        output_files={"output_dir": output_dir},
                    ),
                )
                return {"status": "success", "output_dir": output_dir, "logs": logs}
            else:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED,
                        logs=logs,
                        error_message=f"FastQC failed with exit code {process.returncode}",
                    ),
                )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except FileNotFoundError:
            error_msg = "FastQC not found. Please install FastQC."
            await workflow_service.update_workflow(
                db,
                workflow_id,
                workflow_schema.WorkflowUpdate(
                    status=WorkflowStatus.FAILED, error_message=error_msg
                ),
            )
            return {"status": "error", "error": error_msg}
        except Exception as e:
            logger.error(f"FastQC workflow {workflow_id} failed: {e}")
            await workflow_service.update_workflow(
                db,
                workflow_id,
                workflow_schema.WorkflowUpdate(
                    status=WorkflowStatus.FAILED, error_message=str(e)
                ),
            )
            return {"status": "error", "error": str(e)}


# Global executor instance
workflow_executor = WorkflowExecutor()
