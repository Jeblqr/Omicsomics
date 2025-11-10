"""
Custom Script Tools Service

Provides functionality for managing and executing custom user scripts.
Supports Python, R, and Bash scripts with parameter validation,
sandboxed execution, and result tracking.
"""

import asyncio
import json
import logging
import os
import subprocess
import tempfile
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import jsonschema
from sqlalchemy import and_, or_
from sqlalchemy.orm import Session

from ..core.config import settings
from ..models.custom_scripts import (
    CustomScript,
    ScriptExecution,
    ScriptLanguage,
    ScriptStatus,
    ScriptVisibility,
)
from ..models.file import File
from ..services.minio_service import MinioService

logger = logging.getLogger(__name__)


class CustomScriptService:
    """Service for managing custom scripts and their execution."""

    def __init__(self, db: Session):
        self.db = db
        self.minio_service = MinioService()

    def create_script(
        self,
        user_id: int,
        name: str,
        description: str,
        language: ScriptLanguage,
        script_content: str,
        parameters_schema: Optional[Dict[str, Any]] = None,
        entry_point: Optional[str] = None,
        requirements: Optional[List[str]] = None,
        timeout: int = 300,
        max_memory: Optional[int] = None,
        visibility: ScriptVisibility = ScriptVisibility.PRIVATE,
        category: Optional[str] = None,
        tags: Optional[List[str]] = None,
    ) -> CustomScript:
        """
        Create a new custom script.

        Args:
            user_id: Owner user ID
            name: Script name
            description: Script description
            language: Script language (python, r, bash)
            script_content: The actual script code
            parameters_schema: JSON Schema for parameters validation
            entry_point: Entry point function/method (for Python/R)
            requirements: List of dependencies (pip packages, R packages, etc.)
            timeout: Execution timeout in seconds (default: 300)
            max_memory: Maximum memory in MB (optional)
            visibility: Script visibility (private, project, public)
            category: Script category (optional)
            tags: Script tags (optional)

        Returns:
            Created CustomScript object
        """
        # Generate unique script key
        script_key = self._generate_script_key(name)

        # Validate parameters schema if provided
        if parameters_schema:
            self._validate_json_schema(parameters_schema)

        # Create script record
        script = CustomScript(
            script_key=script_key,
            name=name,
            description=description,
            user_id=user_id,
            language=language,
            script_content=script_content,
            entry_point=entry_point,
            parameters_schema=parameters_schema,
            requirements=requirements or [],
            timeout=timeout,
            max_memory=max_memory,
            visibility=visibility,
            category=category,
            tags=tags or [],
        )

        self.db.add(script)
        self.db.commit()
        self.db.refresh(script)

        logger.info(f"Created custom script: {script.name} (ID: {script.id})")
        return script

    def get_script(
        self, script_id: int, user_id: Optional[int] = None
    ) -> Optional[CustomScript]:
        """
        Get a script by ID with visibility check.

        Args:
            script_id: Script ID
            user_id: Current user ID (for visibility check)

        Returns:
            CustomScript object or None
        """
        query = self.db.query(CustomScript).filter(CustomScript.id == script_id)

        # Apply visibility filter
        if user_id:
            query = query.filter(
                or_(
                    CustomScript.user_id == user_id,
                    CustomScript.visibility == ScriptVisibility.PUBLIC,
                )
            )

        return query.first()

    def list_scripts(
        self,
        user_id: Optional[int] = None,
        language: Optional[ScriptLanguage] = None,
        category: Optional[str] = None,
        visibility: Optional[ScriptVisibility] = None,
        tags: Optional[List[str]] = None,
        search: Optional[str] = None,
        verified_only: bool = False,
        skip: int = 0,
        limit: int = 50,
    ) -> Tuple[List[CustomScript], int]:
        """
        List scripts with filters.

        Args:
            user_id: Filter by owner or accessible scripts
            language: Filter by language
            category: Filter by category
            visibility: Filter by visibility
            tags: Filter by tags (any match)
            search: Search in name and description
            verified_only: Show only verified scripts
            skip: Number of records to skip
            limit: Maximum number of records to return

        Returns:
            Tuple of (list of scripts, total count)
        """
        query = self.db.query(CustomScript)

        # Visibility filter
        if user_id:
            query = query.filter(
                or_(
                    CustomScript.user_id == user_id,
                    CustomScript.visibility == ScriptVisibility.PUBLIC,
                )
            )
        else:
            query = query.filter(CustomScript.visibility == ScriptVisibility.PUBLIC)

        # Language filter
        if language:
            query = query.filter(CustomScript.language == language)

        # Category filter
        if category:
            query = query.filter(CustomScript.category == category)

        # Visibility filter
        if visibility:
            query = query.filter(CustomScript.visibility == visibility)

        # Tags filter
        if tags:
            for tag in tags:
                query = query.filter(CustomScript.tags.contains([tag]))

        # Search filter
        if search:
            search_term = f"%{search}%"
            query = query.filter(
                or_(
                    CustomScript.name.ilike(search_term),
                    CustomScript.description.ilike(search_term),
                )
            )

        # Verified filter
        if verified_only:
            query = query.filter(CustomScript.is_verified == True)

        # Get total count
        total = query.count()

        # Apply pagination and ordering
        scripts = (
            query.order_by(CustomScript.created_at.desc())
            .offset(skip)
            .limit(limit)
            .all()
        )

        return scripts, total

    def update_script(
        self,
        script_id: int,
        user_id: int,
        **update_fields,
    ) -> Optional[CustomScript]:
        """
        Update a script (owner only).

        Args:
            script_id: Script ID
            user_id: Current user ID
            **update_fields: Fields to update

        Returns:
            Updated CustomScript or None
        """
        script = (
            self.db.query(CustomScript)
            .filter(
                CustomScript.id == script_id,
                CustomScript.user_id == user_id,
            )
            .first()
        )

        if not script:
            return None

        # Validate parameters schema if being updated
        if "parameters_schema" in update_fields:
            self._validate_json_schema(update_fields["parameters_schema"])

        # Update fields
        for key, value in update_fields.items():
            if hasattr(script, key) and value is not None:
                setattr(script, key, value)

        script.updated_at = datetime.utcnow()
        self.db.commit()
        self.db.refresh(script)

        logger.info(f"Updated script: {script.name} (ID: {script.id})")
        return script

    def delete_script(self, script_id: int, user_id: int) -> bool:
        """
        Delete a script (owner only).

        Args:
            script_id: Script ID
            user_id: Current user ID

        Returns:
            True if deleted, False otherwise
        """
        script = (
            self.db.query(CustomScript)
            .filter(
                CustomScript.id == script_id,
                CustomScript.user_id == user_id,
            )
            .first()
        )

        if not script:
            return False

        # Delete associated executions
        self.db.query(ScriptExecution).filter(
            ScriptExecution.script_id == script_id
        ).delete()

        self.db.delete(script)
        self.db.commit()

        logger.info(f"Deleted script: {script.name} (ID: {script.id})")
        return True

    def execute_script(
        self,
        script_id: int,
        user_id: int,
        parameters: Optional[Dict[str, Any]] = None,
        input_file_ids: Optional[List[int]] = None,
        description: Optional[str] = None,
    ) -> ScriptExecution:
        """
        Execute a script with given parameters.

        Args:
            script_id: Script ID
            user_id: User ID
            parameters: Script parameters
            input_file_ids: Input file IDs
            description: Execution description

        Returns:
            ScriptExecution object
        """
        # Get script
        script = self.get_script(script_id, user_id)
        if not script:
            raise ValueError(f"Script not found: {script_id}")

        # Validate parameters
        if script.parameters_schema:
            self._validate_parameters(parameters or {}, script.parameters_schema)

        # Create execution record
        execution_key = self._generate_execution_key(script.name)
        execution = ScriptExecution(
            execution_key=execution_key,
            script_id=script_id,
            user_id=user_id,
            parameters=parameters or {},
            input_file_ids=input_file_ids or [],
            status=ScriptStatus.PENDING,
            description=description,
        )

        self.db.add(execution)
        self.db.commit()
        self.db.refresh(execution)

        # Execute script asynchronously (in real implementation, use Celery)
        try:
            self._run_script_execution(execution, script)
        except Exception as e:
            logger.error(f"Script execution failed: {e}", exc_info=True)
            execution.status = ScriptStatus.FAILED
            execution.error_text = str(e)
            execution.completed_at = datetime.utcnow()
            self.db.commit()

        return execution

    def get_execution(
        self, execution_id: int, user_id: Optional[int] = None
    ) -> Optional[ScriptExecution]:
        """Get execution by ID."""
        query = self.db.query(ScriptExecution).filter(
            ScriptExecution.id == execution_id
        )
        if user_id:
            query = query.filter(ScriptExecution.user_id == user_id)
        return query.first()

    def list_executions(
        self,
        script_id: Optional[int] = None,
        user_id: Optional[int] = None,
        status: Optional[ScriptStatus] = None,
        skip: int = 0,
        limit: int = 50,
    ) -> Tuple[List[ScriptExecution], int]:
        """List executions with filters."""
        query = self.db.query(ScriptExecution)

        if script_id:
            query = query.filter(ScriptExecution.script_id == script_id)
        if user_id:
            query = query.filter(ScriptExecution.user_id == user_id)
        if status:
            query = query.filter(ScriptExecution.status == status)

        total = query.count()
        executions = (
            query.order_by(ScriptExecution.created_at.desc())
            .offset(skip)
            .limit(limit)
            .all()
        )

        return executions, total

    def cancel_execution(self, execution_id: int, user_id: int) -> bool:
        """Cancel a running execution."""
        execution = (
            self.db.query(ScriptExecution)
            .filter(
                ScriptExecution.id == execution_id,
                ScriptExecution.user_id == user_id,
            )
            .first()
        )

        if not execution or execution.status not in [
            ScriptStatus.PENDING,
            ScriptStatus.RUNNING,
        ]:
            return False

        execution.status = ScriptStatus.CANCELLED
        execution.completed_at = datetime.utcnow()
        self.db.commit()

        logger.info(f"Cancelled execution: {execution.execution_key}")
        return True

    # Helper methods

    def _generate_script_key(self, name: str) -> str:
        """Generate unique script key."""
        timestamp = int(time.time() * 1000)
        clean_name = "".join(c if c.isalnum() else "_" for c in name.lower())[:30]
        return f"script_{clean_name}_{timestamp}"

    def _generate_execution_key(self, script_name: str) -> str:
        """Generate unique execution key."""
        timestamp = int(time.time() * 1000)
        clean_name = "".join(c if c.isalnum() else "_" for c in script_name.lower())[
            :20
        ]
        return f"exec_{clean_name}_{timestamp}"

    def _validate_json_schema(self, schema: Dict[str, Any]) -> None:
        """Validate that schema is a valid JSON Schema."""
        try:
            # Try to create a validator to check schema validity
            jsonschema.Draft7Validator.check_schema(schema)
        except jsonschema.SchemaError as e:
            raise ValueError(f"Invalid JSON Schema: {e}")

    def _validate_parameters(
        self, parameters: Dict[str, Any], schema: Dict[str, Any]
    ) -> None:
        """Validate parameters against JSON Schema."""
        try:
            jsonschema.validate(instance=parameters, schema=schema)
        except jsonschema.ValidationError as e:
            raise ValueError(f"Invalid parameters: {e.message}")

    def _run_script_execution(
        self, execution: ScriptExecution, script: CustomScript
    ) -> None:
        """
        Run script execution in sandboxed environment.

        Note: This is a simplified implementation. In production, you should:
        - Use Celery for async execution
        - Use Docker containers for sandboxing
        - Implement proper resource limits (CPU, memory, disk)
        - Handle file I/O properly
        - Implement proper logging and monitoring
        """
        execution.status = ScriptStatus.RUNNING
        execution.started_at = datetime.utcnow()
        self.db.commit()

        start_time = time.time()
        temp_dir = None

        try:
            # Create temporary directory
            temp_dir = tempfile.mkdtemp(prefix=f"script_{execution.id}_")

            # Download input files
            input_files = []
            if execution.input_file_ids:
                files = (
                    self.db.query(File)
                    .filter(File.id.in_(execution.input_file_ids))
                    .all()
                )
                for file in files:
                    local_path = os.path.join(temp_dir, file.filename)
                    self.minio_service.download_file(file.minio_key, local_path)
                    input_files.append(local_path)

            # Create script file
            script_ext = {
                ScriptLanguage.PYTHON: ".py",
                ScriptLanguage.R: ".R",
                ScriptLanguage.BASH: ".sh",
            }
            script_file = os.path.join(temp_dir, f"script{script_ext[script.language]}")
            with open(script_file, "w") as f:
                f.write(script.script_content)

            # Prepare command
            if script.language == ScriptLanguage.PYTHON:
                cmd = ["python3", script_file]
            elif script.language == ScriptLanguage.R:
                cmd = ["Rscript", script_file]
            else:  # BASH
                cmd = ["bash", script_file]

            # Add parameters as JSON file
            if execution.parameters:
                params_file = os.path.join(temp_dir, "parameters.json")
                with open(params_file, "w") as f:
                    json.dump(execution.parameters, f)
                cmd.extend(["--params", params_file])

            # Run script
            result = subprocess.run(
                cmd,
                cwd=temp_dir,
                capture_output=True,
                text=True,
                timeout=script.timeout,
            )

            # Process results
            execution.output_text = result.stdout
            execution.error_text = result.stderr if result.returncode != 0 else None
            execution.exit_code = result.returncode

            # Check for output files
            output_files = []
            output_dir = os.path.join(temp_dir, "output")
            if os.path.exists(output_dir):
                for filename in os.listdir(output_dir):
                    file_path = os.path.join(output_dir, filename)
                    if os.path.isfile(file_path):
                        # Upload to MinIO
                        file_key = (
                            f"script_outputs/{execution.execution_key}/{filename}"
                        )
                        self.minio_service.upload_file(file_path, file_key)

                        # Create File record
                        file_size = os.path.getsize(file_path)
                        file_record = File(
                            filename=filename,
                            file_type="application/octet-stream",
                            file_size=file_size,
                            minio_key=file_key,
                            user_id=execution.user_id,
                        )
                        self.db.add(file_record)
                        self.db.flush()
                        output_files.append(file_record.id)

            execution.output_file_ids = output_files

            # Update status
            if result.returncode == 0:
                execution.status = ScriptStatus.COMPLETED
            else:
                execution.status = ScriptStatus.FAILED

            # Update statistics
            script.total_executions += 1
            if execution.status == ScriptStatus.COMPLETED:
                script.successful_executions += 1
            else:
                script.failed_executions += 1
            script.last_executed_at = datetime.utcnow()

        except subprocess.TimeoutExpired:
            execution.status = ScriptStatus.FAILED
            execution.error_text = f"Execution timeout after {script.timeout} seconds"
            logger.error(f"Script execution timeout: {execution.execution_key}")

        except Exception as e:
            execution.status = ScriptStatus.FAILED
            execution.error_text = str(e)
            logger.error(f"Script execution error: {e}", exc_info=True)

        finally:
            # Calculate metrics
            execution.duration = time.time() - start_time
            execution.completed_at = datetime.utcnow()

            # Cleanup temporary directory
            if temp_dir and os.path.exists(temp_dir):
                import shutil

                shutil.rmtree(temp_dir, ignore_errors=True)

            self.db.commit()
            self.db.refresh(execution)

        logger.info(
            f"Script execution completed: {execution.execution_key} "
            f"(status: {execution.status}, duration: {execution.duration:.2f}s)"
        )
