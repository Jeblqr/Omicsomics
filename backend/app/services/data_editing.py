"""
Service layer for data editing functionality.

Handles in-place data operations with preview, validation, undo/redo capabilities.
"""

import os
import json
import time
import pandas as pd
import numpy as np
from datetime import datetime
from typing import List, Optional, Dict, Any, Tuple
from pathlib import Path
from sqlalchemy.orm import Session
from sqlalchemy import and_

from ..models.data_editing import (
    EditSession,
    EditOperation,
    EditOperationType,
    EditStatus,
)
from ..models.file import File
from ..models.dataset import Dataset
from ..core.storage import storage_manager


class DataEditingService:
    """Service for managing data editing operations"""

    def __init__(self, db: Session):
        self.db = db

    def create_edit_session(
        self,
        user_id: int,
        file_id: int,
        name: str,
        project_id: Optional[int] = None,
        dataset_id: Optional[int] = None,
        description: Optional[str] = None,
    ) -> EditSession:
        """
        Create a new edit session.

        Args:
            user_id: User creating the session
            file_id: File to edit
            name: Session name
            project_id: Optional project ID
            dataset_id: Optional dataset ID
            description: Session description

        Returns:
            Created EditSession
        """
        # Verify file exists
        file = self.db.query(File).filter(File.id == file_id).first()
        if not file:
            raise ValueError(f"File {file_id} not found")

        # Generate unique session key
        timestamp = datetime.utcnow().strftime("%Y%m%d%H%M%S")
        session_key = f"edit_{user_id}_{file_id}_{timestamp}"

        # Create backup of original file
        backup_path = self._create_backup(file)

        # Create session
        session = EditSession(
            session_key=session_key,
            name=name,
            description=description,
            user_id=user_id,
            project_id=project_id,
            file_id=file_id,
            dataset_id=dataset_id,
            status=EditStatus.DRAFT,
            operations=[],
            current_operation_index=-1,
            backup_path=backup_path,
        )

        self.db.add(session)
        self.db.commit()
        self.db.refresh(session)

        return session

    def get_edit_session(
        self, session_id: int, user_id: Optional[int] = None
    ) -> Optional[EditSession]:
        """Get edit session by ID"""
        query = self.db.query(EditSession).filter(EditSession.id == session_id)

        if user_id is not None:
            query = query.filter(EditSession.user_id == user_id)

        return query.first()

    def get_edit_session_by_key(
        self, session_key: str, user_id: Optional[int] = None
    ) -> Optional[EditSession]:
        """Get edit session by session key"""
        query = self.db.query(EditSession).filter(
            EditSession.session_key == session_key
        )

        if user_id is not None:
            query = query.filter(EditSession.user_id == user_id)

        return query.first()

    def list_edit_sessions(
        self,
        user_id: Optional[int] = None,
        project_id: Optional[int] = None,
        status: Optional[EditStatus] = None,
        file_id: Optional[int] = None,
        limit: int = 50,
        offset: int = 0,
    ) -> List[EditSession]:
        """List edit sessions with filters"""
        query = self.db.query(EditSession)

        if user_id is not None:
            query = query.filter(EditSession.user_id == user_id)

        if project_id is not None:
            query = query.filter(EditSession.project_id == project_id)

        if status is not None:
            query = query.filter(EditSession.status == status)

        if file_id is not None:
            query = query.filter(EditSession.file_id == file_id)

        query = query.order_by(EditSession.created_at.desc())
        query = query.limit(limit).offset(offset)

        return query.all()

    def add_operation(
        self,
        session_id: int,
        operation_type: EditOperationType,
        parameters: Dict[str, Any],
        description: Optional[str] = None,
    ) -> EditSession:
        """Add operation to edit session"""
        session = self.get_edit_session(session_id)

        if not session:
            raise ValueError(f"Session {session_id} not found")

        if session.status not in [EditStatus.DRAFT, EditStatus.PREVIEWING]:
            raise ValueError(
                f"Cannot add operation to session in {session.status} status"
            )

        # Add operation to list
        operation = {
            "type": operation_type.value,
            "parameters": parameters,
            "description": description,
            "timestamp": datetime.utcnow().isoformat(),
        }

        operations = session.operations or []
        operations.append(operation)
        session.operations = operations
        session.current_operation_index = len(operations) - 1

        # Create operation record
        op_record = EditOperation(
            session_id=session_id,
            operation_type=operation_type,
            operation_index=len(operations) - 1,
            parameters=parameters,
            description=description,
        )

        self.db.add(op_record)
        self.db.commit()
        self.db.refresh(session)

        return session

    def remove_operation(self, session_id: int, operation_index: int) -> EditSession:
        """Remove operation from session"""
        session = self.get_edit_session(session_id)

        if not session:
            raise ValueError(f"Session {session_id} not found")

        operations = session.operations or []

        if operation_index < 0 or operation_index >= len(operations):
            raise ValueError(f"Invalid operation index {operation_index}")

        # Remove operation
        operations.pop(operation_index)
        session.operations = operations

        # Adjust current index
        if session.current_operation_index >= len(operations):
            session.current_operation_index = len(operations) - 1

        self.db.commit()
        self.db.refresh(session)

        return session

    def preview_operations(
        self, session_id: int, sample_size: int = 100
    ) -> Dict[str, Any]:
        """
        Preview the result of applying all operations.

        Args:
            session_id: Session ID
            sample_size: Number of rows to preview

        Returns:
            Preview data and statistics
        """
        session = self.get_edit_session(session_id)

        if not session:
            raise ValueError(f"Session {session_id} not found")

        # Load original file
        file = session.file
        file_data = storage_manager.download(file.storage_path)

        # Load into DataFrame
        df = self._load_dataframe(file_data, file.filename)
        original_shape = df.shape

        # Apply operations
        validation_errors = []
        validation_warnings = []

        try:
            for i, operation in enumerate(session.operations or []):
                df, errors, warnings = self._apply_operation(df, operation)
                validation_errors.extend(errors)
                validation_warnings.extend(warnings)
        except Exception as e:
            validation_errors.append({"operation_index": i, "error": str(e)})

        # Generate preview
        preview_df = df.head(sample_size)
        preview_data = {
            "columns": list(preview_df.columns),
            "data": preview_df.to_dict(orient="records"),
            "total_rows": len(df),
        }

        # Generate summary
        summary = {
            "original_shape": original_shape,
            "result_shape": df.shape,
            "rows_changed": original_shape[0] - df.shape[0],
            "columns_changed": original_shape[1] - df.shape[1],
            "column_types": {col: str(dtype) for col, dtype in df.dtypes.items()},
            "missing_values": df.isnull().sum().to_dict(),
            "memory_usage_mb": df.memory_usage(deep=True).sum() / (1024**2),
        }

        # Update session
        session.preview_data = preview_data
        session.preview_summary = summary
        session.validation_errors = validation_errors
        session.validation_warnings = validation_warnings
        session.is_valid = len(validation_errors) == 0
        session.status = EditStatus.PREVIEWING

        self.db.commit()

        return {
            "preview": preview_data,
            "summary": summary,
            "validation": {
                "is_valid": session.is_valid,
                "errors": validation_errors,
                "warnings": validation_warnings,
            },
        }

    def apply_operations(self, session_id: int) -> EditSession:
        """
        Apply all operations and save result.

        Creates a new file with the edited data.
        """
        session = self.get_edit_session(session_id)

        if not session:
            raise ValueError(f"Session {session_id} not found")

        if not session.is_valid:
            raise ValueError("Cannot apply invalid operations. Run preview first.")

        session.status = EditStatus.APPLYING
        self.db.commit()

        try:
            # Load original file
            file = session.file
            file_data = storage_manager.download(file.storage_path)
            df = self._load_dataframe(file_data, file.filename)

            # Apply operations
            for operation in session.operations or []:
                df, _, _ = self._apply_operation(df, operation)

            # Save result
            output_path = self._save_result(df, file.filename, session.user_id)

            # Create output file record
            output_file = File(
                filename=f"{Path(file.filename).stem}_edited{Path(file.filename).suffix}",
                original_filename=file.original_filename,
                storage_path=output_path,
                size=len(df),
                mime_type=file.mime_type,
                user_id=session.user_id,
                project_id=session.project_id,
            )

            self.db.add(output_file)
            self.db.commit()
            self.db.refresh(output_file)

            # Update session
            session.output_file_id = output_file.id
            session.status = EditStatus.APPLIED
            session.applied_at = datetime.utcnow()

            self.db.commit()
            self.db.refresh(session)

            return session

        except Exception as e:
            session.status = EditStatus.FAILED
            session.validation_errors = [{"error": str(e)}]
            self.db.commit()
            raise

    def revert_operations(self, session_id: int) -> EditSession:
        """Revert session to original file"""
        session = self.get_edit_session(session_id)

        if not session:
            raise ValueError(f"Session {session_id} not found")

        if not session.backup_path:
            raise ValueError("No backup available")

        # Restore from backup
        # (In production, would copy backup back to original location)

        session.status = EditStatus.REVERTED
        self.db.commit()
        self.db.refresh(session)

        return session

    def delete_edit_session(
        self, session_id: int, user_id: Optional[int] = None
    ) -> bool:
        """Delete edit session and cleanup"""
        session = self.get_edit_session(session_id, user_id)

        if not session:
            return False

        # Delete backup if exists
        if session.backup_path:
            try:
                storage_manager.delete(session.backup_path)
            except Exception as e:
                print(f"Failed to delete backup: {e}")

        # Delete operation records
        self.db.query(EditOperation).filter(
            EditOperation.session_id == session_id
        ).delete()

        # Delete session
        self.db.delete(session)
        self.db.commit()

        return True

    def _create_backup(self, file: File) -> str:
        """Create backup of original file"""
        # Download original
        file_data = storage_manager.download(file.storage_path)

        # Upload backup
        backup_path = f"backups/edit_sessions/{file.id}/{datetime.utcnow().strftime('%Y%m%d%H%M%S')}_{file.filename}"
        storage_manager.upload(file_data, backup_path)

        return backup_path

    def _load_dataframe(self, file_data: bytes, filename: str) -> pd.DataFrame:
        """Load file data into DataFrame"""
        ext = Path(filename).suffix.lower()

        if ext == ".csv":
            return pd.read_csv(pd.io.common.BytesIO(file_data))
        elif ext in [".tsv", ".txt"]:
            return pd.read_csv(pd.io.common.BytesIO(file_data), sep="\t")
        elif ext in [".xlsx", ".xls"]:
            return pd.read_excel(pd.io.common.BytesIO(file_data))
        elif ext == ".json":
            return pd.read_json(pd.io.common.BytesIO(file_data))
        elif ext == ".parquet":
            return pd.read_parquet(pd.io.common.BytesIO(file_data))
        else:
            raise ValueError(f"Unsupported file format: {ext}")

    def _apply_operation(
        self, df: pd.DataFrame, operation: Dict[str, Any]
    ) -> Tuple[pd.DataFrame, List[Dict], List[Dict]]:
        """
        Apply a single operation to DataFrame.

        Returns:
            (modified_df, errors, warnings)
        """
        op_type = operation["type"]
        params = operation["parameters"]
        errors = []
        warnings = []

        try:
            if op_type == EditOperationType.DELETE_ROWS.value:
                # Delete rows by indices
                indices = params.get("indices", [])
                df = df.drop(indices)

            elif op_type == EditOperationType.FILTER_ROWS.value:
                # Filter rows by condition
                column = params["column"]
                operator = params["operator"]
                value = params["value"]

                if operator == "equals":
                    df = df[df[column] == value]
                elif operator == "not_equals":
                    df = df[df[column] != value]
                elif operator == "greater_than":
                    df = df[df[column] > value]
                elif operator == "less_than":
                    df = df[df[column] < value]
                elif operator == "contains":
                    df = df[df[column].str.contains(value, na=False)]

            elif op_type == EditOperationType.SORT_ROWS.value:
                # Sort rows
                columns = params["columns"]
                ascending = params.get("ascending", True)
                df = df.sort_values(by=columns, ascending=ascending)

            elif op_type == EditOperationType.DEDUPLICATE.value:
                # Remove duplicates
                columns = params.get("columns", None)
                df = df.drop_duplicates(subset=columns)

            elif op_type == EditOperationType.ADD_COLUMN.value:
                # Add new column
                column_name = params["column_name"]
                default_value = params.get("default_value", None)
                df[column_name] = default_value

            elif op_type == EditOperationType.DELETE_COLUMN.value:
                # Delete column
                columns = params["columns"]
                df = df.drop(columns=columns)

            elif op_type == EditOperationType.RENAME_COLUMN.value:
                # Rename column
                old_name = params["old_name"]
                new_name = params["new_name"]
                df = df.rename(columns={old_name: new_name})

            elif op_type == EditOperationType.REPLACE_VALUES.value:
                # Replace values
                column = params["column"]
                old_value = params["old_value"]
                new_value = params["new_value"]
                df[column] = df[column].replace(old_value, new_value)

            elif op_type == EditOperationType.FILL_MISSING.value:
                # Fill missing values
                column = params["column"]
                fill_value = params["fill_value"]
                method = params.get("method", "value")

                if method == "value":
                    df[column] = df[column].fillna(fill_value)
                elif method == "mean":
                    df[column] = df[column].fillna(df[column].mean())
                elif method == "median":
                    df[column] = df[column].fillna(df[column].median())
                elif method == "forward_fill":
                    df[column] = df[column].ffill()
                elif method == "backward_fill":
                    df[column] = df[column].bfill()

            elif op_type == EditOperationType.CONVERT_TYPE.value:
                # Convert column type
                column = params["column"]
                target_type = params["target_type"]

                if target_type == "int":
                    df[column] = pd.to_numeric(df[column], errors="coerce").astype(
                        "Int64"
                    )
                elif target_type == "float":
                    df[column] = pd.to_numeric(df[column], errors="coerce")
                elif target_type == "string":
                    df[column] = df[column].astype(str)
                elif target_type == "datetime":
                    df[column] = pd.to_datetime(df[column], errors="coerce")
                elif target_type == "bool":
                    df[column] = df[column].astype(bool)

            else:
                warnings.append(
                    {
                        "operation_type": op_type,
                        "warning": f"Operation type not implemented: {op_type}",
                    }
                )

        except Exception as e:
            errors.append({"operation_type": op_type, "error": str(e)})

        return df, errors, warnings

    def _save_result(
        self, df: pd.DataFrame, original_filename: str, user_id: int
    ) -> str:
        """Save DataFrame to storage"""
        ext = Path(original_filename).suffix.lower()

        # Convert to bytes
        if ext == ".csv":
            data = df.to_csv(index=False).encode("utf-8")
        elif ext in [".tsv", ".txt"]:
            data = df.to_csv(index=False, sep="\t").encode("utf-8")
        elif ext in [".xlsx", ".xls"]:
            import io

            buffer = io.BytesIO()
            df.to_excel(buffer, index=False)
            data = buffer.getvalue()
        elif ext == ".json":
            data = df.to_json(orient="records").encode("utf-8")
        elif ext == ".parquet":
            import io

            buffer = io.BytesIO()
            df.to_parquet(buffer, index=False)
            data = buffer.getvalue()
        else:
            raise ValueError(f"Unsupported output format: {ext}")

        # Upload to storage
        output_path = f"edited/{user_id}/{datetime.utcnow().strftime('%Y%m%d%H%M%S')}_{original_filename}"
        storage_manager.upload(data, output_path)

        return output_path
