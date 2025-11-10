"""
Celery tasks for file processing
"""

from celery import current_task
from app.celery_app import celery_app
from app.services.file_processor import FileProcessor
from app.database import AsyncSessionMaker
from app.models.datafile import DataFile
from sqlalchemy.orm.attributes import flag_modified
import logging
import asyncio

logger = logging.getLogger(__name__)


def run_async(coro):
    """Helper to run async code in sync Celery task"""
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    try:
        return loop.run_until_complete(coro)
    finally:
        loop.close()


@celery_app.task(bind=True, name="app.tasks.file_tasks.process_file_async")
def process_file_async(self, raw_file_id: int, sample_id: str = None):
    """
    Asynchronously process a file and convert to unified format

    Args:
        raw_file_id: ID of the raw file to process
        sample_id: Optional sample ID

    Returns:
        dict: Processing result with processed_file_id and metadata
    """
    try:
        # Update task state to STARTED
        self.update_state(
            state="STARTED", meta={"status": "Starting file processing", "progress": 0}
        )

        logger.info(f"Starting async processing for file {raw_file_id}")

        # Run the async processing
        result = run_async(
            _process_file_impl(raw_file_id=raw_file_id, sample_id=sample_id, task=self)
        )

        logger.info(f"Completed processing for file {raw_file_id}")
        return result

    except Exception as e:
        logger.error(f"Error processing file {raw_file_id}: {e}", exc_info=True)
        self.update_state(
            state="FAILURE", meta={"status": "Processing failed", "error": str(e)}
        )
        raise


async def _process_file_impl(raw_file_id: int, sample_id: str, task):
    """Internal async implementation of file processing"""

    async with AsyncSessionMaker() as db:
        # Get the raw file
        task.update_state(
            state="PROGRESS", meta={"status": "Loading file metadata", "progress": 10}
        )

        result = await db.execute(
            "SELECT * FROM datafiles WHERE id = :id", {"id": raw_file_id}
        )
        raw_file_row = result.first()

        if not raw_file_row:
            raise ValueError(f"File {raw_file_id} not found")

        # Convert row to dict
        raw_file_dict = dict(raw_file_row._mapping)

        task.update_state(
            state="PROGRESS", meta={"status": "Detecting file format", "progress": 20}
        )

        # Get file from storage and process
        from app.services.storage import ObjectStorage

        storage = ObjectStorage()

        # Download file content
        file_content = await storage.download_object(raw_file_dict["object_key"])

        task.update_state(
            state="PROGRESS",
            meta={"status": "Converting to unified format", "progress": 40},
        )

        # Process the file
        processor = FileProcessor()
        processing_result = await processor.process_file_async(
            file_content=file_content,
            filename=raw_file_dict["filename"],
            sample_id=sample_id,
        )

        task.update_state(
            state="PROGRESS", meta={"status": "Saving processed data", "progress": 70}
        )

        # Save processed file
        if processing_result["success"]:
            import json
            from datetime import datetime, timezone

            # Create processed file record
            processed_data_json = json.dumps(processing_result["unified_data"].dict())
            encrypted_data, checksum = await storage.encrypt_and_upload(
                processed_data_json.encode(),
                f"{raw_file_dict['filename']}.processed.json",
            )

            # Insert processed file
            insert_query = """
            INSERT INTO datafiles 
            (filename, object_key, size, checksum, project_id, uploaded_by_id, metadata_, created_at, updated_at)
            VALUES 
            (:filename, :object_key, :size, :checksum, :project_id, :uploaded_by_id, :metadata_, :created_at, :updated_at)
            RETURNING id
            """

            processed_metadata = {
                "omics_type": processing_result["omics_type"].value,
                "processing_info": processing_result.get("processing_info", {}),
            }

            result = await db.execute(
                insert_query,
                {
                    "filename": f"{raw_file_dict['filename']}.processed.json",
                    "object_key": encrypted_data["object_key"],
                    "size": len(processed_data_json),
                    "checksum": checksum,
                    "project_id": raw_file_dict["project_id"],
                    "uploaded_by_id": raw_file_dict["uploaded_by_id"],
                    "metadata_": json.dumps(processed_metadata),
                    "created_at": datetime.now(timezone.utc),
                    "updated_at": datetime.now(timezone.utc),
                },
            )

            processed_file_id = result.scalar()

            # Update raw file metadata with processed_file_id
            update_query = """
            UPDATE datafiles 
            SET metadata_ = jsonb_set(
                COALESCE(metadata_, '{}'::jsonb),
                '{processed_file_id}',
                to_jsonb(:processed_id::integer)
            ),
            metadata_ = jsonb_set(
                metadata_,
                '{omics_type}',
                to_jsonb(:omics_type::text)
            ),
            updated_at = :updated_at
            WHERE id = :raw_id
            """

            await db.execute(
                update_query,
                {
                    "processed_id": processed_file_id,
                    "omics_type": processing_result["omics_type"].value,
                    "updated_at": datetime.now(timezone.utc),
                    "raw_id": raw_file_id,
                },
            )

            await db.commit()

            task.update_state(
                state="PROGRESS",
                meta={"status": "Processing complete", "progress": 100},
            )

            return {
                "success": True,
                "processed_file_id": processed_file_id,
                "omics_type": processing_result["omics_type"].value,
                "record_count": len(processing_result["unified_data"].records),
            }
        else:
            task.update_state(
                state="FAILURE",
                meta={
                    "status": "Conversion failed",
                    "error": processing_result.get("error"),
                },
            )
            raise ValueError(f"Processing failed: {processing_result.get('error')}")


@celery_app.task(bind=True, name="app.tasks.file_tasks.process_large_file_async")
def process_large_file_async(self, raw_file_id: int, sample_id: str = None):
    """
    Process large files with streaming and progress updates
    Same as process_file_async but with optimizations for large files
    """
    # For now, use the same implementation
    # In production, this could use chunked processing
    return process_file_async(raw_file_id, sample_id)
