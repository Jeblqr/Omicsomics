import io
import mimetypes
import json

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile, status
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.services import datafiles as datafile_service
from app.services import projects as project_service
from app.services.file_processor import FileProcessor

router = APIRouter()


@router.post("/upload", status_code=status.HTTP_201_CREATED)
async def upload_datafile(
    project_id: int = Form(...),
    file: UploadFile = File(...),
    process_file: bool = Form(True),  # Enable/disable file processing
    async_processing: bool = Form(False),  # Use async Celery task for large files
    sample_id: str = Form(None),  # Optional sample ID
    organism: str = Form(None),  # Optional organism
    reference_genome: str = Form(None),  # Optional reference genome
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Upload a data file with optional processing to unified format.

    Args:
        project_id: Project to upload to
        file: File to upload
        process_file: Whether to process file to unified format (default: True)
        async_processing: Process file asynchronously using Celery (default: False)
        sample_id: Optional sample identifier
        organism: Optional organism name (e.g., "Homo sapiens")
        reference_genome: Optional reference genome (e.g., "hg38")
    """
    # Verify project exists and user is owner
    project = await project_service.get_project(db, project_id)
    if project is None:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    filename = file.filename or "unnamed"
    mime_type = (
        file.content_type
        or mimetypes.guess_type(filename)[0]
        or "application/octet-stream"
    )
    content = await file.read()
    raw_file_obj = io.BytesIO(content)

    # Store raw file first
    df = await datafile_service.create_datafile(
        db=db,
        file_data=raw_file_obj,
        filename=filename,
        project_id=project_id,
        user_id=current_user.id,
        content_type=mime_type,
    )

    # Process file if requested
    processing_info = {"processed": False}
    if process_file and async_processing:
        # Async processing using Celery
        try:
            from app.tasks.file_tasks import process_file_async

            # Submit task to Celery
            task = process_file_async.delay(raw_file_id=df.id, sample_id=sample_id)

            processing_info = {
                "processed": False,
                "async": True,
                "task_id": task.id,
                "status": "submitted",
                "message": "File processing submitted to background queue",
            }

        except Exception as e:
            import logging

            logger = logging.getLogger(__name__)
            logger.error(f"Failed to submit async task: {str(e)}")
            processing_info = {
                "processed": False,
                "async": False,
                "error": "Failed to submit processing task",
            }

    elif process_file:
        # Synchronous processing (original behavior)
        try:
            # Process file to unified format
            file_obj_for_processing = io.BytesIO(content)
            metadata = {}
            if organism:
                metadata["organism"] = organism
            if reference_genome:
                metadata["reference_genome"] = reference_genome

            unified_data, proc_info = await FileProcessor.process_file(
                file_data=file_obj_for_processing,
                filename=filename,
                project_id=project_id,
                sample_id=sample_id,
                metadata=metadata,
            )

            # Store processed version if conversion was successful
            if proc_info.get("converted", False):
                processed_filename = f"{filename}.processed.json"
                processed_content = FileProcessor.serialize_unified_data(unified_data)
                processed_file_obj = io.BytesIO(processed_content)

                # Store processed version
                processed_df = await datafile_service.create_datafile(
                    db=db,
                    file_data=processed_file_obj,
                    filename=processed_filename,
                    project_id=project_id,
                    user_id=current_user.id,
                    content_type="application/json",
                )

                # Update metadata to link raw and processed files
                from sqlalchemy.orm.attributes import flag_modified

                df.metadata_ = df.metadata_ or {}
                df.metadata_["processed_file_id"] = processed_df.id
                df.metadata_["processing_info"] = proc_info
                df.metadata_["omics_type"] = proc_info.get("omics_type")
                flag_modified(df, "metadata_")

                processed_df.metadata_ = processed_df.metadata_ or {}
                processed_df.metadata_["raw_file_id"] = df.id
                processed_df.metadata_["is_processed"] = True
                flag_modified(processed_df, "metadata_")

                await db.commit()
                await db.refresh(df)

                processing_info = {
                    "processed": True,
                    "processed_file_id": processed_df.id,
                    "processed_filename": processed_filename,
                    **proc_info,
                }
            else:
                # Store processing info even if conversion failed
                from sqlalchemy.orm.attributes import flag_modified

                df.metadata_ = df.metadata_ or {}
                df.metadata_["processing_info"] = proc_info
                df.metadata_["omics_type"] = proc_info.get("omics_type")
                flag_modified(df, "metadata_")
                await db.commit()
                await db.refresh(df)

                processing_info = {
                    "processed": False,
                    "reason": proc_info.get("reason")
                    or proc_info.get("error", "Unknown"),
                }
        except Exception as e:
            # Log error but don't fail the upload
            import logging

            logger = logging.getLogger(__name__)
            logger.error(f"File processing failed for {filename}: {str(e)}")
            processing_info = {
                "processed": False,
                "error": str(e),
            }

    return {
        "id": df.id,
        "filename": df.filename,
        "object_key": df.object_key,
        "metadata_": df.metadata_,
        "size": df.size,
        "checksum": df.checksum,
        "project_id": df.project_id,
        "run_id": df.run_id,
        "uploaded_by_id": df.uploaded_by_id,
        "created_at": df.created_at,
        "updated_at": df.updated_at,
        "processing": processing_info,
    }


@router.get("/")
async def list_datafiles(
    project_id: int | None = None,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    if project_id is not None:
        project = await project_service.get_project(db, project_id)
        if project is None:
            raise HTTPException(status_code=404, detail="Project not found")
        if project.owner_id != current_user.id:
            raise HTTPException(status_code=403, detail="Not authorized")

    rows = await datafile_service.list_datafiles(db, project_id)
    return [
        {
            "id": r.id,
            "filename": r.filename,
            "object_key": r.object_key,
            "metadata_": r.metadata_,
            "size": r.size,
            "checksum": r.checksum,
            "project_id": r.project_id,
            "run_id": r.run_id,
            "uploaded_by_id": r.uploaded_by_id,
            "created_at": r.created_at,
            "updated_at": r.updated_at,
        }
        for r in rows
    ]


@router.get("/{datafile_id}/download")
async def download_datafile(
    datafile_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Return presigned URL to encrypted object (client must decrypt or use download-decrypted endpoint)."""
    df = await datafile_service.get_datafile(db, datafile_id)
    if df is None:
        raise HTTPException(status_code=404, detail="DataFile not found")

    project = await project_service.get_project(db, df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    from app.storage.s3 import s3_client

    url = await s3_client.get_presigned_url(df.object_key)
    return {"download_url": url}


@router.get("/{datafile_id}/download-decrypted")
async def download_datafile_decrypted(
    datafile_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Server-side decrypt and stream the file content."""
    from fastapi.responses import StreamingResponse
    from app.services import storage_service

    df = await datafile_service.get_datafile(db, datafile_id)
    if df is None:
        raise HTTPException(status_code=404, detail="DataFile not found")

    project = await project_service.get_project(db, df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Download and decrypt
    plaintext = await storage_service.download_and_decrypt(
        current_user.id, df.object_key
    )

    # Stream back to client
    return StreamingResponse(
        io.BytesIO(plaintext),
        media_type="application/octet-stream",
        headers={"Content-Disposition": f'attachment; filename="{df.filename}"'},
    )


@router.get("/{datafile_id}/processed")
async def get_processed_data(
    datafile_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Get the processed (unified format) version of a data file.
    Returns the unified format JSON data.
    """
    from app.services import storage_service

    df = await datafile_service.get_datafile(db, datafile_id)
    if df is None:
        raise HTTPException(status_code=404, detail="DataFile not found")

    project = await project_service.get_project(db, df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Check if this file has a processed version
    metadata = df.metadata_ or {}
    processed_file_id = metadata.get("processed_file_id")

    if not processed_file_id:
        raise HTTPException(
            status_code=404, detail="No processed version available for this file"
        )

    # Get processed file
    processed_df = await datafile_service.get_datafile(db, processed_file_id)
    if not processed_df:
        raise HTTPException(status_code=404, detail="Processed file not found")

    # Download and decrypt processed file
    plaintext = await storage_service.download_and_decrypt(
        current_user.id, processed_df.object_key
    )

    # Parse and return unified data
    try:
        unified_data = FileProcessor.deserialize_unified_data(plaintext)
        return {
            "datafile_id": datafile_id,
            "processed_file_id": processed_file_id,
            "unified_data": unified_data.model_dump(mode="json"),
            "processing_info": metadata.get("processing_info", {}),
        }
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to parse processed data: {str(e)}"
        )


@router.delete("/{datafile_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_datafile(
    datafile_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Delete a data file and its processed version if it exists."""
    df = await datafile_service.get_datafile(db, datafile_id)
    if df is None:
        raise HTTPException(status_code=404, detail="DataFile not found")

    project = await project_service.get_project(db, df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Delete processed version if it exists
    metadata = df.metadata_ or {}
    processed_file_id = metadata.get("processed_file_id")
    if processed_file_id:
        await datafile_service.delete_datafile(db, processed_file_id)

    # Delete the main file
    await datafile_service.delete_datafile(db, datafile_id)
    return None


@router.get("/task/{task_id}/status")
async def get_task_status(
    task_id: str,
    current_user: User = Depends(get_current_user),
):
    """
    Get the status of an async processing task.

    Args:
        task_id: Celery task ID

    Returns:
        Task status including state, progress, and result
    """
    from app.celery_app import celery_app
    from celery.result import AsyncResult

    task = AsyncResult(task_id, app=celery_app)

    if task.state == "PENDING":
        response = {
            "state": task.state,
            "status": "Task is waiting to be processed",
            "progress": 0,
        }
    elif task.state == "STARTED":
        response = {
            "state": task.state,
            "status": task.info.get("status", "Processing..."),
            "progress": task.info.get("progress", 0),
        }
    elif task.state == "PROGRESS":
        response = {
            "state": task.state,
            "status": task.info.get("status", "Processing..."),
            "progress": task.info.get("progress", 0),
        }
    elif task.state == "SUCCESS":
        response = {
            "state": task.state,
            "status": "Processing complete",
            "progress": 100,
            "result": task.result,
        }
    elif task.state == "FAILURE":
        response = {
            "state": task.state,
            "status": "Processing failed",
            "progress": 0,
            "error": str(task.info),
        }
    else:
        response = {"state": task.state, "status": str(task.info)}

    return response
