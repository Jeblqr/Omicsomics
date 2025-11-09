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


@router.get("/{datafile_id}/preview")
async def preview_archive(
    datafile_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Preview contents of an archive file (ZIP, TAR, TAR.GZ, TAR.BZ2).

    Returns list of files in the archive with metadata.

    Returns:
        {
            "is_archive": bool,
            "files": [{"name": str, "path": str, "size": int, "mime_type": str}],
            "total_files": int,
            "total_size": int
        }
    """
    from app.services import storage_service
    from app.utils.archive import list_archive_contents, is_archive_file, ArchiveError
    import tempfile
    from pathlib import Path

    # Get datafile
    df = await datafile_service.get_datafile(db, datafile_id)
    if df is None:
        raise HTTPException(status_code=404, detail="DataFile not found")

    # Check authorization
    project = await project_service.get_project(db, df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Check if file is an archive
    file_path = Path(df.filename)
    if not is_archive_file(file_path):
        return {
            "is_archive": False,
            "message": f"File '{df.filename}' is not a supported archive format",
        }

    # Download and decrypt file
    try:
        plaintext = await storage_service.download_and_decrypt(
            current_user.id, df.object_key
        )
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to download file: {str(e)}"
        )

    # Write to temporary file for archive processing
    with tempfile.NamedTemporaryFile(delete=False, suffix=file_path.suffix) as tmp_file:
        tmp_path = Path(tmp_file.name)
        tmp_file.write(plaintext)

    try:
        # List archive contents
        files = list_archive_contents(tmp_path)

        # Calculate totals
        total_size = sum(f.size for f in files)

        return {
            "is_archive": True,
            "files": [f.to_dict() for f in files],
            "total_files": len(files),
            "total_size": total_size,
            "archive_filename": df.filename,
        }

    except ArchiveError as e:
        raise HTTPException(status_code=400, detail=f"Failed to read archive: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Unexpected error: {str(e)}")
    finally:
        # Clean up temporary file
        if tmp_path.exists():
            tmp_path.unlink()


@router.post("/{datafile_id}/process-from-archive", status_code=status.HTTP_201_CREATED)
async def process_file_from_archive(
    datafile_id: int,
    file_path: str = Form(...),
    sample_id: str = Form(None),
    organism: str = Form(None),
    reference_genome: str = Form(None),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Extract and process a specific file from an archive.

    Args:
        datafile_id: ID of the archive file
        file_path: Path of the file within the archive to extract and process
        sample_id: Optional sample identifier
        organism: Optional organism name
        reference_genome: Optional reference genome

    Returns:
        {
            "datafile_id": int,  # New datafile created from extracted file
            "processed_file_id": int,  # Processed unified format file
            "extracted_filename": str,
            "processing_info": dict
        }
    """
    from app.services import storage_service
    from app.utils.archive import extract_archive, is_archive_file, ArchiveError
    import tempfile
    from pathlib import Path

    # Get archive datafile
    archive_df = await datafile_service.get_datafile(db, datafile_id)
    if archive_df is None:
        raise HTTPException(status_code=404, detail="Archive file not found")

    # Check authorization
    project = await project_service.get_project(db, archive_df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Verify it's an archive
    archive_path = Path(archive_df.filename)
    if not is_archive_file(archive_path):
        raise HTTPException(
            status_code=400,
            detail=f"File '{archive_df.filename}' is not a supported archive format",
        )

    # Download and decrypt archive
    try:
        plaintext = await storage_service.download_and_decrypt(
            current_user.id, archive_df.object_key
        )
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to download archive: {str(e)}"
        )

    # Create temporary directory for extraction
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Write archive to temp file
        archive_tmp = tmp_path / archive_df.filename
        archive_tmp.write_bytes(plaintext)

        # Extract specific file
        extract_dir = tmp_path / "extracted"
        try:
            extracted_file = extract_archive(archive_tmp, extract_dir, file_path)
        except ArchiveError as e:
            raise HTTPException(
                status_code=400, detail=f"Failed to extract file: {str(e)}"
            )

        if not extracted_file.exists():
            raise HTTPException(
                status_code=404, detail=f"File '{file_path}' not found in archive"
            )

        # Read extracted file
        extracted_content = extracted_file.read_bytes()
        extracted_filename = extracted_file.name

        # Guess MIME type
        mime_type = (
            mimetypes.guess_type(extracted_filename)[0] or "application/octet-stream"
        )

    # Create new datafile for extracted content
    extracted_df = await datafile_service.create_datafile(
        db=db,
        file_data=io.BytesIO(extracted_content),
        filename=extracted_filename,
        project_id=archive_df.project_id,
        user_id=current_user.id,
        content_type=mime_type,
    )

    # Update metadata to link to parent archive
    extracted_metadata = extracted_df.metadata_ or {}
    extracted_metadata["extracted_from_archive"] = {
        "archive_id": datafile_id,
        "archive_filename": archive_df.filename,
        "file_path": file_path,
    }
    extracted_df.metadata_ = extracted_metadata
    await db.commit()

    # Process the extracted file
    try:
        # Process file
        result = await FileProcessor.process_file_async(
            file_content=extracted_content,
            filename=extracted_filename,
            sample_id=sample_id,
        )

        if not result.get("success"):
            raise Exception(result.get("error", "Unknown processing error"))

        unified_data = result["unified_data"]

        # Serialize unified data
        unified_json = FileProcessor.serialize_unified_data(unified_data)

        # Store processed file
        processed_df = await datafile_service.create_datafile(
            db=db,
            file_data=io.BytesIO(unified_json),
            filename=f"processed_{extracted_filename}.json",
            project_id=archive_df.project_id,
            user_id=current_user.id,
            content_type="application/json",
        )

        # Link files
        metadata = extracted_df.metadata_ or {}
        metadata["processed_file_id"] = processed_df.id
        metadata["processing_info"] = {
            "success": True,
            "format_detected": unified_data.format,
            "data_type": unified_data.data_type,
            "features_count": len(unified_data.features),
        }
        extracted_df.metadata_ = metadata
        await db.commit()

        return {
            "datafile_id": extracted_df.id,
            "processed_file_id": processed_df.id,
            "extracted_filename": extracted_filename,
            "archive_source": {
                "archive_id": datafile_id,
                "archive_filename": archive_df.filename,
                "file_path": file_path,
            },
            "processing_info": metadata["processing_info"],
        }

    except Exception as e:
        # Update metadata with error
        metadata = extracted_df.metadata_ or {}
        metadata["processing_info"] = {
            "success": False,
            "error": str(e),
        }
        extracted_df.metadata_ = metadata
        await db.commit()

        raise HTTPException(
            status_code=500, detail=f"Failed to process extracted file: {str(e)}"
        )


@router.get("/{datafile_id}/export")
async def export_processed_data(
    datafile_id: int,
    format: str = "csv",  # Export format: csv, tsv, excel, json
    columns: (
        str | None
    ) = None,  # Comma-separated list of columns to include (None = all)
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Export processed (unified format) data to various formats.

    Args:
        datafile_id: ID of the datafile to export
        format: Export format (csv, tsv, excel, json)
        columns: Comma-separated list of columns to include (optional)

    Returns:
        Exported data in specified format
    """
    from fastapi.responses import StreamingResponse
    from app.services import storage_service
    from app.utils.export import DataExporter

    # Get datafile
    df = await datafile_service.get_datafile(db, datafile_id)
    if df is None:
        raise HTTPException(status_code=404, detail="DataFile not found")

    # Check authorization
    project = await project_service.get_project(db, df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Get processed file ID
    metadata = df.metadata_ or {}
    processed_file_id = metadata.get("processed_file_id")

    if not processed_file_id:
        raise HTTPException(
            status_code=404,
            detail="No processed data available. Please process the file first.",
        )

    # Get processed file
    processed_df = await datafile_service.get_datafile(db, processed_file_id)
    if processed_df is None:
        raise HTTPException(status_code=404, detail="Processed file not found")

    # Download and decrypt processed file
    try:
        plaintext = await storage_service.download_and_decrypt(
            current_user.id, processed_df.object_key
        )
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to retrieve processed data: {str(e)}"
        )

    # Deserialize unified data
    try:
        unified_data = FileProcessor.deserialize_unified_data(plaintext)
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to parse processed data: {str(e)}"
        )

    # Parse columns parameter
    include_columns = None
    if columns:
        include_columns = [col.strip() for col in columns.split(",")]

    # Export data
    try:
        exported_data = DataExporter.export(
            unified_data, format=format, include_columns=include_columns
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except ImportError as e:
        raise HTTPException(
            status_code=500, detail=f"Export format not available: {str(e)}"
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to export data: {str(e)}")

    # Get MIME type and file extension
    mime_type = DataExporter.get_mime_type(format)
    file_extension = DataExporter.get_file_extension(format)

    # Generate filename
    base_filename = df.filename.rsplit(".", 1)[0] if "." in df.filename else df.filename
    export_filename = f"{base_filename}_exported.{file_extension}"

    # Return as streaming response
    return StreamingResponse(
        io.BytesIO(exported_data),
        media_type=mime_type,
        headers={"Content-Disposition": f'attachment; filename="{export_filename}"'},
    )


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


@router.post("/search")
async def search_data(
    query: str | None = None,
    feature_id: str | None = None,
    field_ranges: str | None = None,  # JSON string: {"field": {"min": 0, "max": 100}}
    exact_matches: str | None = None,  # JSON string: {"field": "value"}
    datafile_id: int | None = None,
    project_id: int | None = None,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Advanced search across processed data files.

    Args:
        query: Full-text search query
        feature_id: Search by feature/gene/protein ID
        field_ranges: JSON string of range queries, e.g. {"expression": {"min": 10, "max": 100}}
        exact_matches: JSON string of exact match filters, e.g. {"chromosome": "chr1"}
        datafile_id: Optional: search specific datafile
        project_id: Optional: search within project

    Returns:
        {
            "results": [...],
            "total_results": int,
            "datafiles_searched": int
        }
    """
    from app.services import storage_service
    from app.services.search import DataSearchService, SearchQuery

    # Build search query
    search_query = SearchQuery()

    if query:
        search_query.with_text(query)

    if feature_id:
        search_query.with_feature_id(feature_id)

    if field_ranges:
        try:
            ranges = json.loads(field_ranges)
            for field, range_spec in ranges.items():
                search_query.with_range(
                    field, min_val=range_spec.get("min"), max_val=range_spec.get("max")
                )
        except json.JSONDecodeError:
            raise HTTPException(status_code=400, detail="Invalid field_ranges JSON")

    if exact_matches:
        try:
            matches = json.loads(exact_matches)
            for field, value in matches.items():
                search_query.with_exact_match(field, value)
        except json.JSONDecodeError:
            raise HTTPException(status_code=400, detail="Invalid exact_matches JSON")

    # Get datafiles to search
    datafiles = []

    if datafile_id:
        # Search single datafile
        df = await datafile_service.get_datafile(db, datafile_id)
        if df is None:
            raise HTTPException(status_code=404, detail="DataFile not found")

        project = await project_service.get_project(db, df.project_id)
        if project is None or project.owner_id != current_user.id:
            raise HTTPException(status_code=403, detail="Not authorized")

        datafiles = [df]

    elif project_id:
        # Search all datafiles in project
        project = await project_service.get_project(db, project_id)
        if project is None:
            raise HTTPException(status_code=404, detail="Project not found")
        if project.owner_id != current_user.id:
            raise HTTPException(status_code=403, detail="Not authorized")

        datafiles = await datafile_service.list_datafiles(db, project_id)

    else:
        # Search all user's datafiles
        all_datafiles = await datafile_service.list_datafiles(db, None)
        # Filter to user's files
        user_projects = await project_service.get_projects(db, current_user.id)
        user_project_ids = {p.id for p in user_projects}
        datafiles = [df for df in all_datafiles if df.project_id in user_project_ids]

    # Search each datafile
    all_results = []
    datafiles_searched = 0

    for df in datafiles:
        # Check if has processed file
        metadata = df.metadata_ or {}
        processed_file_id = metadata.get("processed_file_id")

        if not processed_file_id:
            continue

        # Get processed file
        processed_df = await datafile_service.get_datafile(db, processed_file_id)
        if processed_df is None:
            continue

        try:
            # Download and decrypt
            plaintext = await storage_service.download_and_decrypt(
                current_user.id, processed_df.object_key
            )

            # Deserialize
            unified_data = FileProcessor.deserialize_unified_data(plaintext)

            # Search
            results = DataSearchService.search(unified_data, search_query)

            # Add source info to results
            for result in results:
                result["_source"] = {
                    "datafile_id": df.id,
                    "filename": df.filename,
                    "project_id": df.project_id,
                }

            all_results.extend(results)
            datafiles_searched += 1

        except Exception as e:
            import logging

            logger = logging.getLogger(__name__)
            logger.warning(f"Failed to search datafile {df.id}: {e}")
            continue

    return {
        "results": all_results,
        "total_results": len(all_results),
        "datafiles_searched": datafiles_searched,
    }


@router.post("/batch-upload", status_code=status.HTTP_201_CREATED)
async def batch_upload_datafiles(
    project_id: int = Form(...),
    files: list[UploadFile] = File(...),
    process_files: bool = Form(True),
    async_processing: bool = Form(False),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Batch upload multiple data files.

    Args:
        project_id: Project to upload to
        files: List of files to upload
        process_files: Whether to process files (default: True)
        async_processing: Use async processing (default: False)

    Returns:
        {
            "uploaded": [...],
            "failed": [...],
            "total": int,
            "success_count": int,
            "failure_count": int
        }
    """
    # Verify project
    project = await project_service.get_project(db, project_id)
    if project is None:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    uploaded = []
    failed = []

    for file in files:
        try:
            filename = file.filename or "unnamed"
            mime_type = (
                file.content_type
                or mimetypes.guess_type(filename)[0]
                or "application/octet-stream"
            )
            content = await file.read()

            # Store raw file
            df = await datafile_service.create_datafile(
                db=db,
                file_data=io.BytesIO(content),
                filename=filename,
                project_id=project_id,
                user_id=current_user.id,
                content_type=mime_type,
            )

            result = {
                "filename": filename,
                "datafile_id": df.id,
                "size": len(content),
                "processed": False,
            }

            # Process if requested
            if process_files:
                if async_processing:
                    try:
                        from app.tasks.file_tasks import process_file_async

                        task = process_file_async.delay(raw_file_id=df.id)
                        result["task_id"] = task.id
                        result["async"] = True
                    except Exception as e:
                        result["processing_error"] = str(e)
                else:
                    try:
                        file_obj = io.BytesIO(content)
                        processor = FileProcessor()
                        unified_data, processing_info = await processor.process_file(
                            file_data=file_obj, filename=filename, project_id=project_id
                        )

                        # Store processed file
                        unified_json = processor.serialize_unified_data(unified_data)
                        processed_df = await datafile_service.create_datafile(
                            db=db,
                            file_data=io.BytesIO(unified_json),
                            filename=f"processed_{filename}.json",
                            project_id=project_id,
                            user_id=current_user.id,
                            content_type="application/json",
                        )

                        # Link files
                        metadata = df.metadata_ or {}
                        metadata["processed_file_id"] = processed_df.id
                        metadata["processing_info"] = processing_info
                        df.metadata_ = metadata
                        await db.commit()

                        result["processed"] = True
                        result["processed_file_id"] = processed_df.id
                    except Exception as e:
                        result["processing_error"] = str(e)

            uploaded.append(result)

        except Exception as e:
            failed.append({"filename": file.filename or "unknown", "error": str(e)})

    return {
        "uploaded": uploaded,
        "failed": failed,
        "total": len(files),
        "success_count": len(uploaded),
        "failure_count": len(failed),
    }


@router.delete("/batch-delete")
async def batch_delete_datafiles(
    datafile_ids: str = Form(...),  # Comma-separated IDs
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Batch delete multiple data files.

    Args:
        datafile_ids: Comma-separated list of datafile IDs

    Returns:
        {
            "deleted": [...],
            "failed": [...],
            "total": int,
            "success_count": int,
            "failure_count": int
        }
    """
    # Parse IDs
    try:
        ids = [int(id.strip()) for id in datafile_ids.split(",")]
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid datafile IDs format")

    deleted = []
    failed = []

    for datafile_id in ids:
        try:
            # Get datafile
            df = await datafile_service.get_datafile(db, datafile_id)
            if df is None:
                failed.append(
                    {"datafile_id": datafile_id, "error": "DataFile not found"}
                )
                continue

            # Check authorization
            project = await project_service.get_project(db, df.project_id)
            if project is None or project.owner_id != current_user.id:
                failed.append({"datafile_id": datafile_id, "error": "Not authorized"})
                continue

            # Delete processed version if exists
            metadata = df.metadata_ or {}
            processed_file_id = metadata.get("processed_file_id")
            if processed_file_id:
                await datafile_service.delete_datafile(db, processed_file_id)

            # Delete main file
            await datafile_service.delete_datafile(db, datafile_id)

            deleted.append({"datafile_id": datafile_id, "filename": df.filename})

        except Exception as e:
            failed.append({"datafile_id": datafile_id, "error": str(e)})

    return {
        "deleted": deleted,
        "failed": failed,
        "total": len(ids),
        "success_count": len(deleted),
        "failure_count": len(failed),
    }


@router.get("/{datafile_id}/quality-report")
async def get_quality_report(
    datafile_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Generate comprehensive quality control report for a data file.

    Returns:
        {
            "record_count": int,
            "field_count": int,
            "missing_values": {...},
            "outliers": {...},
            "distributions": {...},
            "completeness": {...},
            "overall_quality_score": {...}
        }
    """
    from app.services.quality_control import QualityControlService
    from app.services import storage_service

    # Get datafile
    df = await datafile_service.get_datafile(db, datafile_id)
    if df is None:
        raise HTTPException(status_code=404, detail="DataFile not found")

    # Check authorization
    project = await project_service.get_project(db, df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Get file content from storage
    try:
        file_bytes = await storage_service.download_encrypted_bytes(df.object_key)
        if file_bytes is None:
            raise HTTPException(
                status_code=404, detail="File content not found in storage"
            )

        # Load processed data if available
        metadata = df.metadata_ or {}
        processed_file_id = metadata.get("processed_file_id")

        if processed_file_id:
            # Use processed file for QC
            processed_df = await datafile_service.get_datafile(db, processed_file_id)
            if processed_df:
                file_bytes = await storage_service.download_encrypted_bytes(
                    processed_df.object_key
                )

        # Parse unified format JSON
        try:
            unified_data_dict = json.loads(file_bytes.decode("utf-8"))
            from app.schemas.unified_format import UnifiedData

            unified_data = UnifiedData(**unified_data_dict)
        except Exception as parse_error:
            raise HTTPException(
                status_code=400,
                detail=f"File is not in valid unified format: {str(parse_error)}",
            )

        # Generate QC report
        qc_service = QualityControlService()
        report = qc_service.generate_report(unified_data)

        # Add file metadata
        report["datafile"] = {
            "id": df.id,
            "filename": df.filename,
            "size": df.size,
            "checksum": df.checksum,
        }

        return report

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to generate quality report: {str(e)}"
        )
