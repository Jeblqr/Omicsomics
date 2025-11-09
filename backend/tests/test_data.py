"""Tests for data files service and API."""

import io
import pytest
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from app.services import datafiles as datafile_service


@pytest.mark.asyncio
async def test_create_datafile(db_session: AsyncSession, test_user, test_project):
    """Test creating a datafile with encryption."""
    content = b"This is test file content"
    file_obj = io.BytesIO(content)

    df = await datafile_service.create_datafile(
        db=db_session,
        file_data=file_obj,
        filename="test.txt",
        project_id=test_project.id,
        user_id=test_user.id,
        content_type="text/plain",
    )

    assert df.id is not None
    assert df.filename == "test.txt"
    assert df.size == len(content)
    assert df.checksum is not None
    assert df.object_key.startswith(f"sandbox/{test_project.id}/{test_user.id}/")


@pytest.mark.asyncio
async def test_list_datafiles(db_session: AsyncSession, test_user, test_project):
    """Test listing datafiles."""
    file_obj = io.BytesIO(b"content1")
    await datafile_service.create_datafile(
        db_session, file_obj, "file1.txt", test_project.id, test_user.id
    )

    file_obj2 = io.BytesIO(b"content2")
    await datafile_service.create_datafile(
        db_session, file_obj2, "file2.txt", test_project.id, test_user.id
    )

    files = await datafile_service.list_datafiles(
        db_session, project_id=test_project.id
    )
    assert len(files) >= 2


@pytest.mark.asyncio
async def test_data_api_upload(async_client: AsyncClient, auth_headers, test_project):
    """Test POST /api/v1/data/upload endpoint."""
    files = {"file": ("testfile.txt", b"Hello World", "text/plain")}
    data = {"project_id": test_project.id}

    response = await async_client.post(
        "/api/v1/data/upload",
        data=data,
        files=files,
        headers=auth_headers,
    )

    assert response.status_code == 201
    json_data = response.json()
    assert json_data["filename"] == "testfile.txt"
    assert "object_key" in json_data


@pytest.mark.asyncio
async def test_data_api_list(async_client: AsyncClient, auth_headers, test_project):
    """Test GET /api/v1/data/ endpoint."""
    # Upload a file first
    files = {"file": ("list_test.txt", b"content", "text/plain")}
    data = {"project_id": test_project.id}
    await async_client.post(
        "/api/v1/data/upload", data=data, files=files, headers=auth_headers
    )

    response = await async_client.get(
        f"/api/v1/data/?project_id={test_project.id}",
        headers=auth_headers,
    )

    assert response.status_code == 200
    json_data = response.json()
    assert isinstance(json_data, list)
    assert len(json_data) >= 1


@pytest.mark.asyncio
async def test_data_api_download_decrypted(
    async_client: AsyncClient, auth_headers, test_project
):
    """Test GET /api/v1/data/{id}/download-decrypted endpoint."""
    original_content = b"Secret encrypted content"
    files = {"file": ("encrypted.bin", original_content, "application/octet-stream")}
    data = {"project_id": test_project.id}

    upload_resp = await async_client.post(
        "/api/v1/data/upload", data=data, files=files, headers=auth_headers
    )
    datafile_id = upload_resp.json()["id"]

    # Download decrypted
    response = await async_client.get(
        f"/api/v1/data/{datafile_id}/download-decrypted",
        headers=auth_headers,
    )

    assert response.status_code == 200
    assert response.content == original_content
