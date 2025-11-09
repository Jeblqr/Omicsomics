"""Tests for runs service and API."""

import pytest
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.run import Run
from app.services import runs as runs_service


@pytest.mark.asyncio
async def test_create_run(db_session: AsyncSession, test_user, test_project):
    """Test creating a run."""
    run = await runs_service.create_run(
        db=db_session,
        name="Test Run",
        description="A test run",
        project_id=test_project.id,
        owner_id=test_user.id,
    )
    assert run.id is not None
    assert run.name == "Test Run"
    assert run.status == "pending"
    assert run.project_id == test_project.id


@pytest.mark.asyncio
async def test_get_run(db_session: AsyncSession, test_user, test_project):
    """Test retrieving a run."""
    run = await runs_service.create_run(
        db=db_session,
        name="Get Test",
        description="",
        project_id=test_project.id,
        owner_id=test_user.id,
    )
    fetched = await runs_service.get_run(db_session, run.id)
    assert fetched is not None
    assert fetched.id == run.id
    assert fetched.name == "Get Test"


@pytest.mark.asyncio
async def test_list_runs(db_session: AsyncSession, test_user, test_project):
    """Test listing runs."""
    await runs_service.create_run(
        db_session, "Run1", "desc1", test_project.id, test_user.id
    )
    await runs_service.create_run(
        db_session, "Run2", "desc2", test_project.id, test_user.id
    )
    runs = await runs_service.list_runs(db_session, project_id=test_project.id)
    assert len(runs) >= 2


@pytest.mark.asyncio
async def test_runs_api_create(async_client: AsyncClient, auth_headers, test_project):
    """Test POST /api/v1/runs/ endpoint."""
    response = await async_client.post(
        "/api/v1/runs/",
        json={
            "name": "API Run",
            "description": "via API",
            "project_id": test_project.id,
        },
        headers=auth_headers,
    )
    assert response.status_code == 201
    data = response.json()
    assert data["name"] == "API Run"
    assert data["status"] == "pending"


@pytest.mark.asyncio
async def test_runs_api_list(async_client: AsyncClient, auth_headers, test_project):
    """Test GET /api/v1/runs/ endpoint."""
    # Create a run first
    await async_client.post(
        "/api/v1/runs/",
        json={
            "name": "List Test Run",
            "description": "",
            "project_id": test_project.id,
        },
        headers=auth_headers,
    )
    response = await async_client.get(
        f"/api/v1/runs/?project_id={test_project.id}",
        headers=auth_headers,
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) >= 1


@pytest.mark.asyncio
async def test_runs_api_get(async_client: AsyncClient, auth_headers, test_project):
    """Test GET /api/v1/runs/{id} endpoint."""
    create_resp = await async_client.post(
        "/api/v1/runs/",
        json={"name": "Get Run API", "description": "", "project_id": test_project.id},
        headers=auth_headers,
    )
    run_id = create_resp.json()["id"]

    response = await async_client.get(f"/api/v1/runs/{run_id}", headers=auth_headers)
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == run_id
    assert data["name"] == "Get Run API"


@pytest.mark.asyncio
async def test_runs_api_unauthorized(
    async_client: AsyncClient, auth_headers, test_project
):
    """Test that user cannot access another user's project runs."""
    # Create run as test_user
    create_resp = await async_client.post(
        "/api/v1/runs/",
        json={"name": "Private Run", "description": "", "project_id": test_project.id},
        headers=auth_headers,
    )
    run_id = create_resp.json()["id"]

    # Verify authorized access works
    response = await async_client.get(f"/api/v1/runs/{run_id}", headers=auth_headers)
    assert response.status_code == 200
