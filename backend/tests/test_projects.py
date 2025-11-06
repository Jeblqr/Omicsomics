"""Test project endpoints."""
from __future__ import annotations

import pytest
from httpx import AsyncClient


@pytest.mark.asyncio
async def test_create_project_success(async_client: AsyncClient, auth_headers):
    """Test creating a project."""
    response = await async_client.post(
        "/api/v1/projects/",
        headers=auth_headers,
        json={
            "name": "New Research Project",
            "description": "A comprehensive omics analysis project",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == "New Research Project"
    assert data["description"] == "A comprehensive omics analysis project"
    assert "id" in data
    assert "owner_id" in data
    assert "created_at" in data


@pytest.mark.asyncio
async def test_create_project_unauthorized(async_client: AsyncClient):
    """Test creating project without auth fails."""
    response = await async_client.post(
        "/api/v1/projects/",
        json={
            "name": "Unauthorized Project",
            "description": "This should fail",
        },
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_list_projects(async_client: AsyncClient, auth_headers, test_project):
    """Test listing user's projects."""
    response = await async_client.get("/api/v1/projects/", headers=auth_headers)
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) >= 1
    assert any(p["id"] == test_project.id for p in data)


@pytest.mark.asyncio
async def test_get_project_by_id(async_client: AsyncClient, auth_headers, test_project):
    """Test getting a specific project."""
    response = await async_client.get(
        f"/api/v1/projects/{test_project.id}",
        headers=auth_headers,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == test_project.id
    assert data["name"] == test_project.name


@pytest.mark.asyncio
async def test_get_nonexistent_project(async_client: AsyncClient, auth_headers):
    """Test getting nonexistent project returns 404."""
    response = await async_client.get("/api/v1/projects/99999", headers=auth_headers)
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_update_project(async_client: AsyncClient, auth_headers, test_project):
    """Test updating a project."""
    response = await async_client.put(
        f"/api/v1/projects/{test_project.id}",
        headers=auth_headers,
        json={
            "name": "Updated Project Name",
            "description": "Updated description",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == "Updated Project Name"
    assert data["description"] == "Updated description"


@pytest.mark.asyncio
async def test_delete_project(async_client: AsyncClient, auth_headers, test_project):
    """Test deleting a project."""
    # Delete the project
    response = await async_client.delete(
        f"/api/v1/projects/{test_project.id}",
        headers=auth_headers,
    )
    assert response.status_code == 200
    
    # Verify it's deleted
    response = await async_client.get(
        f"/api/v1/projects/{test_project.id}",
        headers=auth_headers,
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_project_isolation(async_client: AsyncClient, db_session, test_user):
    """Test that users can only see their own projects."""
    from app.models.user import User
    from app.models.project import Project
    from app.core.security import get_password_hash, create_access_token
    
    # Create another user
    other_user = User(
        email="other@example.com",
        hashed_password=get_password_hash("otherpassword"),
        full_name="Other User",
    )
    db_session.add(other_user)
    await db_session.commit()
    await db_session.refresh(other_user)
    
    # Create project for other user
    other_project = Project(
        name="Other's Project",
        description="Should not be visible",
        owner_id=other_user.id,
    )
    db_session.add(other_project)
    await db_session.commit()
    
    # Try to access other user's project with test_user's token
    test_user_token = create_access_token(subject=str(test_user.id))
    test_user_headers = {"Authorization": f"Bearer {test_user_token}"}
    
    response = await async_client.get(
        f"/api/v1/projects/{other_project.id}",
        headers=test_user_headers,
    )
    assert response.status_code == 404

