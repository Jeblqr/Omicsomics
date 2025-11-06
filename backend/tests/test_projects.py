import pytest


@pytest.mark.asyncio
async def test_create_and_list_projects(async_client) -> None:
  create_response = await async_client.post(
    "/api/v1/projects",
    json={"name": "Test Project", "description": "Integration test"},
  )
  assert create_response.status_code == 201
  created = create_response.json()
  assert created["name"] == "Test Project"
  assert created["description"] == "Integration test"

  list_response = await async_client.get("/api/v1/projects")
  assert list_response.status_code == 200
  items = list_response.json()
  assert any(item["name"] == "Test Project" for item in items)
