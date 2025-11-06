from __future__ import annotations

import asyncio
import os
from collections.abc import AsyncGenerator, Generator
from typing import AsyncIterator

import pytest
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import NullPool

from app.api.dependencies.database import get_async_db
from app.main import app
from app.models.base import Base
from app.core.security import create_access_token, get_password_hash

# Use PostgreSQL for testing in dev container, SQLite otherwise
TEST_DATABASE_URL = os.getenv(
    "TEST_DATABASE_URL",
    "postgresql+asyncpg://postgres:postgres@db:5432/omicsomics_test"
)


@pytest.fixture(scope="session")
def event_loop() -> Generator[asyncio.AbstractEventLoop, None, None]:
  loop = asyncio.new_event_loop()
  yield loop
  loop.close()


@pytest.fixture(autouse=True)
async def setup_database() -> AsyncIterator[None]:
  """Create and drop database tables for each test."""
  engine = create_async_engine(TEST_DATABASE_URL, poolclass=NullPool, echo=False)
  
  async with engine.begin() as conn:
    await conn.run_sync(Base.metadata.drop_all)
    await conn.run_sync(Base.metadata.create_all)
  
  yield
  
  async with engine.begin() as conn:
    await conn.run_sync(Base.metadata.drop_all)
  
  await engine.dispose()


@pytest.fixture()
async def db_session() -> AsyncIterator[AsyncSession]:
  """Provide a test database session."""
  engine = create_async_engine(TEST_DATABASE_URL, poolclass=NullPool, echo=False)
  session_maker = async_sessionmaker(engine, expire_on_commit=False)

  async def override_get_async_db() -> AsyncGenerator[AsyncSession, None]:
    async with session_maker() as session:
      yield session

  app.dependency_overrides[get_async_db] = override_get_async_db

  try:
    async with session_maker() as session:
      yield session
  finally:
    app.dependency_overrides.pop(get_async_db, None)
    await engine.dispose()


@pytest.fixture()
async def async_client(
  db_session: AsyncSession,
) -> AsyncIterator[AsyncClient]:
  """Provide an async HTTP client for testing."""
  async with AsyncClient(app=app, base_url="http://testserver") as client:
    yield client


@pytest.fixture()
async def test_user(db_session: AsyncSession):
  """Create a test user."""
  from app.models.user import User
  
  user = User(
    email="test@example.com",
    hashed_password=get_password_hash("testpassword123"),
    full_name="Test User",
  )
  db_session.add(user)
  await db_session.commit()
  await db_session.refresh(user)
  return user


@pytest.fixture()
def test_user_token(test_user):
  """Create authentication token for test user."""
  return create_access_token(subject=str(test_user.id))


@pytest.fixture()
def auth_headers(test_user_token):
  """Create authorization headers with test token."""
  return {"Authorization": f"Bearer {test_user_token}"}


@pytest.fixture()
async def test_project(db_session: AsyncSession, test_user):
  """Create a test project."""
  from app.models.project import Project
  
  project = Project(
    name="Test Project",
    description="A test project for unit tests",
    owner_id=test_user.id,
  )
  db_session.add(project)
  await db_session.commit()
  await db_session.refresh(project)
  return project


@pytest.fixture()
async def test_sample(db_session: AsyncSession, test_project):
  """Create a test sample."""
  from app.models.sample import Sample
  
  sample = Sample(
    name="Test Sample",
    sample_type="RNA-seq",
    project_id=test_project.id,
    metadata_={"organism": "Homo sapiens", "tissue": "brain"},
  )
  db_session.add(sample)
  await db_session.commit()
  await db_session.refresh(sample)
  return sample
