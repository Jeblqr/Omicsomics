from __future__ import annotations

import asyncio
from collections.abc import AsyncGenerator, Generator
from typing import AsyncIterator

import pytest
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine

from app.api.dependencies.database import get_db
from app.main import app
from app.models.base import Base

TEST_DATABASE_URL = "sqlite+aiosqlite:///:memory:"


@pytest.fixture(scope="session")
def event_loop() -> Generator[asyncio.AbstractEventLoop, None, None]:
  loop = asyncio.new_event_loop()
  yield loop
  loop.close()


@pytest.fixture()
async def db_session() -> AsyncIterator[async_sessionmaker[AsyncSession]]:
  engine = create_async_engine(TEST_DATABASE_URL, future=True)
  async with engine.begin() as conn:
    await conn.run_sync(Base.metadata.create_all)

  session_maker = async_sessionmaker(engine, expire_on_commit=False)

  async def override_get_db() -> AsyncGenerator[AsyncSession, None]:
    async with session_maker() as session:
      yield session

  app.dependency_overrides[get_db] = override_get_db

  try:
    yield session_maker
  finally:
    app.dependency_overrides.pop(get_db, None)
    await engine.dispose()


@pytest.fixture()
async def async_client(
  db_session: async_sessionmaker[AsyncSession],
) -> AsyncIterator[AsyncClient]:
  async with AsyncClient(app=app, base_url="http://testserver") as client:
    yield client
