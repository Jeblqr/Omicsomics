from collections.abc import AsyncGenerator

from sqlalchemy.ext.asyncio import AsyncSession

from app.database import AsyncSessionMaker


async def get_db() -> AsyncGenerator[AsyncSession, None]:
  async with AsyncSessionMaker() as session:
    yield session
