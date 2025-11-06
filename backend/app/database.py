from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine

from app.settings import settings

engine = create_async_engine(settings.database_url, echo=settings.db_echo, future=True)
AsyncSessionMaker = async_sessionmaker(engine, expire_on_commit=False)

__all__ = ["engine", "AsyncSessionMaker", "AsyncSession"]
