from collections.abc import AsyncGenerator, Callable

from sqlalchemy.ext.asyncio import AsyncSession, create_async_engine
from sqlalchemy.orm import sessionmaker

from app.settings import settings

_engine = create_async_engine(settings.database_url, echo=False, future=True)
_async_session_factory: Callable[[], AsyncSession] = sessionmaker(
  bind=_engine,
  class_=AsyncSession,
  expire_on_commit=False,
)


def async_session() -> AsyncGenerator[AsyncSession, None]:
  async def _session() -> AsyncGenerator[AsyncSession, None]:
    async with _async_session_factory() as session:
      yield session

  return _session()
