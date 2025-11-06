from __future__ import annotations

import asyncio
from logging.config import fileConfig

from alembic import context
from sqlalchemy import pool
from sqlalchemy.engine import Connection
from sqlalchemy.ext.asyncio import AsyncEngine

from app.database import engine
from app.models.base import Base
from app.settings import settings

# Ensure models are imported so Alembic can discover them
from app import models  # noqa: F401

config = context.config

if config.config_file_name is not None:
  fileConfig(config.config_file_name)

config.set_main_option("sqlalchemy.url", settings.database_url)

target_metadata = Base.metadata


def run_migrations_offline() -> None:
  url = config.get_main_option("sqlalchemy.url")
  context.configure(
    url=url,
    target_metadata=target_metadata,
    literal_binds=True,
    dialect_opts={"paramstyle": "named"},
  )

  with context.begin_transaction():
    context.run_migrations()


def do_run_migrations(connection: Connection) -> None:
  context.configure(connection=connection, target_metadata=target_metadata)

  with context.begin_transaction():
    context.run_migrations()


def run_migrations_online() -> None:
  connectable = AsyncEngine(engine)

  async def process_migrations() -> None:
    async with connectable.connect() as connection:
      await connection.run_sync(do_run_migrations)

  asyncio.run(process_migrations())


if context.is_offline_mode():
  run_migrations_offline()
else:
  run_migrations_online()
