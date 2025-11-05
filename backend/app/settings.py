from functools import lru_cache

from pydantic import BaseSettings, Field


class Settings(BaseSettings):
  database_url: str = Field(..., env="DATABASE_URL")
  object_storage_endpoint: str = Field(..., env="OBJECT_STORAGE_ENDPOINT")
  object_storage_access_key: str = Field(..., env="OBJECT_STORAGE_ACCESS_KEY")
  object_storage_secret_key: str = Field(..., env="OBJECT_STORAGE_SECRET_KEY")
  object_storage_bucket: str = Field(..., env="OBJECT_STORAGE_BUCKET")

  class Config:
    env_file = ".env"
    env_file_encoding = "utf-8"


@lru_cache()
def get_settings() -> "Settings":
  return Settings()


settings = get_settings()
