from functools import lru_cache

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
  model_config = SettingsConfigDict(env_file=".env", env_file_encoding="utf-8", extra="ignore")

  api_host: str = Field("0.0.0.0", env="API_HOST")
  api_port: int = Field(8001, env="API_PORT")
  database_url: str = Field(
    "postgresql+asyncpg://postgres:postgres@db:5432/omicsomics",
    env="DATABASE_URL",
  )
  db_echo: bool = Field(False, env="DB_ECHO")
  object_storage_endpoint: str = Field("http://localhost:9000", env="OBJECT_STORAGE_ENDPOINT")
  object_storage_access_key: str = Field("minio", env="OBJECT_STORAGE_ACCESS_KEY")
  object_storage_secret_key: str = Field("minio123", env="OBJECT_STORAGE_SECRET_KEY")
  object_storage_bucket: str = Field("omicsomics", env="OBJECT_STORAGE_BUCKET")

  # Security
  SECRET_KEY: str = Field(
      "5b1b472e91ed11fb37ac3452877cc5aa6beda975a6d0f61fe1d80e715ed1a83f", env="SECRET_KEY"
  )
  ACCESS_TOKEN_EXPIRE_MINUTES: int = Field(60 * 24 * 8, env="ACCESS_TOKEN_EXPIRE_MINUTES")  # 8 days


@lru_cache()
def get_settings() -> "Settings":
  return Settings()


settings = get_settings()
