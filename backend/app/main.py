from fastapi import FastAPI

from app.api.routers import api_router

app = FastAPI(title="Omicsomics API", version="0.1.0")

app.include_router(api_router)


@app.get("/healthz", tags=["health"])
def healthcheck() -> dict[str, str]:
  """Lightweight health endpoint for probes."""
  return {"status": "ok"}
