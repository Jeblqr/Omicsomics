import os
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.api.routers import api_router

app = FastAPI(
    title="Omicsomics API",
    version="0.1.0",
    # Increase file upload size limit to 5GB
    # Note: Also need to configure nginx/reverse proxy if using one
)

# Define allowed origins for local development
allowed_origins = [
    "http://localhost:5173",
    "http://localhost:5174",
    "http://127.0.0.1:5173",
    "http://127.0.0.1:5174",
]

# Regex to allow all GitHub Codespaces frontend URLs
allow_origin_regex = r"https://.*-5174\.app\.github\.dev"

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_origin_regex=allow_origin_regex,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(api_router, prefix="/api/v1")
app.include_router(api_router, prefix="/api")


@app.get("/healthz", tags=["health"])
def healthcheck() -> dict[str, str]:
    """Lightweight health endpoint for probes."""
    return {"status": "ok"}
