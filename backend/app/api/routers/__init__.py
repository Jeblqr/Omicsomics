from fastapi import APIRouter

from . import auth, files, projects, qc, samples, workflows

api_router = APIRouter()

api_router.include_router(auth.router, tags=["auth"])
api_router.include_router(projects.router, prefix="/projects", tags=["projects"])
api_router.include_router(samples.router, prefix="/samples", tags=["samples"])
api_router.include_router(files.router, prefix="/files", tags=["files"])
api_router.include_router(workflows.router, prefix="/workflows", tags=["workflows"])
api_router.include_router(qc.router, prefix="/qc", tags=["qc"])
api_router.include_router(files.router, prefix="/files", tags=["files"])
