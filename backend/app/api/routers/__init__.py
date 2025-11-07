from fastapi import APIRouter

from . import (
    auth,
    epigenomics,
    files,
    genomics,
    metabolomics,
    projects,
    proteomics,
    qc,
    samples,
    singlecell,
    transcriptomics,
    visualizations,
    workflows,
)

api_router = APIRouter()

api_router.include_router(auth.router, tags=["auth"])
api_router.include_router(projects.router, prefix="/projects", tags=["projects"])
api_router.include_router(samples.router, prefix="/samples", tags=["samples"])
api_router.include_router(files.router, prefix="/files", tags=["files"])
api_router.include_router(workflows.router, prefix="/workflows", tags=["workflows"])
api_router.include_router(qc.router, prefix="/qc", tags=["qc"])
api_router.include_router(genomics.router, prefix="/genomics", tags=["genomics"])
api_router.include_router(
    transcriptomics.router, prefix="/transcriptomics", tags=["transcriptomics"]
)
api_router.include_router(singlecell.router, prefix="/singlecell", tags=["singlecell"])
api_router.include_router(
    epigenomics.router, prefix="/epigenomics", tags=["epigenomics"]
)
api_router.include_router(
    proteomics.router, prefix="/proteomics", tags=["proteomics"]
)
api_router.include_router(
    metabolomics.router, prefix="/metabolomics", tags=["metabolomics"]
)
api_router.include_router(
    visualizations.router, prefix="/visualizations", tags=["visualizations"]
)
