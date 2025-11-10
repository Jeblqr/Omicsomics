from fastapi import APIRouter

from . import (
    auth,
    converters,
    custom_pipelines,
    epigenomics,
    files,
    formats,
    genomics,
    gwas,
    metabolomics,
    multiomics,
    pipelines,
    projects,
    data,
    runs,
    proteomics,
    qc,
    samples,
    singlecell,
    tools,
    transcriptomics,
    users,
    visualizations,
    workflows,
)

# Import additional API modules from parent package
from .. import (
    data_merger,
    quick_visualizer,
    interactive_conversion,
    datasets,
    pipeline_datasets,
    visualization_workspace,
    data_export,
    data_editing,
    custom_scripts,
)

api_router = APIRouter()

api_router.include_router(auth.router, prefix="/auth", tags=["auth"])
api_router.include_router(converters.router, prefix="/converters", tags=["converters"])
api_router.include_router(
    formats.router, tags=["formats"]
)  # No prefix, already has /api/formats
api_router.include_router(projects.router, prefix="/projects", tags=["projects"])
api_router.include_router(samples.router, prefix="/samples", tags=["samples"])
api_router.include_router(files.router, prefix="/files", tags=["files"])
api_router.include_router(data.router, prefix="/data", tags=["data"])
api_router.include_router(runs.router, prefix="/runs", tags=["runs"])
api_router.include_router(pipelines.router, prefix="/pipelines", tags=["pipelines"])
api_router.include_router(
    custom_pipelines.router, prefix="/custom-pipelines", tags=["custom-pipelines"]
)
api_router.include_router(tools.router, prefix="/tools", tags=["tools"])
api_router.include_router(workflows.router, prefix="/workflows", tags=["workflows"])
api_router.include_router(qc.router, prefix="/qc", tags=["qc"])
api_router.include_router(genomics.router, prefix="/genomics", tags=["genomics"])
api_router.include_router(gwas.router, prefix="/gwas", tags=["gwas"])
api_router.include_router(
    transcriptomics.router, prefix="/transcriptomics", tags=["transcriptomics"]
)
api_router.include_router(singlecell.router, prefix="/singlecell", tags=["singlecell"])
api_router.include_router(
    epigenomics.router, prefix="/epigenomics", tags=["epigenomics"]
)
api_router.include_router(proteomics.router, prefix="/proteomics", tags=["proteomics"])
api_router.include_router(
    metabolomics.router, prefix="/metabolomics", tags=["metabolomics"]
)
api_router.include_router(multiomics.router, prefix="/multiomics", tags=["multiomics"])
api_router.include_router(
    visualizations.router, prefix="/visualizations", tags=["visualizations"]
)
api_router.include_router(users.router, prefix="/users", tags=["users"])

# Additional routers from parent api package
api_router.include_router(
    data_merger.router, tags=["data-merger"]
)  # Already has /api/data-merger prefix
api_router.include_router(
    quick_visualizer.router, tags=["quick-visualizer"]
)  # Already has /api/quick-visualizer prefix
api_router.include_router(
    interactive_conversion.router, tags=["interactive-conversion"]
)  # Already has /api/interactive-conversion prefix
api_router.include_router(
    datasets.router, tags=["datasets"]
)  # Already has /api/datasets prefix
api_router.include_router(
    pipeline_datasets.router, tags=["pipeline-datasets"]
)  # Already has /api/pipeline-datasets prefix
api_router.include_router(
    visualization_workspace.router, tags=["visualization-workspace"]
)  # Already has /api/viz-workspace prefix
api_router.include_router(
    data_export.router, tags=["data-export"]
)  # Already has /api/data-export prefix
api_router.include_router(
    data_editing.router, tags=["data-editing"]
)  # Already has /api/data-editing prefix
api_router.include_router(
    custom_scripts.router, tags=["custom-scripts"]
)  # Already has /api/custom-scripts prefix
