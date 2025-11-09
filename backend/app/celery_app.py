"""
Celery application configuration for async task processing
"""

from celery import Celery
import os

# Get broker URL from environment
CELERY_BROKER_URL = os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0")
CELERY_RESULT_BACKEND = os.getenv("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

# Create Celery app
celery_app = Celery(
    "omicsomics",
    broker=CELERY_BROKER_URL,
    backend=CELERY_RESULT_BACKEND,
    include=["app.tasks.file_tasks"],
)

# Configure Celery
celery_app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="UTC",
    enable_utc=True,
    task_track_started=True,
    task_time_limit=3600,  # 1 hour max
    task_soft_time_limit=3300,  # 55 minutes soft limit
    worker_prefetch_multiplier=1,
    worker_max_tasks_per_child=50,
)

# Task routing (可以为不同类型的任务设置不同的队列)
celery_app.conf.task_routes = {
    "app.tasks.file_tasks.process_file_async": {"queue": "file_processing"},
    "app.tasks.file_tasks.process_large_file_async": {"queue": "large_files"},
}

# 设置默认队列
celery_app.conf.task_default_queue = "file_processing"

if __name__ == "__main__":
    celery_app.start()
