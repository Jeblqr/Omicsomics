FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Copy and install Python dependencies
COPY pyproject.toml ./
RUN pip install --upgrade pip \
    && pip install "setuptools>=68" "wheel" \
    && pip install -e .[dev]

# Copy application code
COPY app ./app
COPY alembic ./alembic
COPY alembic.ini ./

# Create necessary directories
RUN mkdir -p /app/data /app/logs

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8001", "--reload"]
