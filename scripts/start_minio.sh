#!/bin/bash
# Start MinIO server

cd "$(dirname "$0")/.."

export MINIO_ROOT_USER=minioadmin
export MINIO_ROOT_PASSWORD=minioadmin123

./bin/minio server local_minio_data \
  --address "127.0.0.1:9002" \
  --console-address "127.0.0.1:9003"
