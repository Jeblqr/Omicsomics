#!/usr/bin/env python3
"""Initialize MinIO buckets for Omicsomics."""

from minio import Minio
from minio.error import S3Error

# MinIO connection settings
endpoint = "localhost:9002"
access_key = "minioadmin"
secret_key = "minioadmin123"
bucket_name = "omicsomics"

# Create MinIO client
client = Minio(endpoint, access_key=access_key, secret_key=secret_key, secure=False)

# Create bucket if it doesn't exist
try:
    if not client.bucket_exists(bucket_name):
        client.make_bucket(bucket_name)
        print(f"✓ Created bucket: {bucket_name}")
    else:
        print(f"✓ Bucket already exists: {bucket_name}")
except S3Error as e:
    print(f"✗ Error: {e}")
    exit(1)

print("✓ MinIO initialization complete!")
