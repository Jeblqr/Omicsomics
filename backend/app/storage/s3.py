"""S3/MinIO storage client."""

import io
from typing import BinaryIO

import aioboto3
from botocore.exceptions import ClientError

from app.settings import settings


class S3Client:
    """Async S3 client for file storage operations."""

    def __init__(self):
        self.session = aioboto3.Session()
        self.endpoint_url = settings.object_storage_endpoint
        self.access_key = settings.object_storage_access_key
        self.secret_key = settings.object_storage_secret_key
        self.bucket_name = settings.object_storage_bucket

    async def ensure_bucket_exists(self):
        """Create bucket if it doesn't exist."""
        async with self.session.client(
            "s3",
            endpoint_url=self.endpoint_url,
            aws_access_key_id=self.access_key,
            aws_secret_access_key=self.secret_key,
        ) as s3:
            try:
                await s3.head_bucket(Bucket=self.bucket_name)
            except ClientError:
                # Bucket doesn't exist, create it
                await s3.create_bucket(Bucket=self.bucket_name)

    async def upload_file(
        self,
        file_data: BinaryIO,
        object_name: str,
        content_type: str = "application/octet-stream",
    ) -> str:
        """
        Upload a file to S3.

        Args:
            file_data: File-like object containing the data
            object_name: S3 object key (path in bucket)
            content_type: MIME type of the file

        Returns:
            The object key of the uploaded file
        """
        await self.ensure_bucket_exists()

        async with self.session.client(
            "s3",
            endpoint_url=self.endpoint_url,
            aws_access_key_id=self.access_key,
            aws_secret_access_key=self.secret_key,
        ) as s3:
            await s3.upload_fileobj(
                file_data,
                self.bucket_name,
                object_name,
                ExtraArgs={"ContentType": content_type},
            )
        return object_name

    async def download_file(self, object_name: str) -> bytes:
        """
        Download a file from S3.

        Args:
            object_name: S3 object key

        Returns:
            File contents as bytes
        """
        async with self.session.client(
            "s3",
            endpoint_url=self.endpoint_url,
            aws_access_key_id=self.access_key,
            aws_secret_access_key=self.secret_key,
        ) as s3:
            response = await s3.get_object(Bucket=self.bucket_name, Key=object_name)
            async with response["Body"] as stream:
                return await stream.read()

    async def delete_file(self, object_name: str):
        """
        Delete a file from S3.

        Args:
            object_name: S3 object key
        """
        async with self.session.client(
            "s3",
            endpoint_url=self.endpoint_url,
            aws_access_key_id=self.access_key,
            aws_secret_access_key=self.secret_key,
        ) as s3:
            await s3.delete_object(Bucket=self.bucket_name, Key=object_name)

    async def get_presigned_url(self, object_name: str, expiration: int = 3600) -> str:
        """
        Generate a presigned URL for file download.

        Args:
            object_name: S3 object key
            expiration: URL expiration time in seconds

        Returns:
            Presigned URL
        """
        async with self.session.client(
            "s3",
            endpoint_url=self.endpoint_url,
            aws_access_key_id=self.access_key,
            aws_secret_access_key=self.secret_key,
        ) as s3:
            url = await s3.generate_presigned_url(
                "get_object",
                Params={"Bucket": self.bucket_name, "Key": object_name},
                ExpiresIn=expiration,
            )
        return url


# Global S3 client instance
s3_client = S3Client()
