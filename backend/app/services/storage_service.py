import base64
import json
import io
from typing import BinaryIO

from app.storage.s3 import s3_client
from app.core import crypto


async def upload_encrypted_bytes(
    user_id: int,
    project_id: int,
    filename: str,
    content: bytes,
    content_type: str = "application/octet-stream",
) -> tuple[str, dict]:
    """Encrypt content for user and upload to object storage.

    Returns (object_key, metadata_dict)
    """
    # Encrypt
    encrypted_blob = crypto.encrypt_for_user(user_id, content)

    # Prepare metadata JSON package
    metadata = {
        "alg": "AES-GCM",
        "version": "1",
        "filename": filename,
        "content_type": content_type,
        "size": len(content),
        # store blob as base64 for metadata convenience OR store raw object body
    }

    # Generate object key in sandbox path
    import uuid

    object_key = f"sandbox/{project_id}/{user_id}/{uuid.uuid4().hex}.enc"

    # Upload encrypted blob as raw bytes
    await s3_client.upload_file(io.BytesIO(encrypted_blob), object_key, content_type)

    # Attach storage key to metadata
    metadata["object_key"] = object_key

    return object_key, metadata


async def download_and_decrypt(user_id: int, object_key: str) -> bytes:
    """Download encrypted blob and decrypt for user."""
    blob = await s3_client.download_file(object_key)
    plaintext = crypto.decrypt_for_user(user_id, blob)
    return plaintext


async def delete_object(object_key: str) -> None:
    """Delete an object from storage."""
    await s3_client.delete_file(object_key)
