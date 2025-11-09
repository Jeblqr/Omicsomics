import os
import binascii
from typing import Tuple

from cryptography.hazmat.primitives.ciphers.aead import AESGCM
from cryptography.hazmat.primitives.kdf.hkdf import HKDF
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.backends import default_backend

from app.settings import settings


def _master_key_bytes() -> bytes:
    """Return master key bytes derived from settings.MASTER_KEY.

    Expect a hex string (64 hex chars => 32 bytes). Fall back to UTF-8
    bytes if not hex.
    """
    mk = settings.MASTER_KEY
    try:
        # support hex-encoded master key
        return binascii.unhexlify(mk)
    except (binascii.Error, TypeError):
        return mk.encode("utf-8")


def derive_user_key(user_id: int, length: int = 32) -> bytes:
    """Derive a per-user symmetric key from the master key using HKDF.

    Deterministic (same user_id -> same key) but secure as long as
    MASTER_KEY remains secret. Returns raw bytes suitable for AESGCM.
    """
    master = _master_key_bytes()
    hkdf = HKDF(
        algorithm=hashes.SHA256(),
        length=length,
        salt=None,
        info=str(user_id).encode("utf-8"),
        backend=default_backend(),
    )
    return hkdf.derive(master)


def encrypt_for_user(
    user_id: int, plaintext: bytes, associated_data: bytes | None = None
) -> bytes:
    """Encrypt bytes with user-derived key using AES-GCM.

    Returns: nonce || ciphertext || tag (raw bytes). Caller can store this
    blob as-is (e.g. in object storage). Nonce is 12 bytes.
    """
    key = derive_user_key(user_id)
    aesgcm = AESGCM(key)
    nonce = os.urandom(12)
    ct = aesgcm.encrypt(nonce, plaintext, associated_data)
    return nonce + ct


def decrypt_for_user(
    user_id: int, blob: bytes, associated_data: bytes | None = None
) -> bytes:
    """Decrypt a blob produced by encrypt_for_user.

    Expects nonce (12 bytes) prefix.
    """
    if len(blob) < 13:
        raise ValueError("Encrypted blob too short")
    nonce = blob[:12]
    ct = blob[12:]
    key = derive_user_key(user_id)
    aesgcm = AESGCM(key)
    return aesgcm.decrypt(nonce, ct, associated_data)


def generate_user_key_hex(user_id: int) -> str:
    """Helper: derive per-user key and return hex (useful for debugging)."""
    return binascii.hexlify(derive_user_key(user_id)).decode("ascii")
