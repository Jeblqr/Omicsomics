"""Tests for encryption/decryption utilities."""

import pytest
from app.core import crypto


def test_encrypt_decrypt_roundtrip():
    """Test that encryption and decryption produce the original plaintext."""
    user_id = 42
    plaintext = b"Hello, secure world!"

    encrypted = crypto.encrypt_for_user(user_id, plaintext)
    decrypted = crypto.decrypt_for_user(user_id, encrypted)

    assert decrypted == plaintext


def test_encrypt_decrypt_with_associated_data():
    """Test encryption with associated data (AEAD)."""
    user_id = 100
    plaintext = b"Sensitive data"
    associated = b"metadata-tag"

    encrypted = crypto.encrypt_for_user(user_id, plaintext, associated_data=associated)
    decrypted = crypto.decrypt_for_user(user_id, encrypted, associated_data=associated)

    assert decrypted == plaintext


def test_decrypt_with_wrong_associated_data_fails():
    """Test that decryption fails if associated data doesn't match."""
    user_id = 50
    plaintext = b"Secret message"
    associated_correct = b"correct-tag"
    associated_wrong = b"wrong-tag"

    encrypted = crypto.encrypt_for_user(
        user_id, plaintext, associated_data=associated_correct
    )

    with pytest.raises(Exception):
        crypto.decrypt_for_user(user_id, encrypted, associated_data=associated_wrong)


def test_different_users_produce_different_keys():
    """Test that different users produce different encrypted outputs."""
    plaintext = b"Same message"
    user1_encrypted = crypto.encrypt_for_user(1, plaintext)
    user2_encrypted = crypto.encrypt_for_user(2, plaintext)

    # Encrypted blobs should differ (different keys)
    assert user1_encrypted != user2_encrypted

    # Each user can decrypt their own
    assert crypto.decrypt_for_user(1, user1_encrypted) == plaintext
    assert crypto.decrypt_for_user(2, user2_encrypted) == plaintext

    # User 1 cannot decrypt user 2's message
    with pytest.raises(Exception):
        crypto.decrypt_for_user(1, user2_encrypted)


def test_generate_user_key_hex():
    """Test that generate_user_key_hex returns a 64-char hex string (32 bytes)."""
    key_hex = crypto.generate_user_key_hex(123)
    assert len(key_hex) == 64
    assert all(c in "0123456789abcdef" for c in key_hex)
