"""
Tests for archive utilities.

Tests safe extraction, path traversal protection, and archive listing.
"""

import tempfile
import zipfile
import tarfile
from pathlib import Path
import pytest

from app.utils.archive import (
    list_archive_contents,
    extract_archive,
    is_archive_file,
    is_safe_path,
    sanitize_archive_member_name,
    ArchiveError,
    UnsafeArchiveError,
)


@pytest.fixture
def temp_dir():
    """Create temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def safe_zip_file(temp_dir):
    """Create a safe ZIP file for testing."""
    zip_path = temp_dir / "safe.zip"

    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("file1.txt", "Content 1")
        zf.writestr("dir/file2.csv", "col1,col2\nval1,val2")
        zf.writestr("data.json", '{"key": "value"}')

    return zip_path


@pytest.fixture
def unsafe_zip_file(temp_dir):
    """Create a ZIP file with path traversal attempt."""
    zip_path = temp_dir / "unsafe.zip"

    with zipfile.ZipFile(zip_path, "w") as zf:
        # This would try to escape the extraction directory
        zf.writestr("../../../etc/passwd", "malicious content")
        zf.writestr("normal_file.txt", "normal content")

    return zip_path


@pytest.fixture
def safe_tar_file(temp_dir):
    """Create a safe TAR.GZ file for testing."""
    tar_path = temp_dir / "safe.tar.gz"

    with tarfile.open(tar_path, "w:gz") as tf:
        # Create temp files to add
        file1 = temp_dir / "temp1.txt"
        file1.write_text("Content 1")
        tf.add(file1, arcname="file1.txt")

        file2 = temp_dir / "temp2.csv"
        file2.write_text("col1,col2\nval1,val2")
        tf.add(file2, arcname="dir/file2.csv")

    return tar_path


class TestPathSafety:
    """Test path safety checks."""

    def test_is_safe_path_within_base(self, temp_dir):
        """Test that path within base is safe."""
        base = temp_dir / "base"
        target = base / "subdir" / "file.txt"

        assert is_safe_path(base, target)

    def test_is_safe_path_outside_base(self, temp_dir):
        """Test that path outside base is unsafe."""
        base = temp_dir / "base"
        target = temp_dir / "other" / "file.txt"

        assert not is_safe_path(base, target)

    def test_is_safe_path_traversal_attempt(self, temp_dir):
        """Test that path traversal is detected."""
        base = temp_dir / "base"
        target = base / ".." / ".." / "etc" / "passwd"

        assert not is_safe_path(base, target)

    def test_sanitize_normal_name(self):
        """Test sanitizing normal file name."""
        name = "dir/subdir/file.txt"
        result = sanitize_archive_member_name(name)
        assert result == name

    def test_sanitize_parent_reference(self):
        """Test that parent directory reference raises error."""
        name = "../etc/passwd"

        with pytest.raises(UnsafeArchiveError):
            sanitize_archive_member_name(name)

    def test_sanitize_absolute_path(self):
        """Test that absolute path raises error."""
        name = "/etc/passwd"

        with pytest.raises(UnsafeArchiveError):
            sanitize_archive_member_name(name)

    def test_sanitize_windows_drive(self):
        """Test sanitizing Windows drive letter."""
        name = "C:/Windows/System32/file.dll"
        result = sanitize_archive_member_name(name)
        assert not result.startswith("C:")
        assert "Windows" in result


class TestArchiveDetection:
    """Test archive file detection."""

    def test_is_archive_zip(self, temp_dir):
        """Test ZIP detection."""
        assert is_archive_file(temp_dir / "file.zip")

    def test_is_archive_tar(self, temp_dir):
        """Test TAR detection."""
        assert is_archive_file(temp_dir / "file.tar")

    def test_is_archive_tar_gz(self, temp_dir):
        """Test TAR.GZ detection."""
        assert is_archive_file(temp_dir / "file.tar.gz")

    def test_is_archive_tar_bz2(self, temp_dir):
        """Test TAR.BZ2 detection."""
        assert is_archive_file(temp_dir / "file.tar.bz2")

    def test_is_archive_tgz(self, temp_dir):
        """Test TGZ detection."""
        assert is_archive_file(temp_dir / "file.tgz")

    def test_not_archive(self, temp_dir):
        """Test non-archive file."""
        assert not is_archive_file(temp_dir / "file.txt")
        assert not is_archive_file(temp_dir / "file.csv")


class TestZipListing:
    """Test ZIP file listing."""

    def test_list_safe_zip(self, safe_zip_file):
        """Test listing safe ZIP contents."""
        contents = list_archive_contents(safe_zip_file)

        assert len(contents) == 3

        # Check file names
        names = {f.name for f in contents}
        assert "file1.txt" in names
        assert "file2.csv" in names
        assert "data.json" in names

        # Check paths
        paths = {f.path for f in contents}
        assert "file1.txt" in paths
        assert "dir/file2.csv" in paths

    def test_list_unsafe_zip_skips_bad_files(self, unsafe_zip_file):
        """Test that unsafe files are skipped in listing."""
        contents = list_archive_contents(unsafe_zip_file)

        # Should only contain the safe file
        assert len(contents) == 1
        assert contents[0].name == "normal_file.txt"


class TestZipExtraction:
    """Test ZIP file extraction."""

    def test_extract_full_zip(self, safe_zip_file, temp_dir):
        """Test extracting entire ZIP."""
        extract_dir = temp_dir / "extracted"
        result = extract_archive(safe_zip_file, extract_dir)

        assert result == extract_dir
        assert (extract_dir / "file1.txt").exists()
        assert (extract_dir / "dir" / "file2.csv").exists()
        assert (extract_dir / "data.json").exists()

        # Verify content
        content = (extract_dir / "file1.txt").read_text()
        assert content == "Content 1"

    def test_extract_specific_file_from_zip(self, safe_zip_file, temp_dir):
        """Test extracting specific file from ZIP."""
        extract_dir = temp_dir / "extracted"
        result = extract_archive(safe_zip_file, extract_dir, "dir/file2.csv")

        assert result == extract_dir / "dir" / "file2.csv"
        assert result.exists()

        # Verify only this file was extracted
        assert not (extract_dir / "file1.txt").exists()
        assert not (extract_dir / "data.json").exists()

    def test_extract_unsafe_zip_raises_error(self, unsafe_zip_file, temp_dir):
        """Test that unsafe ZIP extraction is blocked."""
        extract_dir = temp_dir / "extracted"

        # Should raise error due to path traversal
        with pytest.raises(UnsafeArchiveError):
            extract_archive(unsafe_zip_file, extract_dir)


class TestTarExtraction:
    """Test TAR file extraction."""

    def test_extract_full_tar(self, safe_tar_file, temp_dir):
        """Test extracting entire TAR."""
        extract_dir = temp_dir / "extracted"
        result = extract_archive(safe_tar_file, extract_dir)

        assert result == extract_dir
        assert (extract_dir / "file1.txt").exists()
        assert (extract_dir / "dir" / "file2.csv").exists()

        # Verify content
        content = (extract_dir / "file1.txt").read_text()
        assert content == "Content 1"

    def test_extract_specific_file_from_tar(self, safe_tar_file, temp_dir):
        """Test extracting specific file from TAR."""
        extract_dir = temp_dir / "extracted"
        result = extract_archive(safe_tar_file, extract_dir, "file1.txt")

        assert result == extract_dir / "file1.txt"
        assert result.exists()

        # Verify only this file was extracted
        assert not (extract_dir / "dir" / "file2.csv").exists()


class TestErrorHandling:
    """Test error handling."""

    def test_invalid_zip_raises_error(self, temp_dir):
        """Test that invalid ZIP raises error."""
        fake_zip = temp_dir / "fake.zip"
        fake_zip.write_text("not a zip file")

        with pytest.raises(ArchiveError):
            list_archive_contents(fake_zip)

    def test_unsupported_format_raises_error(self, temp_dir):
        """Test that unsupported format raises error."""
        unsupported = temp_dir / "file.rar"

        with pytest.raises(ArchiveError):
            list_archive_contents(unsupported)

    def test_extract_nonexistent_member_raises_error(self, safe_zip_file, temp_dir):
        """Test extracting non-existent member raises error."""
        extract_dir = temp_dir / "extracted"

        with pytest.raises(ArchiveError):
            extract_archive(safe_zip_file, extract_dir, "nonexistent.txt")


class TestArchiveFileInfo:
    """Test ArchiveFile information."""

    def test_archive_file_mime_types(self, safe_zip_file):
        """Test that MIME types are detected correctly."""
        contents = list_archive_contents(safe_zip_file)

        mime_types = {f.name: f.mime_type for f in contents}

        assert mime_types["file1.txt"] == "text/plain"
        assert mime_types["file2.csv"] == "text/csv"
        assert mime_types["data.json"] == "application/json"

    def test_archive_file_to_dict(self, safe_zip_file):
        """Test converting ArchiveFile to dictionary."""
        contents = list_archive_contents(safe_zip_file)
        file_dict = contents[0].to_dict()

        assert "name" in file_dict
        assert "path" in file_dict
        assert "size" in file_dict
        assert "is_dir" in file_dict
        assert "mime_type" in file_dict
