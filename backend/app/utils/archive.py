"""
Archive file handling utilities.

Provides safe extraction and preview functionality for compressed archives
(ZIP, TAR, TAR.GZ, TAR.BZ2) with path traversal protection.
"""

import logging
import mimetypes
import tarfile
import zipfile
from pathlib import Path
from typing import List, Dict, Any, Optional
import os

logger = logging.getLogger(__name__)


class ArchiveError(Exception):
    """Base exception for archive operations."""
    pass


class UnsafeArchiveError(ArchiveError):
    """Raised when archive contains unsafe paths (path traversal attempt)."""
    pass


class ArchiveFile:
    """Information about a file in an archive."""
    
    def __init__(
        self,
        name: str,
        path: str,
        size: int,
        is_dir: bool = False,
        mime_type: Optional[str] = None
    ):
        self.name = name
        self.path = path  # Relative path within archive
        self.size = size
        self.is_dir = is_dir
        self.mime_type = mime_type or self._guess_mime_type()
    
    def _guess_mime_type(self) -> str:
        """Guess MIME type from file extension."""
        if self.is_dir:
            return "inode/directory"
        
        mime_type, _ = mimetypes.guess_type(self.name)
        return mime_type or "application/octet-stream"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "path": self.path,
            "size": self.size,
            "is_dir": self.is_dir,
            "mime_type": self.mime_type,
        }


def is_safe_path(base_path: Path, target_path: Path) -> bool:
    """
    Check if target path is safely within base path.
    
    Prevents path traversal attacks like "../../../etc/passwd".
    
    Args:
        base_path: The base directory that should contain the target
        target_path: The path to check
    
    Returns:
        True if target is safely within base, False otherwise
    """
    # Resolve to absolute paths
    base_abs = base_path.resolve()
    target_abs = target_path.resolve()
    
    # Check if target is within base
    try:
        target_abs.relative_to(base_abs)
        return True
    except ValueError:
        return False


def sanitize_archive_member_name(name: str) -> str:
    """
    Sanitize archive member name to prevent path traversal.
    
    Args:
        name: Original member name from archive
    
    Returns:
        Sanitized name safe for extraction
    
    Raises:
        UnsafeArchiveError: If name contains unsafe patterns
    """
    # Remove leading slashes and drive letters (Windows)
    name = name.lstrip("/\\")
    if len(name) > 1 and name[1] == ":":
        name = name[2:].lstrip("/\\")
    
    # Check for parent directory references
    parts = Path(name).parts
    if ".." in parts:
        raise UnsafeArchiveError(f"Archive contains unsafe path: {name}")
    
    # Check for absolute paths
    if Path(name).is_absolute():
        raise UnsafeArchiveError(f"Archive contains absolute path: {name}")
    
    return name


def list_zip_contents(zip_path: Path) -> List[ArchiveFile]:
    """
    List contents of a ZIP file.
    
    Args:
        zip_path: Path to ZIP file
    
    Returns:
        List of ArchiveFile objects
    
    Raises:
        ArchiveError: If ZIP is invalid or cannot be read
    """
    files = []
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as zf:
            for info in zf.infolist():
                # Skip directories if they don't contain files
                if info.is_dir():
                    continue
                
                try:
                    safe_name = sanitize_archive_member_name(info.filename)
                except UnsafeArchiveError as e:
                    logger.warning(f"Skipping unsafe file in ZIP: {e}")
                    continue
                
                files.append(ArchiveFile(
                    name=Path(safe_name).name,
                    path=safe_name,
                    size=info.file_size,
                    is_dir=info.is_dir()
                ))
    except zipfile.BadZipFile as e:
        raise ArchiveError(f"Invalid ZIP file: {e}")
    except Exception as e:
        raise ArchiveError(f"Error reading ZIP file: {e}")
    
    return files


def list_tar_contents(tar_path: Path) -> List[ArchiveFile]:
    """
    List contents of a TAR file (including .tar.gz, .tar.bz2).
    
    Args:
        tar_path: Path to TAR file
    
    Returns:
        List of ArchiveFile objects
    
    Raises:
        ArchiveError: If TAR is invalid or cannot be read
    """
    files = []
    
    try:
        with tarfile.open(tar_path, 'r:*') as tf:
            for member in tf.getmembers():
                # Skip directories
                if member.isdir():
                    continue
                
                try:
                    safe_name = sanitize_archive_member_name(member.name)
                except UnsafeArchiveError as e:
                    logger.warning(f"Skipping unsafe file in TAR: {e}")
                    continue
                
                files.append(ArchiveFile(
                    name=Path(safe_name).name,
                    path=safe_name,
                    size=member.size,
                    is_dir=member.isdir()
                ))
    except tarfile.TarError as e:
        raise ArchiveError(f"Invalid TAR file: {e}")
    except Exception as e:
        raise ArchiveError(f"Error reading TAR file: {e}")
    
    return files


def list_archive_contents(archive_path: Path) -> List[ArchiveFile]:
    """
    List contents of an archive file.
    
    Automatically detects archive type from extension.
    
    Args:
        archive_path: Path to archive file
    
    Returns:
        List of ArchiveFile objects
    
    Raises:
        ArchiveError: If archive type is unsupported or cannot be read
    """
    suffix = archive_path.suffix.lower()
    
    if suffix == '.zip':
        return list_zip_contents(archive_path)
    elif suffix in ['.tar', '.gz', '.bz2', '.tgz', '.tbz']:
        # Handle .tar.gz, .tar.bz2, etc.
        if archive_path.name.endswith('.tar.gz') or archive_path.name.endswith('.tar.bz2'):
            return list_tar_contents(archive_path)
        elif suffix in ['.tar', '.tgz', '.tbz']:
            return list_tar_contents(archive_path)
        else:
            raise ArchiveError(f"Unsupported archive format: {suffix}")
    else:
        raise ArchiveError(f"Unsupported archive format: {suffix}")


def extract_zip_file(
    zip_path: Path,
    extract_to: Path,
    member_name: Optional[str] = None
) -> Path:
    """
    Extract ZIP file or specific member from ZIP.
    
    Args:
        zip_path: Path to ZIP file
        extract_to: Directory to extract to
        member_name: Optional specific file to extract. If None, extracts all.
    
    Returns:
        Path to extraction directory
    
    Raises:
        UnsafeArchiveError: If archive contains unsafe paths
        ArchiveError: If extraction fails
    """
    extract_to.mkdir(parents=True, exist_ok=True)
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as zf:
            if member_name:
                # Extract specific file
                safe_name = sanitize_archive_member_name(member_name)
                
                # Verify member exists
                if member_name not in zf.namelist():
                    raise ArchiveError(f"File not found in archive: {member_name}")
                
                # Extract with safe path
                target_path = extract_to / safe_name
                target_path.parent.mkdir(parents=True, exist_ok=True)
                
                # Extract member
                with zf.open(member_name) as source, open(target_path, 'wb') as target:
                    target.write(source.read())
                
                return target_path
            else:
                # Extract all files
                for member in zf.namelist():
                    try:
                        safe_name = sanitize_archive_member_name(member)
                    except UnsafeArchiveError as e:
                        logger.warning(f"Skipping unsafe file: {e}")
                        continue
                    
                    target_path = extract_to / safe_name
                    
                    # Verify safety
                    if not is_safe_path(extract_to, target_path):
                        raise UnsafeArchiveError(f"Path traversal attempt: {member}")
                    
                    # Extract
                    target_path.parent.mkdir(parents=True, exist_ok=True)
                    with zf.open(member) as source, open(target_path, 'wb') as target:
                        target.write(source.read())
                
                return extract_to
    except zipfile.BadZipFile as e:
        raise ArchiveError(f"Invalid ZIP file: {e}")
    except UnsafeArchiveError:
        raise
    except Exception as e:
        raise ArchiveError(f"Error extracting ZIP: {e}")


def extract_tar_file(
    tar_path: Path,
    extract_to: Path,
    member_name: Optional[str] = None
) -> Path:
    """
    Extract TAR file or specific member from TAR.
    
    Args:
        tar_path: Path to TAR file
        extract_to: Directory to extract to
        member_name: Optional specific file to extract. If None, extracts all.
    
    Returns:
        Path to extraction directory or extracted file
    
    Raises:
        UnsafeArchiveError: If archive contains unsafe paths
        ArchiveError: If extraction fails
    """
    extract_to.mkdir(parents=True, exist_ok=True)
    
    try:
        with tarfile.open(tar_path, 'r:*') as tf:
            if member_name:
                # Extract specific file
                safe_name = sanitize_archive_member_name(member_name)
                
                # Find member
                member = None
                for m in tf.getmembers():
                    if m.name == member_name:
                        member = m
                        break
                
                if member is None:
                    raise ArchiveError(f"File not found in archive: {member_name}")
                
                # Extract
                target_path = extract_to / safe_name
                target_path.parent.mkdir(parents=True, exist_ok=True)
                
                with tf.extractfile(member) as source:
                    if source:
                        with open(target_path, 'wb') as target:
                            target.write(source.read())
                
                return target_path
            else:
                # Extract all files with safety checks
                for member in tf.getmembers():
                    try:
                        safe_name = sanitize_archive_member_name(member.name)
                    except UnsafeArchiveError as e:
                        logger.warning(f"Skipping unsafe file: {e}")
                        continue
                    
                    target_path = extract_to / safe_name
                    
                    # Verify safety
                    if not is_safe_path(extract_to, target_path):
                        raise UnsafeArchiveError(f"Path traversal attempt: {member.name}")
                    
                    # Extract
                    if member.isfile():
                        target_path.parent.mkdir(parents=True, exist_ok=True)
                        with tf.extractfile(member) as source:
                            if source:
                                with open(target_path, 'wb') as target:
                                    target.write(source.read())
                
                return extract_to
    except tarfile.TarError as e:
        raise ArchiveError(f"Invalid TAR file: {e}")
    except UnsafeArchiveError:
        raise
    except Exception as e:
        raise ArchiveError(f"Error extracting TAR: {e}")


def extract_archive(
    archive_path: Path,
    extract_to: Path,
    member_name: Optional[str] = None
) -> Path:
    """
    Extract archive file or specific member.
    
    Automatically detects archive type from extension.
    
    Args:
        archive_path: Path to archive file
        extract_to: Directory to extract to
        member_name: Optional specific file to extract. If None, extracts all.
    
    Returns:
        Path to extraction directory or extracted file
    
    Raises:
        UnsafeArchiveError: If archive contains unsafe paths
        ArchiveError: If extraction fails
    """
    suffix = archive_path.suffix.lower()
    
    if suffix == '.zip':
        return extract_zip_file(archive_path, extract_to, member_name)
    elif suffix in ['.tar', '.gz', '.bz2', '.tgz', '.tbz']:
        if archive_path.name.endswith('.tar.gz') or archive_path.name.endswith('.tar.bz2'):
            return extract_tar_file(archive_path, extract_to, member_name)
        elif suffix in ['.tar', '.tgz', '.tbz']:
            return extract_tar_file(archive_path, extract_to, member_name)
        else:
            raise ArchiveError(f"Unsupported archive format: {suffix}")
    else:
        raise ArchiveError(f"Unsupported archive format: {suffix}")


def is_archive_file(file_path: Path) -> bool:
    """
    Check if file is a supported archive format.
    
    Args:
        file_path: Path to file
    
    Returns:
        True if file is a supported archive
    """
    suffix = file_path.suffix.lower()
    name_lower = file_path.name.lower()
    
    return (
        suffix in ['.zip', '.tar', '.tgz', '.tbz'] or
        name_lower.endswith('.tar.gz') or
        name_lower.endswith('.tar.bz2')
    )
