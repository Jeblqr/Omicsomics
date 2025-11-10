"""
Compression format handler for various compression formats.

Supports:
- gzip (.gz) - Standard compression
- bgzip (.gz) - Block-gzip for genomics (with tabix indexing)
- zip (.zip) - Archive compression
- tar (.tar, .tar.gz, .tar.bz2) - Archive bundling

Dependencies:
- Built-in: gzip, zipfile, tarfile
- Optional: pysam (for bgzip)
"""
from pathlib import Path
from typing import Dict, Optional, List
import gzip
import zipfile
import tarfile
import shutil


class CompressionHandler:
    """
    Handler for compression and archive formats.
    
    Singleton pattern for efficiency.
    """
    
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def __init__(self):
        """Initialize compression handler."""
        if not hasattr(self, '_initialized'):
            self._initialized = True
            self._check_dependencies()
    
    def _check_dependencies(self):
        """Check if optional tools are available."""
        try:
            import pysam
            self.pysam = pysam
            self.has_pysam = True
        except ImportError:
            self.has_pysam = False
    
    # ============================================================
    # Gzip Compression
    # ============================================================
    
    def compress_gzip(self, source: str, target: str):
        """
        Compress file with gzip.
        
        Args:
            source: Input file path
            target: Output .gz file path
        """
        source_path = Path(source)
        target_path = Path(target)
        
        with open(source_path, 'rb') as f_in:
            with gzip.open(target_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    
    def decompress_gzip(self, source: str, target: str):
        """
        Decompress gzip file.
        
        Args:
            source: Input .gz file path
            target: Output file path
        """
        source_path = Path(source)
        target_path = Path(target)
        
        with gzip.open(source_path, 'rb') as f_in:
            with open(target_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    
    # ============================================================
    # BGzip Compression (Block-gzip for genomics)
    # ============================================================
    
    def compress_bgzip(self, source: str, target: str):
        """
        Compress file with bgzip (block-gzip).
        
        Requires pysam.
        
        Args:
            source: Input file path
            target: Output .gz file path
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for bgzip compression")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Use pysam to bgzip
        self.pysam.tabix_compress(str(source_path), str(target_path), force=True)
    
    def decompress_bgzip(self, source: str, target: str):
        """
        Decompress bgzip file.
        
        Can use standard gzip decompression.
        
        Args:
            source: Input .gz file path
            target: Output file path
        """
        self.decompress_gzip(source, target)
    
    def index_bgzip(self, file_path: str, preset: str = 'gff'):
        """
        Create tabix index for bgzipped file.
        
        Requires pysam.
        
        Args:
            file_path: Bgzipped file path
            preset: Index preset (gff, bed, vcf, sam)
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for tabix indexing")
        
        # Create tabix index
        self.pysam.tabix_index(str(file_path), preset=preset, force=True)
    
    # ============================================================
    # Zip Archive
    # ============================================================
    
    def compress_zip(self, sources: List[str], target: str):
        """
        Compress files into zip archive.
        
        Args:
            sources: List of input file paths
            target: Output .zip file path
        """
        target_path = Path(target)
        
        with zipfile.ZipFile(target_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for source in sources:
                source_path = Path(source)
                if source_path.is_file():
                    zipf.write(source_path, source_path.name)
                elif source_path.is_dir():
                    # Add directory recursively
                    for file in source_path.rglob('*'):
                        if file.is_file():
                            arcname = file.relative_to(source_path.parent)
                            zipf.write(file, arcname)
    
    def decompress_zip(self, source: str, target_dir: str):
        """
        Decompress zip archive.
        
        Args:
            source: Input .zip file path
            target_dir: Output directory path
        """
        source_path = Path(source)
        target_path = Path(target_dir)
        
        target_path.mkdir(parents=True, exist_ok=True)
        
        with zipfile.ZipFile(source_path, 'r') as zipf:
            zipf.extractall(target_path)
    
    def list_zip_contents(self, file_path: str) -> List[Dict]:
        """
        List contents of zip archive.
        
        Args:
            file_path: Zip file path
            
        Returns:
            List of file information dictionaries
        """
        source_path = Path(file_path)
        contents = []
        
        with zipfile.ZipFile(source_path, 'r') as zipf:
            for info in zipf.infolist():
                contents.append({
                    'filename': info.filename,
                    'file_size': info.file_size,
                    'compress_size': info.compress_size,
                    'compress_type': info.compress_type,
                    'is_dir': info.is_dir()
                })
        
        return contents
    
    # ============================================================
    # Tar Archive
    # ============================================================
    
    def compress_tar(self, sources: List[str], target: str,
                    compression: Optional[str] = None):
        """
        Compress files into tar archive.
        
        Args:
            sources: List of input file paths
            target: Output .tar file path
            compression: Compression type (None, 'gz', 'bz2', 'xz')
        """
        target_path = Path(target)
        
        # Determine mode
        if compression == 'gz':
            mode = 'w:gz'
        elif compression == 'bz2':
            mode = 'w:bz2'
        elif compression == 'xz':
            mode = 'w:xz'
        else:
            mode = 'w'
        
        with tarfile.open(target_path, mode) as tar:
            for source in sources:
                source_path = Path(source)
                if source_path.exists():
                    tar.add(source_path, arcname=source_path.name)
    
    def decompress_tar(self, source: str, target_dir: str):
        """
        Decompress tar archive.
        
        Args:
            source: Input .tar file path
            target_dir: Output directory path
        """
        source_path = Path(source)
        target_path = Path(target_dir)
        
        target_path.mkdir(parents=True, exist_ok=True)
        
        with tarfile.open(source_path, 'r:*') as tar:
            tar.extractall(target_path)
    
    def list_tar_contents(self, file_path: str) -> List[Dict]:
        """
        List contents of tar archive.
        
        Args:
            file_path: Tar file path
            
        Returns:
            List of file information dictionaries
        """
        source_path = Path(file_path)
        contents = []
        
        with tarfile.open(source_path, 'r:*') as tar:
            for member in tar.getmembers():
                contents.append({
                    'name': member.name,
                    'size': member.size,
                    'is_file': member.isfile(),
                    'is_dir': member.isdir(),
                    'is_link': member.issym() or member.islnk(),
                    'mode': oct(member.mode),
                    'mtime': member.mtime
                })
        
        return contents
    
    # ============================================================
    # Auto-detect Compression Type
    # ============================================================
    
    def detect_compression(self, file_path: str) -> str:
        """
        Auto-detect compression type from file.
        
        Args:
            file_path: File path
            
        Returns:
            Compression type ('gzip', 'zip', 'tar', 'tar.gz', 'tar.bz2', 'none')
        """
        source_path = Path(file_path)
        
        # Check by extension first
        if source_path.suffix == '.gz':
            if source_path.stem.endswith('.tar'):
                return 'tar.gz'
            else:
                return 'gzip'
        elif source_path.suffix == '.bz2':
            if source_path.stem.endswith('.tar'):
                return 'tar.bz2'
            else:
                return 'bzip2'
        elif source_path.suffix == '.xz':
            if source_path.stem.endswith('.tar'):
                return 'tar.xz'
            else:
                return 'xz'
        elif source_path.suffix == '.zip':
            return 'zip'
        elif source_path.suffix == '.tar':
            return 'tar'
        
        # Check by magic bytes
        with open(source_path, 'rb') as f:
            magic = f.read(10)
            
            # Gzip magic: 1f 8b
            if magic[:2] == b'\x1f\x8b':
                return 'gzip'
            
            # Zip magic: 50 4b
            elif magic[:2] == b'\x50\x4b':
                return 'zip'
            
            # Tar magic: "ustar" at offset 257
            elif b'ustar' in magic:
                return 'tar'
        
        return 'none'
    
    def decompress_auto(self, source: str, target: str):
        """
        Auto-detect compression type and decompress.
        
        Args:
            source: Input compressed file path
            target: Output file/directory path
        """
        compression = self.detect_compression(source)
        
        if compression == 'gzip':
            self.decompress_gzip(source, target)
        elif compression == 'zip':
            self.decompress_zip(source, target)
        elif compression in ['tar', 'tar.gz', 'tar.bz2', 'tar.xz']:
            self.decompress_tar(source, target)
        elif compression == 'none':
            raise ValueError(f"File {source} is not compressed")
        else:
            raise ValueError(f"Unsupported compression type: {compression}")
    
    # ============================================================
    # Statistics
    # ============================================================
    
    def get_compression_stats(self, file_path: str) -> Dict:
        """
        Get compression statistics.
        
        Args:
            file_path: Compressed file path
            
        Returns:
            Dictionary with statistics
        """
        source_path = Path(file_path)
        compression = self.detect_compression(file_path)
        
        stats = {
            'compression_type': compression,
            'compressed_size': source_path.stat().st_size
        }
        
        # Get detailed stats for zip
        if compression == 'zip':
            with zipfile.ZipFile(source_path, 'r') as zipf:
                total_uncompressed = sum(info.file_size for info in zipf.infolist())
                stats['uncompressed_size'] = total_uncompressed
                stats['num_files'] = len(zipf.infolist())
                stats['compression_ratio'] = total_uncompressed / stats['compressed_size'] if stats['compressed_size'] > 0 else 0
        
        # Get detailed stats for tar
        elif compression in ['tar', 'tar.gz', 'tar.bz2', 'tar.xz']:
            with tarfile.open(source_path, 'r:*') as tar:
                members = tar.getmembers()
                total_uncompressed = sum(m.size for m in members if m.isfile())
                stats['uncompressed_size'] = total_uncompressed
                stats['num_files'] = len([m for m in members if m.isfile()])
                stats['num_dirs'] = len([m for m in members if m.isdir()])
                stats['compression_ratio'] = total_uncompressed / stats['compressed_size'] if stats['compressed_size'] > 0 else 0
        
        return stats


def get_compression_handler() -> CompressionHandler:
    """Get singleton instance of CompressionHandler."""
    return CompressionHandler()
