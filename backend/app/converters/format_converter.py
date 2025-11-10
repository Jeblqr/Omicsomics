"""
Data format converter for cross-runtime tool integration.

Supports automatic and manual conversion between:
- CSV, TSV, Excel (tabular formats)
- RDS (R data format)
- h5ad (Python AnnData format)
- pickle (Python object serialization)
- JSON (universal format)
"""

import os
import json
import pandas as pd
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


class FormatConverter:
    """
    Handle data format conversions between different tool runtimes.

    Conversion Matrix:
    CSV ↔ TSV ↔ Excel ↔ JSON
    CSV ↔ RDS (via R)
    CSV ↔ h5ad (via Python)
    CSV ↔ pickle (via Python)
    """

    # Supported format mappings
    FORMATS = {
        # Tabular formats
        "csv": {"extensions": [".csv"], "mime_type": "text/csv"},
        "tsv": {
            "extensions": [".tsv", ".txt"],
            "mime_type": "text/tab-separated-values",
        },
        "excel": {
            "extensions": [".xlsx", ".xls"],
            "mime_type": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        },
        "json": {"extensions": [".json"], "mime_type": "application/json"},
        # R/Python data formats
        "rds": {
            "extensions": [".rds", ".RDS"],
            "mime_type": "application/octet-stream",
        },
        "h5ad": {"extensions": [".h5ad"], "mime_type": "application/octet-stream"},
        "pickle": {
            "extensions": [".pkl", ".pickle"],
            "mime_type": "application/octet-stream",
        },
        # Sequence formats
        "fasta": {"extensions": [".fa", ".fasta", ".fna"], "mime_type": "text/plain"},
        "fastq": {"extensions": [".fq", ".fastq"], "mime_type": "text/plain"},
        "fastq.gz": {
            "extensions": [".fastq.gz", ".fq.gz"],
            "mime_type": "application/gzip",
        },
        # Genomic interval formats
        "bed": {"extensions": [".bed"], "mime_type": "text/plain"},
        "bed3": {"extensions": [".bed"], "mime_type": "text/plain"},
        "bed6": {"extensions": [".bed"], "mime_type": "text/plain"},
        "bed12": {"extensions": [".bed"], "mime_type": "text/plain"},
        "bedgraph": {"extensions": [".bg", ".bedgraph"], "mime_type": "text/plain"},
        "bigwig": {
            "extensions": [".bw", ".bigwig"],
            "mime_type": "application/octet-stream",
        },
        # Alignment formats
        "sam": {"extensions": [".sam"], "mime_type": "text/plain"},
        "bam": {"extensions": [".bam"], "mime_type": "application/octet-stream"},
        # Variant formats
        "vcf": {"extensions": [".vcf"], "mime_type": "text/plain"},
        "vcf.gz": {"extensions": [".vcf.gz"], "mime_type": "application/gzip"},
        "bcf": {"extensions": [".bcf"], "mime_type": "application/octet-stream"},
        # Annotation formats
        "gtf": {"extensions": [".gtf"], "mime_type": "text/plain"},
        "gff3": {"extensions": [".gff", ".gff3"], "mime_type": "text/plain"},
        # Expression matrix formats
        "mtx": {"extensions": [".mtx"], "mime_type": "text/plain"},
        "10x_h5": {"extensions": [".h5"], "mime_type": "application/octet-stream"},
        # Python data formats
        "h5": {"extensions": [".h5", ".hdf5"], "mime_type": "application/octet-stream"},
        "npy": {"extensions": [".npy"], "mime_type": "application/octet-stream"},
        "npz": {"extensions": [".npz"], "mime_type": "application/octet-stream"},
        # R data formats
        "rdata": {"extensions": [".RData", ".rda"], "mime_type": "application/octet-stream"},
        # Compression formats
        "gzip": {"extensions": [".gz"], "mime_type": "application/gzip"},
        "bgzip": {"extensions": [".gz"], "mime_type": "application/gzip"},
        "zip": {"extensions": [".zip"], "mime_type": "application/zip"},
        "tar": {"extensions": [".tar"], "mime_type": "application/x-tar"},
        "tar.gz": {"extensions": [".tar.gz", ".tgz"], "mime_type": "application/gzip"},
    }

    # Conversion time estimates (seconds per GB)
    CONVERSION_TIME_ESTIMATES = {
        # Tabular conversions
        ("csv", "tsv"): 2,
        ("csv", "excel"): 15,
        ("csv", "json"): 10,
        ("csv", "rds"): 20,
        ("csv", "h5ad"): 25,
        ("csv", "pickle"): 8,
        ("rds", "csv"): 20,
        ("h5ad", "csv"): 25,
        ("pickle", "csv"): 8,
        # Sequence format conversions
        ("fastq.gz", "fastq"): 5,  # Decompression
        ("fastq", "fastq.gz"): 10,  # Compression
        ("fastq", "fasta"): 3,  # Remove quality
        ("fasta", "tsv"): 2,  # Parse to table
        ("fastq", "tsv"): 3,  # Parse to table with quality
        # Genomic interval conversions
        ("bed", "csv"): 2,
        ("csv", "bed"): 2,
        ("bed", "bedgraph"): 2,
        ("bedgraph", "bed"): 2,
        ("bigwig", "bedgraph"): 15,
        ("bedgraph", "bigwig"): 20,
        # Alignment format conversions
        ("sam", "bam"): 8,  # SAM to BAM compression
        ("bam", "sam"): 6,  # BAM to SAM decompression
        ("sam", "tsv"): 4,  # Parse to table
        ("bam", "tsv"): 5,  # Decompress and parse
        # Variant format conversions
        ("vcf", "bcf"): 10,  # VCF to BCF compression
        ("bcf", "vcf"): 8,  # BCF to VCF decompression
        ("vcf", "tsv"): 5,  # Parse to table
        ("vcf", "bed"): 3,  # Extract positions
        # Annotation format conversions
        ("gtf", "gff3"): 3,  # Format conversion
        ("gff3", "gtf"): 3,  # Format conversion
        ("gtf", "tsv"): 4,  # Parse to table
        ("gff3", "tsv"): 4,  # Parse to table
        ("gtf", "bed"): 3,  # Extract features
        ("gff3", "bed"): 3,  # Extract features
        # Expression matrix conversions
        ("mtx", "csv"): 20,  # Sparse to dense matrix
        ("csv", "mtx"): 15,  # Dense to sparse matrix
        ("10x_h5", "csv"): 25,  # 10X HDF5 to CSV
        ("10x_h5", "h5ad"): 15,  # 10X HDF5 to AnnData
        ("10x_h5", "mtx"): 10,  # 10X HDF5 to MTX
        ("h5ad", "csv"): 20,  # AnnData to CSV
        ("mtx", "h5ad"): 18,  # MTX to AnnData
        # Python data conversions
        ("h5", "csv"): 15,  # HDF5 to CSV
        ("csv", "h5"): 15,  # CSV to HDF5
        ("h5", "npy"): 10,  # HDF5 to NumPy
        ("npy", "csv"): 8,  # NumPy to CSV
        ("csv", "npy"): 8,  # CSV to NumPy
        ("npy", "npz"): 5,  # NumPy to compressed
        ("npz", "npy"): 5,  # Compressed to NumPy
        ("npz", "csv"): 10,  # Compressed NumPy to CSV
        # R data conversions
        ("rdata", "csv"): 20,  # RData to CSV
        ("rdata", "json"): 18,  # RData to JSON
        # Compression operations
        ("gzip", "decompress"): 5,  # Gzip decompression
        ("bgzip", "decompress"): 5,  # Bgzip decompression
        ("zip", "decompress"): 8,  # Zip decompression
        ("tar", "decompress"): 10,  # Tar extraction
        ("tar.gz", "decompress"): 12,  # Tar.gz extraction
    }

    def __init__(self, storage_path: str = "/tmp/format_conversions"):
        self.storage_path = Path(storage_path)
        self.storage_path.mkdir(parents=True, exist_ok=True)

        # Import specialized converters
        from .sequence_converter import get_sequence_converter
        from .interval_converter import get_interval_converter
        from .alignment_converter import get_alignment_converter
        from .variant_converter import get_variant_converter
        from .annotation_converter import get_annotation_converter
        from .expression_converter import get_expression_converter
        from .python_data_converter import get_python_data_converter
        from .r_data_converter import get_r_data_converter
        from .compression_handler import get_compression_handler

        self.sequence_converter = get_sequence_converter()
        self.interval_converter = get_interval_converter()
        self.alignment_converter = get_alignment_converter()
        self.variant_converter = get_variant_converter()
        self.annotation_converter = get_annotation_converter()
        self.expression_converter = get_expression_converter()
        self.python_data_converter = get_python_data_converter()
        self.r_data_converter = get_r_data_converter()
        self.compression_handler = get_compression_handler()

    def detect_format(self, file_path: str) -> Optional[str]:
        """
        Detect file format from extension or content.

        Args:
            file_path: Path to the file

        Returns:
            Format identifier (csv, rds, h5ad, etc.) or None if unknown
        """
        path = Path(file_path)
        ext = path.suffix.lower()

        for format_name, format_info in self.FORMATS.items():
            if ext in format_info["extensions"]:
                return format_name

        return None

    def get_conversion_path(self, from_format: str, to_format: str) -> List[str]:
        """
        Find conversion path between two formats.

        Args:
            from_format: Source format
            to_format: Target format

        Returns:
            List of intermediate formats (including source and target)
        """
        if from_format == to_format:
            return [from_format]

        # Direct conversions
        direct_conversions = {
            # Tabular formats
            "csv": ["tsv", "excel", "json", "rds", "h5ad", "pickle"],
            "tsv": ["csv", "excel", "json"],
            "excel": ["csv", "tsv", "json"],
            "json": ["csv", "tsv", "excel"],
            "rds": ["csv"],
            "h5ad": ["csv"],
            "pickle": ["csv"],
            # Sequence formats
            "fastq.gz": ["fastq"],
            "fastq": ["fasta", "tsv"],
            "fasta": ["tsv"],
            # Interval formats
            "bigwig": ["bedgraph"],
            "bedgraph": ["bed", "bigwig"],
            "bed": ["csv", "bedgraph"],
            "bed3": ["csv", "bedgraph"],
            "bed6": ["csv", "bedgraph"],
            "bed12": ["csv", "bedgraph"],
            # Alignment formats
            "sam": ["bam", "tsv"],
            "bam": ["sam", "tsv"],
            # Variant formats
            "vcf": ["bcf", "tsv", "bed"],
            "vcf.gz": ["vcf"],
            "bcf": ["vcf"],
            # Annotation formats
            "gtf": ["gff3", "tsv", "bed"],
            "gff3": ["gtf", "tsv", "bed"],
        }

        # Check direct conversion
        if to_format in direct_conversions.get(from_format, []):
            return [from_format, to_format]

        # Try via CSV as intermediary (tabular formats)
        if "csv" in direct_conversions.get(
            from_format, []
        ) and to_format in direct_conversions.get("csv", []):
            return [from_format, "csv", to_format]

        # Try multi-step conversions for sequence formats
        # fastq.gz → fastq → fasta → tsv
        if from_format == "fastq.gz" and to_format == "fasta":
            return ["fastq.gz", "fastq", "fasta"]
        if from_format == "fastq.gz" and to_format == "tsv":
            return ["fastq.gz", "fastq", "tsv"]

        # Try multi-step conversions for interval formats
        # bigwig → bedgraph → bed → csv
        if from_format == "bigwig" and to_format == "bed":
            return ["bigwig", "bedgraph", "bed"]
        if from_format == "bigwig" and to_format == "csv":
            return ["bigwig", "bedgraph", "bed", "csv"]
        if from_format == "bedgraph" and to_format == "csv":
            return ["bedgraph", "bed", "csv"]

        raise ValueError(f"No conversion path from {from_format} to {to_format}")

    def estimate_conversion_time(
        self, file_path: str, from_format: str, to_format: str
    ) -> float:
        """
        Estimate conversion time in seconds.

        Args:
            file_path: Path to source file
            from_format: Source format
            to_format: Target format

        Returns:
            Estimated time in seconds
        """
        file_size_gb = Path(file_path).stat().st_size / (1024**3)

        conversion_path = self.get_conversion_path(from_format, to_format)
        total_time = 0

        for i in range(len(conversion_path) - 1):
            pair = (conversion_path[i], conversion_path[i + 1])
            # Get estimate, default to 10s/GB if not found
            time_per_gb = self.CONVERSION_TIME_ESTIMATES.get(pair, 10)
            # For symmetric conversions
            if pair not in self.CONVERSION_TIME_ESTIMATES:
                reverse_pair = (pair[1], pair[0])
                time_per_gb = self.CONVERSION_TIME_ESTIMATES.get(reverse_pair, 10)

            total_time += time_per_gb * file_size_gb

        # Minimum 1 second
        return max(total_time, 1.0)

    def convert(
        self,
        source_path: str,
        target_path: str,
        from_format: Optional[str] = None,
        to_format: Optional[str] = None,
        **kwargs,
    ) -> Dict:
        """
        Convert file from one format to another.

        Args:
            source_path: Path to source file
            target_path: Path to target file
            from_format: Source format (auto-detect if None)
            to_format: Target format (auto-detect if None)
            **kwargs: Additional conversion parameters

        Returns:
            Dict with conversion result info
        """
        start_time = datetime.now()

        # Auto-detect formats if not provided
        if from_format is None:
            from_format = self.detect_format(source_path)
            if from_format is None:
                raise ValueError(f"Cannot detect format for {source_path}")

        if to_format is None:
            to_format = self.detect_format(target_path)
            if to_format is None:
                raise ValueError(f"Cannot detect format for {target_path}")

        logger.info(f"Converting {source_path} from {from_format} to {to_format}")

        # Get conversion path
        conversion_path = self.get_conversion_path(from_format, to_format)

        # Perform conversions
        current_path = source_path
        for i in range(len(conversion_path) - 1):
            src_fmt = conversion_path[i]
            dst_fmt = conversion_path[i + 1]

            # Last step uses target_path, intermediate steps use temp files
            if i == len(conversion_path) - 2:
                next_path = target_path
            else:
                next_path = str(self.storage_path / f"temp_{i}.{dst_fmt}")

            self._convert_direct(current_path, next_path, src_fmt, dst_fmt, **kwargs)
            current_path = next_path

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        result = {
            "success": True,
            "source_path": source_path,
            "target_path": target_path,
            "from_format": from_format,
            "to_format": to_format,
            "conversion_path": conversion_path,
            "duration_seconds": duration,
            "target_size_bytes": Path(target_path).stat().st_size,
        }

        logger.info(f"Conversion completed in {duration:.2f}s")
        return result

    def _convert_direct(
        self,
        source_path: str,
        target_path: str,
        from_format: str,
        to_format: str,
        **kwargs,
    ):
        """
        Perform direct conversion between two formats.

        Args:
            source_path: Source file path
            target_path: Target file path
            from_format: Source format
            to_format: Target format
            **kwargs: Additional parameters
        """
        # Tabular format conversions
        if from_format in ["csv", "tsv", "excel", "json"] and to_format in [
            "csv",
            "tsv",
            "excel",
            "json",
        ]:
            self._convert_tabular(source_path, target_path, from_format, to_format)

        # CSV to R (RDS)
        elif from_format == "csv" and to_format == "rds":
            self._convert_csv_to_rds(source_path, target_path)

        # R (RDS) to CSV
        elif from_format == "rds" and to_format == "csv":
            self._convert_rds_to_csv(source_path, target_path)

        # CSV to Python (h5ad)
        elif from_format == "csv" and to_format == "h5ad":
            self._convert_csv_to_h5ad(source_path, target_path, **kwargs)

        # Python (h5ad) to CSV
        elif from_format == "h5ad" and to_format == "csv":
            self._convert_h5ad_to_csv(source_path, target_path)

        # CSV to pickle
        elif from_format == "csv" and to_format == "pickle":
            self._convert_csv_to_pickle(source_path, target_path)

        # Pickle to CSV
        elif from_format == "pickle" and to_format == "csv":
            self._convert_pickle_to_csv(source_path, target_path)

        # === Sequence format conversions ===
        elif from_format == "fastq.gz" and to_format == "fastq":
            self.sequence_converter.decompress_fastq_gz(source_path, target_path)

        elif from_format == "fastq" and to_format == "fasta":
            self.sequence_converter.convert_fastq_to_fasta(source_path, target_path)

        elif from_format == "fastq" and to_format == "tsv":
            include_quality = kwargs.get("include_quality", True)
            self.sequence_converter.convert_fastq_to_tsv(
                source_path, target_path, include_quality
            )

        elif from_format == "fasta" and to_format == "tsv":
            self.sequence_converter.convert_fasta_to_tsv(source_path, target_path)

        # === Interval format conversions ===
        elif from_format in ["bed", "bed3", "bed6", "bed12"] and to_format == "csv":
            bed_type = from_format.replace("bed", "") if from_format != "bed" else "3"
            self.interval_converter.convert_bed_to_csv(
                source_path, target_path, bed_type
            )

        elif from_format == "csv" and to_format in ["bed", "bed3", "bed6", "bed12"]:
            bed_type = to_format.replace("bed", "") if to_format != "bed" else "3"
            self.interval_converter.convert_csv_to_bed(
                source_path, target_path, bed_type
            )

        elif from_format == "bedgraph" and to_format == "bed":
            self.interval_converter.convert_bedgraph_to_bed(source_path, target_path)

        elif (
            from_format in ["bed", "bed3", "bed6", "bed12"] and to_format == "bedgraph"
        ):
            self.interval_converter.convert_bed_to_bedgraph(source_path, target_path)

        elif from_format == "bigwig" and to_format == "bedgraph":
            chrom = kwargs.get("chromosome")
            self.interval_converter.convert_bigwig_to_bedgraph(
                source_path, target_path, chrom
            )

        elif from_format == "bedgraph" and to_format == "bigwig":
            chrom_sizes = kwargs.get("chrom_sizes")
            if not chrom_sizes:
                raise ValueError(
                    "chrom_sizes required for bedGraph to BigWig conversion"
                )
            self.interval_converter.convert_bedgraph_to_bigwig(
                source_path, target_path, chrom_sizes
            )

        # === Alignment format conversions ===
        elif from_format == "sam" and to_format == "bam":
            sort = kwargs.get("sort", False)
            index = kwargs.get("index", False)
            self.alignment_converter.convert_sam_to_bam(
                source_path, target_path, sort, index
            )

        elif from_format == "bam" and to_format == "sam":
            include_header = kwargs.get("include_header", True)
            self.alignment_converter.convert_bam_to_sam(
                source_path, target_path, include_header
            )

        elif from_format == "sam" and to_format == "tsv":
            max_reads = kwargs.get("max_reads")
            include_tags = kwargs.get("include_tags", False)
            self.alignment_converter.convert_sam_to_tsv(
                source_path, target_path, max_reads, include_tags
            )

        elif from_format == "bam" and to_format == "tsv":
            max_reads = kwargs.get("max_reads")
            include_tags = kwargs.get("include_tags", False)
            self.alignment_converter.convert_bam_to_tsv(
                source_path, target_path, max_reads, include_tags
            )

        # === Variant format conversions ===
        elif from_format == "vcf" and to_format == "bcf":
            index = kwargs.get("index", False)
            self.variant_converter.convert_vcf_to_bcf(source_path, target_path, index)

        elif from_format == "bcf" and to_format == "vcf":
            compress = kwargs.get("compress", False)
            self.variant_converter.convert_bcf_to_vcf(
                source_path, target_path, compress
            )

        elif from_format == "vcf" and to_format == "tsv":
            max_variants = kwargs.get("max_variants")
            include_info = kwargs.get("include_info", True)
            include_samples = kwargs.get("include_samples", False)
            self.variant_converter.convert_vcf_to_tsv(
                source_path, target_path, max_variants, include_info, include_samples
            )

        elif from_format == "vcf" and to_format == "bed":
            extend = kwargs.get("extend", 0)
            variant_type = kwargs.get("variant_type")
            self.variant_converter.convert_vcf_to_bed(
                source_path, target_path, extend, variant_type
            )

        # === Annotation format conversions ===
        elif from_format == "gtf" and to_format == "gff3":
            self.annotation_converter.convert_gtf_to_gff3(source_path, target_path)

        elif from_format == "gff3" and to_format == "gtf":
            self.annotation_converter.convert_gff3_to_gtf(source_path, target_path)

        elif from_format == "gtf" and to_format == "tsv":
            max_features = kwargs.get("max_features")
            feature_types = kwargs.get("feature_types")
            self.annotation_converter.convert_gtf_to_tsv(
                source_path, target_path, max_features, feature_types
            )

        elif from_format == "gff3" and to_format == "tsv":
            max_features = kwargs.get("max_features")
            feature_types = kwargs.get("feature_types")
            self.annotation_converter.convert_gff3_to_tsv(
                source_path, target_path, max_features, feature_types
            )

        elif from_format == "gtf" and to_format == "bed":
            feature_type = kwargs.get("feature_type", "exon")
            name_field = kwargs.get("name_field", "gene_name")
            self.annotation_converter.convert_gtf_to_bed(
                source_path, target_path, feature_type, name_field
            )

        elif from_format == "gff3" and to_format == "bed":
            feature_type = kwargs.get("feature_type", "exon")
            name_field = kwargs.get("name_field", "Name")
            self.annotation_converter.convert_gff3_to_bed(
                source_path, target_path, feature_type, name_field
            )

        # === Expression matrix format conversions ===
        elif from_format == "mtx" and to_format == "csv":
            genes_file = kwargs.get("genes_file")
            barcodes_file = kwargs.get("barcodes_file")
            self.expression_converter.convert_mtx_to_csv(
                source_path, target_path, genes_file, barcodes_file
            )

        elif from_format == "csv" and to_format == "mtx":
            genes_file = kwargs.get("genes_file")
            barcodes_file = kwargs.get("barcodes_file")
            self.expression_converter.convert_csv_to_mtx(
                source_path, target_path, genes_file, barcodes_file
            )

        elif from_format == "10x_h5" and to_format == "csv":
            self.expression_converter.convert_10x_h5_to_csv(source_path, target_path)

        elif from_format == "10x_h5" and to_format == "h5ad":
            self.expression_converter.convert_10x_h5_to_h5ad(source_path, target_path)

        elif from_format == "h5ad" and to_format == "csv":
            layer = kwargs.get("layer")
            self.expression_converter.convert_h5ad_to_csv(
                source_path, target_path, layer
            )

        elif from_format == "csv" and to_format == "h5ad":
            self.expression_converter.convert_csv_to_h5ad(source_path, target_path)

        # === Python data format conversions ===
        elif from_format == "h5" and to_format == "csv":
            dataset_path = kwargs.get("dataset_path", "data")
            self.python_data_converter.convert_h5_to_csv(
                source_path, target_path, dataset_path
            )

        elif from_format == "csv" and to_format == "h5":
            dataset_path = kwargs.get("dataset_path", "data")
            self.python_data_converter.convert_csv_to_h5(
                source_path, target_path, dataset_path
            )

        elif from_format == "h5" and to_format == "npy":
            dataset_path = kwargs.get("dataset_path", "data")
            self.python_data_converter.convert_h5_to_npy(
                source_path, target_path, dataset_path
            )

        elif from_format == "npy" and to_format == "csv":
            self.python_data_converter.convert_npy_to_csv(source_path, target_path)

        elif from_format == "csv" and to_format == "npy":
            self.python_data_converter.convert_csv_to_npy(source_path, target_path)

        elif from_format == "npy" and to_format == "npz":
            array_name = kwargs.get("array_name", "data")
            self.python_data_converter.convert_npy_to_npz(
                source_path, target_path, array_name
            )

        elif from_format == "npz" and to_format == "npy":
            array_name = kwargs.get("array_name")
            self.python_data_converter.convert_npz_to_npy(
                source_path, target_path, array_name
            )

        elif from_format == "npz" and to_format == "csv":
            array_name = kwargs.get("array_name")
            self.python_data_converter.convert_npz_to_csv(
                source_path, target_path, array_name
            )

        # === R data format conversions ===
        elif from_format == "rdata" and to_format == "csv":
            object_name = kwargs.get("object_name")
            self.r_data_converter.convert_rdata_to_csv(
                source_path, target_path, object_name
            )

        elif from_format == "rdata" and to_format == "json":
            object_name = kwargs.get("object_name")
            self.r_data_converter.convert_rdata_to_json(
                source_path, target_path, object_name
            )

        # === Compression operations ===
        elif from_format == "gzip" and to_format == "decompress":
            self.compression_handler.decompress_gzip(source_path, target_path)

        elif from_format == "bgzip" and to_format == "decompress":
            self.compression_handler.decompress_bgzip(source_path, target_path)

        elif from_format == "zip" and to_format == "decompress":
            self.compression_handler.decompress_zip(source_path, target_path)

        elif from_format == "tar" and to_format == "decompress":
            self.compression_handler.decompress_tar(source_path, target_path)

        elif from_format == "tar.gz" and to_format == "decompress":
            self.compression_handler.decompress_tar(source_path, target_path)

        else:
            raise ValueError(f"Unsupported conversion: {from_format} -> {to_format}")

    def _convert_tabular(
        self, source_path: str, target_path: str, from_format: str, to_format: str
    ):
        """Convert between tabular formats using pandas."""
        # Read source
        if from_format == "csv":
            df = pd.read_csv(source_path)
        elif from_format == "tsv":
            df = pd.read_csv(source_path, sep="\t")
        elif from_format == "excel":
            df = pd.read_excel(source_path)
        elif from_format == "json":
            df = pd.read_json(source_path)
        else:
            raise ValueError(f"Unknown source format: {from_format}")

        # Write target
        if to_format == "csv":
            df.to_csv(target_path, index=False)
        elif to_format == "tsv":
            df.to_csv(target_path, sep="\t", index=False)
        elif to_format == "excel":
            df.to_excel(target_path, index=False)
        elif to_format == "json":
            df.to_json(target_path, orient="records", indent=2)
        else:
            raise ValueError(f"Unknown target format: {to_format}")

    def _convert_csv_to_rds(self, source_path: str, target_path: str):
        """Convert CSV to R RDS format using Rscript."""
        r_script = f"""
        data <- read.csv("{source_path}", stringsAsFactors = FALSE)
        saveRDS(data, "{target_path}")
        """

        # Write temporary R script
        script_path = self.storage_path / "convert_to_rds.R"
        script_path.write_text(r_script)

        # Execute R script
        result = subprocess.run(
            ["Rscript", str(script_path)], capture_output=True, text=True
        )

        if result.returncode != 0:
            raise RuntimeError(f"R conversion failed: {result.stderr}")

    def _convert_rds_to_csv(self, source_path: str, target_path: str):
        """Convert R RDS to CSV using Rscript."""
        r_script = f"""
        data <- readRDS("{source_path}")
        write.csv(data, "{target_path}", row.names = FALSE)
        """

        # Write temporary R script
        script_path = self.storage_path / "convert_from_rds.R"
        script_path.write_text(r_script)

        # Execute R script
        result = subprocess.run(
            ["Rscript", str(script_path)], capture_output=True, text=True
        )

        if result.returncode != 0:
            raise RuntimeError(f"R conversion failed: {result.stderr}")

    def _convert_csv_to_h5ad(self, source_path: str, target_path: str, **kwargs):
        """Convert CSV to h5ad (AnnData) format."""
        try:
            import anndata as ad
        except ImportError:
            raise ImportError("anndata package required for h5ad conversion")

        # Read CSV
        df = pd.read_csv(source_path, index_col=0)

        # Create AnnData object
        adata = ad.AnnData(X=df.values)
        adata.obs_names = df.index.astype(str)
        adata.var_names = df.columns.astype(str)

        # Add metadata if provided
        if "obs" in kwargs:
            adata.obs = kwargs["obs"]
        if "var" in kwargs:
            adata.var = kwargs["var"]

        # Write h5ad
        adata.write_h5ad(target_path)

    def _convert_h5ad_to_csv(self, source_path: str, target_path: str):
        """Convert h5ad (AnnData) to CSV."""
        try:
            import anndata as ad
        except ImportError:
            raise ImportError("anndata package required for h5ad conversion")

        # Read h5ad
        adata = ad.read_h5ad(source_path)

        # Convert to DataFrame
        df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

        # Write CSV
        df.to_csv(target_path)

    def _convert_csv_to_pickle(self, source_path: str, target_path: str):
        """Convert CSV to Python pickle format."""
        df = pd.read_csv(source_path)
        df.to_pickle(target_path)

    def _convert_pickle_to_csv(self, source_path: str, target_path: str):
        """Convert Python pickle to CSV."""
        df = pd.read_pickle(source_path)
        df.to_csv(target_path, index=False)


# Singleton instance
_converter_instance = None


def get_format_converter() -> FormatConverter:
    """Get singleton format converter instance."""
    global _converter_instance
    if _converter_instance is None:
        _converter_instance = FormatConverter()
    return _converter_instance
