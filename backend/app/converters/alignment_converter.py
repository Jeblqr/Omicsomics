"""
Alignment format converter for SAM/BAM files.

Supports conversions between:
- SAM (Sequence Alignment/Map format)
- BAM (Binary SAM)
- TSV (Tabular export)

Dependencies:
- pysam: SAM/BAM file handling
"""

from pathlib import Path
from typing import Dict, Optional, List
import subprocess
import tempfile


class AlignmentConverter:
    """
    Converter for alignment file formats (SAM, BAM).

    Singleton pattern for efficiency.
    """

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize alignment converter."""
        if not hasattr(self, "_initialized"):
            self._initialized = True
            self._check_dependencies()

    def _check_dependencies(self):
        """Check if required tools are available."""
        try:
            import pysam

            self.pysam = pysam
            self.has_pysam = True
        except ImportError:
            self.has_pysam = False

    # ============================================================
    # SAM/BAM Conversions
    # ============================================================

    def convert_sam_to_bam(
        self, source: str, target: str, sort: bool = False, index: bool = False
    ):
        """
        Convert SAM to BAM format.

        Args:
            source: Input SAM file path
            target: Output BAM file path
            sort: Sort BAM by coordinates
            index: Create BAM index (.bai)
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for SAM/BAM conversion")

        source_path = Path(source)
        target_path = Path(target)

        # Read SAM and write BAM
        with self.pysam.AlignmentFile(source_path, "r") as sam_in:
            with self.pysam.AlignmentFile(
                target_path, "wb", header=sam_in.header
            ) as bam_out:
                for read in sam_in:
                    bam_out.write(read)

        # Sort BAM if requested
        if sort:
            sorted_path = target_path.with_suffix(".sorted.bam")
            self.pysam.sort("-o", str(sorted_path), str(target_path))
            sorted_path.replace(target_path)

        # Create index if requested
        if index:
            self.pysam.index(str(target_path))

    def convert_bam_to_sam(self, source: str, target: str, include_header: bool = True):
        """
        Convert BAM to SAM format.

        Args:
            source: Input BAM file path
            target: Output SAM file path
            include_header: Include SAM header
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for SAM/BAM conversion")

        source_path = Path(source)
        target_path = Path(target)

        # Read BAM and write SAM
        mode = "r" if include_header else "r"
        with self.pysam.AlignmentFile(source_path, "rb") as bam_in:
            with self.pysam.AlignmentFile(
                target_path, "w", header=bam_in.header
            ) as sam_out:
                for read in bam_in:
                    sam_out.write(read)

    def convert_sam_to_tsv(
        self,
        source: str,
        target: str,
        max_reads: Optional[int] = None,
        include_tags: bool = False,
    ):
        """
        Convert SAM to TSV format for analysis.

        Args:
            source: Input SAM file path
            target: Output TSV file path
            max_reads: Maximum reads to export (None for all)
            include_tags: Include optional tags
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for SAM/BAM conversion")

        source_path = Path(source)
        target_path = Path(target)

        with self.pysam.AlignmentFile(source_path, "r") as sam_in:
            with open(target_path, "w") as tsv_out:
                # Write header
                header = [
                    "qname",
                    "flag",
                    "rname",
                    "pos",
                    "mapq",
                    "cigar",
                    "rnext",
                    "pnext",
                    "tlen",
                    "seq",
                    "qual",
                ]
                if include_tags:
                    header.append("tags")

                tsv_out.write("\t".join(header) + "\n")

                # Write alignments
                count = 0
                for read in sam_in:
                    if max_reads and count >= max_reads:
                        break

                    row = [
                        read.query_name,
                        str(read.flag),
                        read.reference_name if read.reference_name else "*",
                        (
                            str(read.reference_start + 1)
                            if read.reference_start is not None
                            else "0"
                        ),
                        str(read.mapping_quality),
                        read.cigarstring if read.cigarstring else "*",
                        read.next_reference_name if read.next_reference_name else "*",
                        (
                            str(read.next_reference_start + 1)
                            if read.next_reference_start is not None
                            else "0"
                        ),
                        str(read.template_length),
                        read.query_sequence if read.query_sequence else "*",
                        (
                            "".join([chr(q + 33) for q in read.query_qualities])
                            if read.query_qualities
                            else "*"
                        ),
                    ]

                    if include_tags:
                        tags = ",".join(
                            [
                                f"{tag}:{typ}:{val}"
                                for tag, val in read.get_tags()
                                for typ in [read.get_tag_type(tag)]
                            ]
                        )
                        row.append(tags)

                    tsv_out.write("\t".join(row) + "\n")
                    count += 1

    def convert_bam_to_tsv(
        self,
        source: str,
        target: str,
        max_reads: Optional[int] = None,
        include_tags: bool = False,
    ):
        """
        Convert BAM to TSV format for analysis.

        Args:
            source: Input BAM file path
            target: Output TSV file path
            max_reads: Maximum reads to export (None for all)
            include_tags: Include optional tags
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for SAM/BAM conversion")

        source_path = Path(source)
        target_path = Path(target)

        with self.pysam.AlignmentFile(source_path, "rb") as bam_in:
            with open(target_path, "w") as tsv_out:
                # Write header
                header = [
                    "qname",
                    "flag",
                    "rname",
                    "pos",
                    "mapq",
                    "cigar",
                    "rnext",
                    "pnext",
                    "tlen",
                    "seq",
                    "qual",
                ]
                if include_tags:
                    header.append("tags")

                tsv_out.write("\t".join(header) + "\n")

                # Write alignments
                count = 0
                for read in bam_in:
                    if max_reads and count >= max_reads:
                        break

                    row = [
                        read.query_name,
                        str(read.flag),
                        read.reference_name if read.reference_name else "*",
                        (
                            str(read.reference_start + 1)
                            if read.reference_start is not None
                            else "0"
                        ),
                        str(read.mapping_quality),
                        read.cigarstring if read.cigarstring else "*",
                        read.next_reference_name if read.next_reference_name else "*",
                        (
                            str(read.next_reference_start + 1)
                            if read.next_reference_start is not None
                            else "0"
                        ),
                        str(read.template_length),
                        read.query_sequence if read.query_sequence else "*",
                        (
                            "".join([chr(q + 33) for q in read.query_qualities])
                            if read.query_qualities
                            else "*"
                        ),
                    ]

                    if include_tags:
                        tags = ",".join(
                            [
                                f"{tag}:{typ}:{val}"
                                for tag, val in read.get_tags()
                                for typ in [read.get_tag_type(tag)]
                            ]
                        )
                        row.append(tags)

                    tsv_out.write("\t".join(row) + "\n")
                    count += 1

    # ============================================================
    # BAM Operations
    # ============================================================

    def sort_bam(self, source: str, target: str, by_name: bool = False):
        """
        Sort BAM file by coordinates or read name.

        Args:
            source: Input BAM file path
            target: Output sorted BAM file path
            by_name: Sort by read name instead of coordinates
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for BAM sorting")

        if by_name:
            self.pysam.sort("-n", "-o", str(target), str(source))
        else:
            self.pysam.sort("-o", str(target), str(source))

    def index_bam(self, bam_file: str):
        """
        Create BAM index file (.bai).

        Args:
            bam_file: BAM file to index
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for BAM indexing")

        self.pysam.index(str(bam_file))

    def filter_bam(
        self,
        source: str,
        target: str,
        min_mapq: int = 0,
        exclude_flags: int = 0,
        require_flags: int = 0,
    ):
        """
        Filter BAM file by mapping quality and flags.

        Args:
            source: Input BAM file path
            target: Output filtered BAM file path
            min_mapq: Minimum mapping quality
            exclude_flags: Flag bits to exclude (e.g., 0x4 for unmapped)
            require_flags: Flag bits that must be set
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for BAM filtering")

        with self.pysam.AlignmentFile(source, "rb") as bam_in:
            with self.pysam.AlignmentFile(
                target, "wb", header=bam_in.header
            ) as bam_out:
                for read in bam_in:
                    # Filter by mapping quality
                    if read.mapping_quality < min_mapq:
                        continue

                    # Filter by flags
                    if exclude_flags and (read.flag & exclude_flags):
                        continue

                    if require_flags and not (read.flag & require_flags):
                        continue

                    bam_out.write(read)

    # ============================================================
    # Statistics
    # ============================================================

    def get_alignment_stats(self, file_path: str, file_format: str = "bam") -> Dict:
        """
        Get alignment statistics.

        Args:
            file_path: SAM or BAM file path
            file_format: 'sam' or 'bam'

        Returns:
            Dictionary with statistics
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for alignment statistics")

        mode = "rb" if file_format == "bam" else "r"

        stats = {
            "total_reads": 0,
            "mapped_reads": 0,
            "unmapped_reads": 0,
            "paired_reads": 0,
            "proper_pairs": 0,
            "duplicates": 0,
            "secondary": 0,
            "supplementary": 0,
            "total_bases": 0,
            "mapq_distribution": {},
            "references": set(),
        }

        with self.pysam.AlignmentFile(file_path, mode) as aln_file:
            for read in aln_file:
                stats["total_reads"] += 1

                # Mapping status
                if read.is_unmapped:
                    stats["unmapped_reads"] += 1
                else:
                    stats["mapped_reads"] += 1
                    stats["references"].add(read.reference_name)

                # Pairing
                if read.is_paired:
                    stats["paired_reads"] += 1
                    if read.is_proper_pair:
                        stats["proper_pairs"] += 1

                # Flags
                if read.is_duplicate:
                    stats["duplicates"] += 1
                if read.is_secondary:
                    stats["secondary"] += 1
                if read.is_supplementary:
                    stats["supplementary"] += 1

                # Sequence
                if read.query_sequence:
                    stats["total_bases"] += len(read.query_sequence)

                # Mapping quality
                mapq = read.mapping_quality
                stats["mapq_distribution"][mapq] = (
                    stats["mapq_distribution"].get(mapq, 0) + 1
                )

        # Convert reference set to list
        stats["references"] = sorted(list(stats["references"]))
        stats["num_references"] = len(stats["references"])

        # Calculate percentages
        if stats["total_reads"] > 0:
            stats["mapped_percent"] = (
                stats["mapped_reads"] / stats["total_reads"]
            ) * 100
            stats["unmapped_percent"] = (
                stats["unmapped_reads"] / stats["total_reads"]
            ) * 100

        return stats


def get_alignment_converter() -> AlignmentConverter:
    """Get singleton instance of AlignmentConverter."""
    return AlignmentConverter()
