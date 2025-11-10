"""
Sequence format converter for genomic data.

Supports conversions between:
- FASTA, FASTQ, FASTQ.gz
- Sequence list (TSV/CSV)
"""

import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime

logger = logging.getLogger(__name__)


class SequenceConverter:
    """
    Handle sequence format conversions (FASTA, FASTQ, etc.).
    """

    def __init__(self):
        try:
            from Bio import SeqIO

            self.SeqIO = SeqIO
            self.biopython_available = True
        except ImportError:
            logger.warning(
                "BioPython not available. Sequence conversion will be limited."
            )
            self.biopython_available = False

    def convert_fastq_to_fasta(self, source_path: str, target_path: str) -> Dict:
        """
        Convert FASTQ to FASTA (removes quality scores).

        Args:
            source_path: Path to input FASTQ file
            target_path: Path to output FASTA file

        Returns:
            Conversion result info
        """
        if not self.biopython_available:
            raise RuntimeError("BioPython is required for sequence conversion")

        start_time = datetime.now()

        # Check if input is gzipped
        is_gzipped = source_path.endswith(".gz")

        if is_gzipped:
            input_handle = gzip.open(source_path, "rt")
        else:
            input_handle = open(source_path, "r")

        try:
            # Convert
            count = self.SeqIO.convert(input_handle, "fastq", target_path, "fasta")

            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()

            return {
                "success": True,
                "sequences_converted": count,
                "duration_seconds": duration,
                "source_format": "fastq",
                "target_format": "fasta",
            }
        finally:
            input_handle.close()

    def convert_fasta_to_tsv(self, source_path: str, target_path: str) -> Dict:
        """
        Convert FASTA to TSV (sequence list).

        Output format: ID, Sequence, Length, Description

        Args:
            source_path: Path to input FASTA file
            target_path: Path to output TSV file

        Returns:
            Conversion result info
        """
        if not self.biopython_available:
            raise RuntimeError("BioPython is required for sequence conversion")

        start_time = datetime.now()
        count = 0

        with open(target_path, "w") as out_f:
            # Write header
            out_f.write("ID\tSequence\tLength\tDescription\n")

            # Parse and write sequences
            for record in self.SeqIO.parse(source_path, "fasta"):
                out_f.write(
                    f"{record.id}\t{str(record.seq)}\t{len(record.seq)}\t{record.description}\n"
                )
                count += 1

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        return {
            "success": True,
            "sequences_converted": count,
            "duration_seconds": duration,
            "source_format": "fasta",
            "target_format": "tsv",
        }

    def convert_fastq_to_tsv(
        self, source_path: str, target_path: str, include_quality: bool = True
    ) -> Dict:
        """
        Convert FASTQ to TSV (sequence list with quality).

        Output format: ID, Sequence, Length, Quality, Avg_Quality, Description

        Args:
            source_path: Path to input FASTQ file
            target_path: Path to output TSV file
            include_quality: Whether to include quality string

        Returns:
            Conversion result info
        """
        if not self.biopython_available:
            raise RuntimeError("BioPython is required for sequence conversion")

        start_time = datetime.now()
        count = 0

        # Check if input is gzipped
        is_gzipped = source_path.endswith(".gz")

        if is_gzipped:
            input_handle = gzip.open(source_path, "rt")
        else:
            input_handle = open(source_path, "r")

        try:
            with open(target_path, "w") as out_f:
                # Write header
                if include_quality:
                    out_f.write(
                        "ID\tSequence\tLength\tQuality\tAvg_Quality\tDescription\n"
                    )
                else:
                    out_f.write("ID\tSequence\tLength\tAvg_Quality\tDescription\n")

                # Parse and write sequences
                for record in self.SeqIO.parse(input_handle, "fastq"):
                    seq_str = str(record.seq)
                    qual_scores = record.letter_annotations.get("phred_quality", [])
                    avg_qual = sum(qual_scores) / len(qual_scores) if qual_scores else 0

                    if include_quality:
                        qual_str = "".join([chr(q + 33) for q in qual_scores])
                        out_f.write(
                            f"{record.id}\t{seq_str}\t{len(seq_str)}\t{qual_str}\t{avg_qual:.2f}\t{record.description}\n"
                        )
                    else:
                        out_f.write(
                            f"{record.id}\t{seq_str}\t{len(seq_str)}\t{avg_qual:.2f}\t{record.description}\n"
                        )

                    count += 1
        finally:
            input_handle.close()

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        return {
            "success": True,
            "sequences_converted": count,
            "duration_seconds": duration,
            "source_format": "fastq",
            "target_format": "tsv",
        }

    def decompress_fastq_gz(self, source_path: str, target_path: str) -> Dict:
        """
        Decompress FASTQ.gz to FASTQ.

        Args:
            source_path: Path to input FASTQ.gz file
            target_path: Path to output FASTQ file

        Returns:
            Conversion result info
        """
        start_time = datetime.now()

        with gzip.open(source_path, "rb") as f_in:
            with open(target_path, "wb") as f_out:
                # Copy in chunks for large files
                chunk_size = 1024 * 1024  # 1MB chunks
                while True:
                    chunk = f_in.read(chunk_size)
                    if not chunk:
                        break
                    f_out.write(chunk)

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        # Get file sizes
        source_size = Path(source_path).stat().st_size
        target_size = Path(target_path).stat().st_size

        return {
            "success": True,
            "source_size_bytes": source_size,
            "target_size_bytes": target_size,
            "compression_ratio": source_size / target_size if target_size > 0 else 0,
            "duration_seconds": duration,
            "source_format": "fastq.gz",
            "target_format": "fastq",
        }

    def compress_fastq(self, source_path: str, target_path: str) -> Dict:
        """
        Compress FASTQ to FASTQ.gz.

        Args:
            source_path: Path to input FASTQ file
            target_path: Path to output FASTQ.gz file

        Returns:
            Conversion result info
        """
        start_time = datetime.now()

        with open(source_path, "rb") as f_in:
            with gzip.open(target_path, "wb", compresslevel=6) as f_out:
                # Copy in chunks
                chunk_size = 1024 * 1024  # 1MB chunks
                while True:
                    chunk = f_in.read(chunk_size)
                    if not chunk:
                        break
                    f_out.write(chunk)

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        # Get file sizes
        source_size = Path(source_path).stat().st_size
        target_size = Path(target_path).stat().st_size

        return {
            "success": True,
            "source_size_bytes": source_size,
            "target_size_bytes": target_size,
            "compression_ratio": target_size / source_size if source_size > 0 else 0,
            "duration_seconds": duration,
            "source_format": "fastq",
            "target_format": "fastq.gz",
        }

    def get_sequence_stats(self, file_path: str, format_type: str = "fasta") -> Dict:
        """
        Get statistics about a sequence file.

        Args:
            file_path: Path to sequence file
            format_type: File format ('fasta' or 'fastq')

        Returns:
            Statistics dictionary
        """
        if not self.biopython_available:
            raise RuntimeError("BioPython is required for sequence analysis")

        # Check if gzipped
        is_gzipped = file_path.endswith(".gz")

        if is_gzipped:
            input_handle = gzip.open(file_path, "rt")
        else:
            input_handle = open(file_path, "r")

        try:
            sequences = list(self.SeqIO.parse(input_handle, format_type))

            if not sequences:
                return {
                    "num_sequences": 0,
                    "total_bases": 0,
                    "avg_length": 0,
                    "min_length": 0,
                    "max_length": 0,
                }

            lengths = [len(seq.seq) for seq in sequences]

            stats = {
                "num_sequences": len(sequences),
                "total_bases": sum(lengths),
                "avg_length": sum(lengths) / len(lengths),
                "min_length": min(lengths),
                "max_length": max(lengths),
                "format": format_type,
                "compressed": is_gzipped,
            }

            # Add quality stats for FASTQ
            if format_type == "fastq":
                qual_scores = []
                for seq in sequences:
                    quals = seq.letter_annotations.get("phred_quality", [])
                    if quals:
                        qual_scores.append(sum(quals) / len(quals))

                if qual_scores:
                    stats["avg_quality"] = sum(qual_scores) / len(qual_scores)
                    stats["min_avg_quality"] = min(qual_scores)
                    stats["max_avg_quality"] = max(qual_scores)

            return stats
        finally:
            input_handle.close()


# Singleton instance
_sequence_converter = None


def get_sequence_converter() -> SequenceConverter:
    """Get singleton sequence converter instance."""
    global _sequence_converter
    if _sequence_converter is None:
        _sequence_converter = SequenceConverter()
    return _sequence_converter
