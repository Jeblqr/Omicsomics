"""
Tests for bioinformatics format converters (Phase 1 Week 1).

Tests sequence formats (FASTA, FASTQ, FASTQ.gz) and
interval formats (BED, bedGraph, BigWig).
"""

import pytest
from pathlib import Path
import gzip
import tempfile
import shutil

from app.converters.sequence_converter import get_sequence_converter
from app.converters.interval_converter import get_interval_converter
from app.converters.format_converter import FormatConverter


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    tmp = tempfile.mkdtemp()
    yield Path(tmp)
    shutil.rmtree(tmp)


@pytest.fixture
def sample_fasta(temp_dir):
    """Create sample FASTA file."""
    fasta_file = temp_dir / "sample.fasta"
    content = """>seq1 Sample sequence 1
ATCGATCGATCG
>seq2 Sample sequence 2
GCTAGCTAGCTA
>seq3 Sample sequence 3
TTAATTAATTAA
"""
    fasta_file.write_text(content)
    return fasta_file


@pytest.fixture
def sample_fastq(temp_dir):
    """Create sample FASTQ file."""
    fastq_file = temp_dir / "sample.fastq"
    content = """@seq1 Sample read 1
ATCGATCGATCG
+
IIIIIIIIIIII
@seq2 Sample read 2
GCTAGCTAGCTA
+
HHHHHHHHHHHH
@seq3 Sample read 3
TTAATTAATTAA
+
JJJJJJJJJJJJ
"""
    fastq_file.write_text(content)
    return fastq_file


@pytest.fixture
def sample_fastq_gz(temp_dir):
    """Create sample gzipped FASTQ file."""
    fastq_gz_file = temp_dir / "sample.fastq.gz"
    content = """@seq1 Sample read 1
ATCGATCGATCG
+
IIIIIIIIIIII
@seq2 Sample read 2
GCTAGCTAGCTA
+
HHHHHHHHHHHH
"""
    with gzip.open(fastq_gz_file, "wt") as f:
        f.write(content)
    return fastq_gz_file


@pytest.fixture
def sample_bed3(temp_dir):
    """Create sample BED3 file."""
    bed_file = temp_dir / "sample.bed"
    content = """chr1\t100\t200
chr1\t300\t400
chr2\t500\t600
"""
    bed_file.write_text(content)
    return bed_file


@pytest.fixture
def sample_bed6(temp_dir):
    """Create sample BED6 file."""
    bed_file = temp_dir / "sample_bed6.bed"
    content = """chr1\t100\t200\tfeature1\t100\t+
chr1\t300\t400\tfeature2\t200\t-
chr2\t500\t600\tfeature3\t150\t+
"""
    bed_file.write_text(content)
    return bed_file


@pytest.fixture
def sample_bedgraph(temp_dir):
    """Create sample bedGraph file."""
    bg_file = temp_dir / "sample.bedgraph"
    content = """chr1\t100\t200\t1.5
chr1\t200\t300\t2.3
chr2\t400\t500\t0.8
"""
    bg_file.write_text(content)
    return bg_file


# ============================================================
# Sequence Converter Tests
# ============================================================


class TestSequenceConverter:
    """Test SequenceConverter class."""

    def test_fastq_to_fasta(self, temp_dir, sample_fastq):
        """Test FASTQ to FASTA conversion."""
        converter = get_sequence_converter()
        output = temp_dir / "output.fasta"

        converter.convert_fastq_to_fasta(str(sample_fastq), str(output))

        assert output.exists()
        content = output.read_text()
        assert ">seq1" in content
        assert "ATCGATCGATCG" in content
        assert "+" not in content  # Quality lines removed
        assert "IIII" not in content  # Quality scores removed

    def test_fasta_to_tsv(self, temp_dir, sample_fasta):
        """Test FASTA to TSV conversion."""
        converter = get_sequence_converter()
        output = temp_dir / "output.tsv"

        converter.convert_fasta_to_tsv(str(sample_fasta), str(output))

        assert output.exists()
        lines = output.read_text().strip().split("\n")
        assert len(lines) == 4  # Header + 3 sequences
        assert lines[0] == "id\tdescription\tsequence\tlength"
        assert "seq1" in lines[1]
        assert "ATCGATCGATCG" in lines[1]

    def test_fastq_to_tsv_with_quality(self, temp_dir, sample_fastq):
        """Test FASTQ to TSV with quality scores."""
        converter = get_sequence_converter()
        output = temp_dir / "output.tsv"

        converter.convert_fastq_to_tsv(
            str(sample_fastq), str(output), include_quality=True
        )

        assert output.exists()
        lines = output.read_text().strip().split("\n")
        assert len(lines) == 4  # Header + 3 sequences
        assert "quality" in lines[0]
        assert "IIII" in lines[1]

    def test_fastq_to_tsv_without_quality(self, temp_dir, sample_fastq):
        """Test FASTQ to TSV without quality scores."""
        converter = get_sequence_converter()
        output = temp_dir / "output.tsv"

        converter.convert_fastq_to_tsv(
            str(sample_fastq), str(output), include_quality=False
        )

        assert output.exists()
        lines = output.read_text().strip().split("\n")
        assert "quality" not in lines[0]

    def test_decompress_fastq_gz(self, temp_dir, sample_fastq_gz):
        """Test FASTQ.gz decompression."""
        converter = get_sequence_converter()
        output = temp_dir / "output.fastq"

        converter.decompress_fastq_gz(str(sample_fastq_gz), str(output))

        assert output.exists()
        content = output.read_text()
        assert "@seq1" in content
        assert "ATCGATCGATCG" in content

    def test_compress_fastq(self, temp_dir, sample_fastq):
        """Test FASTQ compression."""
        converter = get_sequence_converter()
        output = temp_dir / "output.fastq.gz"

        converter.compress_fastq(str(sample_fastq), str(output))

        assert output.exists()
        with gzip.open(output, "rt") as f:
            content = f.read()
            assert "@seq1" in content

    def test_get_sequence_stats_fasta(self, temp_dir, sample_fasta):
        """Test sequence statistics for FASTA."""
        converter = get_sequence_converter()

        stats = converter.get_sequence_stats(str(sample_fasta), "fasta")

        assert stats["count"] == 3
        assert stats["total_length"] == 36
        assert stats["avg_length"] == 12.0
        assert stats["min_length"] == 12
        assert stats["max_length"] == 12

    def test_get_sequence_stats_fastq(self, temp_dir, sample_fastq):
        """Test sequence statistics for FASTQ."""
        converter = get_sequence_converter()

        stats = converter.get_sequence_stats(str(sample_fastq), "fastq")

        assert stats["count"] == 3
        assert "avg_quality" in stats
        assert stats["avg_quality"] > 30  # High quality scores


# ============================================================
# Interval Converter Tests
# ============================================================


class TestIntervalConverter:
    """Test IntervalConverter class."""

    def test_bed3_to_csv(self, temp_dir, sample_bed3):
        """Test BED3 to CSV conversion."""
        converter = get_interval_converter()
        output = temp_dir / "output.csv"

        converter.convert_bed_to_csv(str(sample_bed3), str(output), bed_type="3")

        assert output.exists()
        content = output.read_text()
        assert "chrom,start,end" in content
        assert "chr1,100,200" in content

    def test_bed6_to_csv(self, temp_dir, sample_bed6):
        """Test BED6 to CSV conversion."""
        converter = get_interval_converter()
        output = temp_dir / "output.csv"

        converter.convert_bed_to_csv(str(sample_bed6), str(output), bed_type="6")

        assert output.exists()
        content = output.read_text()
        assert "chrom,start,end,name,score,strand" in content
        assert "feature1" in content

    def test_csv_to_bed3(self, temp_dir):
        """Test CSV to BED3 conversion."""
        converter = get_interval_converter()

        # Create CSV file
        csv_file = temp_dir / "input.csv"
        csv_file.write_text("chrom,start,end\nchr1,100,200\nchr2,300,400\n")

        output = temp_dir / "output.bed"
        converter.convert_csv_to_bed(str(csv_file), str(output), bed_type="3")

        assert output.exists()
        content = output.read_text()
        assert "chr1\t100\t200" in content

    def test_bedgraph_to_bed(self, temp_dir, sample_bedgraph):
        """Test bedGraph to BED conversion."""
        converter = get_interval_converter()
        output = temp_dir / "output.bed"

        converter.convert_bedgraph_to_bed(str(sample_bedgraph), str(output))

        assert output.exists()
        lines = output.read_text().strip().split("\n")
        assert len(lines) == 3
        assert "chr1\t100\t200\t1.5" in lines[0]

    def test_bed_to_bedgraph(self, temp_dir):
        """Test BED to bedGraph conversion."""
        converter = get_interval_converter()

        # Create BED4 file with scores
        bed_file = temp_dir / "input.bed"
        bed_file.write_text("chr1\t100\t200\t1.5\nchr1\t200\t300\t2.3\n")

        output = temp_dir / "output.bedgraph"
        converter.convert_bed_to_bedgraph(str(bed_file), str(output))

        assert output.exists()
        content = output.read_text()
        assert "chr1\t100\t200\t1.5" in content

    def test_get_bed_stats(self, temp_dir, sample_bed3):
        """Test BED statistics calculation."""
        converter = get_interval_converter()

        stats = converter.get_bed_stats(str(sample_bed3))

        assert stats["count"] == 3
        assert stats["total_length"] == 300
        assert stats["avg_length"] == 100.0
        assert "chr1" in stats["chromosomes"]
        assert "chr2" in stats["chromosomes"]


# ============================================================
# FormatConverter Integration Tests
# ============================================================


class TestFormatConverterIntegration:
    """Test integration of new converters into FormatConverter."""

    def test_detect_sequence_formats(self, temp_dir, sample_fasta, sample_fastq):
        """Test format detection for sequence files."""
        converter = FormatConverter()

        assert converter.detect_format(str(sample_fasta)) == "fasta"
        assert converter.detect_format(str(sample_fastq)) == "fastq"

    def test_detect_interval_formats(self, temp_dir, sample_bed3, sample_bedgraph):
        """Test format detection for interval files."""
        converter = FormatConverter()

        assert converter.detect_format(str(sample_bed3)) == "bed"
        assert converter.detect_format(str(sample_bedgraph)) == "bedgraph"

    def test_conversion_path_sequence(self):
        """Test conversion path calculation for sequence formats."""
        converter = FormatConverter()

        # Direct conversions
        path = converter.get_conversion_path("fastq", "fasta")
        assert path == ["fastq", "fasta"]

        # Multi-step conversions
        path = converter.get_conversion_path("fastq.gz", "fasta")
        assert path == ["fastq.gz", "fastq", "fasta"]

        path = converter.get_conversion_path("fastq.gz", "tsv")
        assert path == ["fastq.gz", "fastq", "tsv"]

    def test_conversion_path_intervals(self):
        """Test conversion path calculation for interval formats."""
        converter = FormatConverter()

        # Direct conversions
        path = converter.get_conversion_path("bed", "csv")
        assert path == ["bed", "csv"]

        path = converter.get_conversion_path("bedgraph", "bed")
        assert path == ["bedgraph", "bed"]

        # Multi-step conversions
        path = converter.get_conversion_path("bigwig", "bed")
        assert path == ["bigwig", "bedgraph", "bed"]

        path = converter.get_conversion_path("bigwig", "csv")
        assert path == ["bigwig", "bedgraph", "bed", "csv"]

    def test_estimate_conversion_time(self, temp_dir, sample_fastq):
        """Test conversion time estimation."""
        converter = FormatConverter()

        # Estimate time for FASTQ to FASTA
        time_estimate = converter.estimate_conversion_time(
            str(sample_fastq), "fastq", "fasta"
        )
        assert time_estimate > 0

    def test_full_conversion_fastq_to_fasta(self, temp_dir, sample_fastq):
        """Test full conversion pipeline: FASTQ to FASTA."""
        converter = FormatConverter()
        output = temp_dir / "output.fasta"

        result = converter.convert(
            str(sample_fastq), str(output), from_format="fastq", to_format="fasta"
        )

        assert output.exists()
        assert result["status"] == "completed"

    def test_full_conversion_bed_to_csv(self, temp_dir, sample_bed3):
        """Test full conversion pipeline: BED to CSV."""
        converter = FormatConverter()
        output = temp_dir / "output.csv"

        result = converter.convert(
            str(sample_bed3), str(output), from_format="bed", to_format="csv"
        )

        assert output.exists()
        assert result["status"] == "completed"


# ============================================================
# Error Handling Tests
# ============================================================


class TestErrorHandling:
    """Test error handling in converters."""

    def test_invalid_fasta_file(self, temp_dir):
        """Test handling of invalid FASTA file."""
        converter = get_sequence_converter()

        invalid_file = temp_dir / "invalid.fasta"
        invalid_file.write_text("This is not a valid FASTA file\nwith random content")

        output = temp_dir / "output.tsv"

        # Should handle gracefully or raise appropriate error
        with pytest.raises(Exception):
            converter.convert_fasta_to_tsv(str(invalid_file), str(output))

    def test_invalid_bed_file(self, temp_dir):
        """Test handling of invalid BED file."""
        converter = get_interval_converter()

        invalid_file = temp_dir / "invalid.bed"
        invalid_file.write_text("chr1\t100\tNOT_A_NUMBER\n")

        output = temp_dir / "output.csv"

        # Should handle gracefully or raise appropriate error
        with pytest.raises(Exception):
            converter.convert_bed_to_csv(str(invalid_file), str(output), bed_type="3")

    def test_unsupported_conversion(self):
        """Test unsupported conversion path."""
        converter = FormatConverter()

        with pytest.raises(ValueError, match="No conversion path"):
            converter.get_conversion_path("fasta", "excel")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
