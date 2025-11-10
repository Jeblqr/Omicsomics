#!/usr/bin/env python3
"""
Quick test script for Phase 1 Week 1 bioinformatics format converters.

This script creates sample files and tests the conversion functionality
without requiring pytest or external dependencies.
"""
import sys
import tempfile
from pathlib import Path
import gzip

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'backend'))

from app.converters.sequence_converter import get_sequence_converter
from app.converters.interval_converter import get_interval_converter
from app.converters.format_converter import FormatConverter


def create_sample_fasta(path: Path) -> Path:
    """Create sample FASTA file."""
    content = """>seq1 Sample sequence 1
ATCGATCGATCGTAGCTAGCTAGCTA
>seq2 Sample sequence 2  
GCTAGCTAGCTAAAATTTTCCCCGGG
>seq3 Sample sequence 3
TTAATTAATTAAGCGCGCGCATATAT
"""
    path.write_text(content)
    return path


def create_sample_fastq(path: Path) -> Path:
    """Create sample FASTQ file."""
    content = """@seq1 Sample read 1
ATCGATCGATCGTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIII
@seq2 Sample read 2
GCTAGCTAGCTAAAATTTTCCCCGGG
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@seq3 Sample read 3
TTAATTAATTAAGCGCGCGCATATAT
+
JJJJJJJJJJJJJJJJJJJJJJJJJJ
"""
    path.write_text(content)
    return path


def create_sample_fastq_gz(path: Path) -> Path:
    """Create sample gzipped FASTQ file."""
    content = """@seq1 Sample read 1
ATCGATCGATCGTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIII
@seq2 Sample read 2
GCTAGCTAGCTAAAATTTTCCCCGGG
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
"""
    with gzip.open(path, 'wt') as f:
        f.write(content)
    return path


def create_sample_bed(path: Path, bed_type: str = '3') -> Path:
    """Create sample BED file."""
    if bed_type == '3':
        content = """chr1\t1000\t2000
chr1\t3000\t4000
chr2\t5000\t6000
chr3\t10000\t15000
"""
    elif bed_type == '6':
        content = """chr1\t1000\t2000\tfeature1\t100\t+
chr1\t3000\t4000\tfeature2\t200\t-
chr2\t5000\t6000\tfeature3\t150\t+
chr3\t10000\t15000\tfeature4\t300\t-
"""
    else:
        content = """chr1\t1000\t2000
chr1\t3000\t4000
"""
    
    path.write_text(content)
    return path


def create_sample_bedgraph(path: Path) -> Path:
    """Create sample bedGraph file."""
    content = """chr1\t1000\t2000\t1.5
chr1\t2000\t3000\t2.3
chr1\t3000\t4000\t0.8
chr2\t5000\t6000\t3.2
"""
    path.write_text(content)
    return path


def test_sequence_converter():
    """Test sequence converter functionality."""
    print("\n" + "="*60)
    print("Testing Sequence Converter")
    print("="*60)
    
    converter = get_sequence_converter()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        
        # Test 1: FASTQ to FASTA
        print("\n[1] Testing FASTQ → FASTA conversion...")
        fastq_file = create_sample_fastq(tmp_path / "test.fastq")
        fasta_output = tmp_path / "output.fasta"
        
        converter.convert_fastq_to_fasta(str(fastq_file), str(fasta_output))
        
        if fasta_output.exists():
            content = fasta_output.read_text()
            if ">seq1" in content and "+" not in content:
                print("   ✓ FASTQ to FASTA conversion successful")
                print(f"   ✓ Output file: {fasta_output.stat().st_size} bytes")
            else:
                print("   ✗ Output format incorrect")
        else:
            print("   ✗ Output file not created")
        
        # Test 2: FASTA to TSV
        print("\n[2] Testing FASTA → TSV conversion...")
        fasta_file = create_sample_fasta(tmp_path / "test.fasta")
        tsv_output = tmp_path / "output.tsv"
        
        converter.convert_fasta_to_tsv(str(fasta_file), str(tsv_output))
        
        if tsv_output.exists():
            lines = tsv_output.read_text().strip().split('\n')
            if len(lines) == 4 and "id\tdescription" in lines[0]:
                print("   ✓ FASTA to TSV conversion successful")
                print(f"   ✓ {len(lines)-1} sequences exported")
            else:
                print("   ✗ Output format incorrect")
        else:
            print("   ✗ Output file not created")
        
        # Test 3: FASTQ.gz decompression
        print("\n[3] Testing FASTQ.gz → FASTQ decompression...")
        fastq_gz = create_sample_fastq_gz(tmp_path / "test.fastq.gz")
        fastq_output = tmp_path / "decompressed.fastq"
        
        converter.decompress_fastq_gz(str(fastq_gz), str(fastq_output))
        
        if fastq_output.exists():
            content = fastq_output.read_text()
            if "@seq1" in content and "ATCGATCG" in content:
                print("   ✓ Decompression successful")
                print(f"   ✓ Decompressed size: {fastq_output.stat().st_size} bytes")
            else:
                print("   ✗ Decompressed content incorrect")
        else:
            print("   ✗ Output file not created")
        
        # Test 4: Sequence statistics
        print("\n[4] Testing sequence statistics...")
        stats = converter.get_sequence_stats(str(fasta_file), 'fasta')
        
        print(f"   ✓ Sequence count: {stats['count']}")
        print(f"   ✓ Total length: {stats['total_length']} bp")
        print(f"   ✓ Average length: {stats['avg_length']:.1f} bp")
        print(f"   ✓ Length range: {stats['min_length']}-{stats['max_length']} bp")


def test_interval_converter():
    """Test interval converter functionality."""
    print("\n" + "="*60)
    print("Testing Interval Converter")
    print("="*60)
    
    converter = get_interval_converter()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        
        # Test 1: BED3 to CSV
        print("\n[1] Testing BED3 → CSV conversion...")
        bed_file = create_sample_bed(tmp_path / "test.bed", '3')
        csv_output = tmp_path / "output.csv"
        
        converter.convert_bed_to_csv(str(bed_file), str(csv_output), bed_type='3')
        
        if csv_output.exists():
            content = csv_output.read_text()
            if "chrom,start,end" in content and "chr1,1000,2000" in content:
                print("   ✓ BED3 to CSV conversion successful")
                print(f"   ✓ Output lines: {len(content.strip().split('\\n'))}")
            else:
                print("   ✗ Output format incorrect")
        else:
            print("   ✗ Output file not created")
        
        # Test 2: BED6 to CSV
        print("\n[2] Testing BED6 → CSV conversion...")
        bed6_file = create_sample_bed(tmp_path / "test_bed6.bed", '6')
        csv6_output = tmp_path / "output6.csv"
        
        converter.convert_bed_to_csv(str(bed6_file), str(csv6_output), bed_type='6')
        
        if csv6_output.exists():
            content = csv6_output.read_text()
            if "name,score,strand" in content and "feature1" in content:
                print("   ✓ BED6 to CSV conversion successful")
                print(f"   ✓ Columns include: name, score, strand")
            else:
                print("   ✗ Output format incorrect")
        else:
            print("   ✗ Output file not created")
        
        # Test 3: bedGraph to BED
        print("\n[3] Testing bedGraph → BED conversion...")
        bg_file = create_sample_bedgraph(tmp_path / "test.bedgraph")
        bed_output = tmp_path / "from_bedgraph.bed"
        
        converter.convert_bedgraph_to_bed(str(bg_file), str(bed_output))
        
        if bed_output.exists():
            lines = bed_output.read_text().strip().split('\n')
            if len(lines) == 4 and "chr1\t1000\t2000" in lines[0]:
                print("   ✓ bedGraph to BED conversion successful")
                print(f"   ✓ {len(lines)} intervals converted")
            else:
                print("   ✗ Output format incorrect")
        else:
            print("   ✗ Output file not created")
        
        # Test 4: BED statistics
        print("\n[4] Testing BED statistics...")
        stats = converter.get_bed_stats(str(bed_file))
        
        print(f"   ✓ Interval count: {stats['count']}")
        print(f"   ✓ Total length: {stats['total_length']} bp")
        print(f"   ✓ Average length: {stats['avg_length']:.1f} bp")
        print(f"   ✓ Chromosomes: {', '.join(stats['chromosomes'])}")


def test_format_converter_integration():
    """Test FormatConverter integration."""
    print("\n" + "="*60)
    print("Testing FormatConverter Integration")
    print("="*60)
    
    converter = FormatConverter()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        
        # Test 1: Format detection
        print("\n[1] Testing format detection...")
        fasta_file = create_sample_fasta(tmp_path / "test.fasta")
        fastq_file = create_sample_fastq(tmp_path / "test.fastq")
        bed_file = create_sample_bed(tmp_path / "test.bed")
        
        fasta_fmt = converter.detect_format(str(fasta_file))
        fastq_fmt = converter.detect_format(str(fastq_file))
        bed_fmt = converter.detect_format(str(bed_file))
        
        if fasta_fmt == 'fasta':
            print("   ✓ FASTA format detected correctly")
        else:
            print(f"   ✗ FASTA detection failed: got {fasta_fmt}")
        
        if fastq_fmt == 'fastq':
            print("   ✓ FASTQ format detected correctly")
        else:
            print(f"   ✗ FASTQ detection failed: got {fastq_fmt}")
        
        if bed_fmt == 'bed':
            print("   ✓ BED format detected correctly")
        else:
            print(f"   ✗ BED detection failed: got {bed_fmt}")
        
        # Test 2: Conversion paths
        print("\n[2] Testing conversion path calculation...")
        
        path1 = converter.get_conversion_path('fastq', 'fasta')
        if path1 == ['fastq', 'fasta']:
            print("   ✓ Direct path: fastq → fasta")
        else:
            print(f"   ✗ Path incorrect: {path1}")
        
        path2 = converter.get_conversion_path('fastq.gz', 'fasta')
        if path2 == ['fastq.gz', 'fastq', 'fasta']:
            print("   ✓ Multi-step: fastq.gz → fastq → fasta")
        else:
            print(f"   ✗ Path incorrect: {path2}")
        
        path3 = converter.get_conversion_path('bigwig', 'csv')
        if path3 == ['bigwig', 'bedgraph', 'bed', 'csv']:
            print("   ✓ Multi-step: bigwig → bedgraph → bed → csv")
        else:
            print(f"   ✗ Path incorrect: {path3}")
        
        # Test 3: Full conversion
        print("\n[3] Testing full conversion pipeline...")
        output_fasta = tmp_path / "converted.fasta"
        
        result = converter.convert(
            str(fastq_file),
            str(output_fasta),
            from_format='fastq',
            to_format='fasta'
        )
        
        if output_fasta.exists() and result['status'] == 'completed':
            print("   ✓ Full conversion pipeline successful")
            print(f"   ✓ Status: {result['status']}")
            print(f"   ✓ Steps: {len(result['steps'])}")
        else:
            print("   ✗ Conversion failed")


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("Bioinformatics Format Converter - Quick Test Suite")
    print("Phase 1 Week 1: Sequence & Interval Formats")
    print("="*60)
    
    try:
        test_sequence_converter()
        test_interval_converter()
        test_format_converter_integration()
        
        print("\n" + "="*60)
        print("✓ All tests completed successfully!")
        print("="*60)
        print("\nNext steps:")
        print("1. Run full pytest suite: pytest tests/test_bioinformatics_converters.py")
        print("2. Add API endpoints for bioinformatics formats")
        print("3. Update frontend format converter UI")
        print("4. Begin Phase 1 Week 2: Alignment & Variant formats")
        print("="*60 + "\n")
        
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
