"""
Variant format converter for VCF/BCF files.

Supports conversions between:
- VCF (Variant Call Format)
- BCF (Binary VCF)
- TSV (Tabular export)
- BED (Genomic intervals)

Dependencies:
- pysam: VCF/BCF file handling
"""

from pathlib import Path
from typing import Dict, Optional, List, Set
import gzip


class VariantConverter:
    """
    Converter for variant file formats (VCF, BCF).

    Singleton pattern for efficiency.
    """

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize variant converter."""
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
    # VCF/BCF Conversions
    # ============================================================

    def convert_vcf_to_bcf(self, source: str, target: str, index: bool = False):
        """
        Convert VCF to BCF format.

        Args:
            source: Input VCF file path
            target: Output BCF file path
            index: Create BCF index (.csi)
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for VCF/BCF conversion")

        source_path = Path(source)
        target_path = Path(target)

        # Read VCF and write BCF
        with self.pysam.VariantFile(source_path) as vcf_in:
            with self.pysam.VariantFile(
                target_path, "wb", header=vcf_in.header
            ) as bcf_out:
                for record in vcf_in:
                    bcf_out.write(record)

        # Create index if requested
        if index:
            self.pysam.tabix_index(str(target_path), preset="vcf", force=True)

    def convert_bcf_to_vcf(self, source: str, target: str, compress: bool = False):
        """
        Convert BCF to VCF format.

        Args:
            source: Input BCF file path
            target: Output VCF file path
            compress: Compress output with gzip
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for VCF/BCF conversion")

        source_path = Path(source)
        target_path = Path(target)

        # Read BCF and write VCF
        mode = "wz" if compress else "w"
        with self.pysam.VariantFile(source_path) as bcf_in:
            with self.pysam.VariantFile(
                target_path, mode, header=bcf_in.header
            ) as vcf_out:
                for record in bcf_in:
                    vcf_out.write(record)

    def convert_vcf_to_tsv(
        self,
        source: str,
        target: str,
        max_variants: Optional[int] = None,
        include_info: bool = True,
        include_samples: bool = False,
    ):
        """
        Convert VCF to TSV format for analysis.

        Args:
            source: Input VCF file path
            target: Output TSV file path
            max_variants: Maximum variants to export (None for all)
            include_info: Include INFO fields
            include_samples: Include sample genotypes
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for VCF conversion")

        source_path = Path(source)
        target_path = Path(target)

        with self.pysam.VariantFile(source_path) as vcf_in:
            with open(target_path, "w") as tsv_out:
                # Build header
                header = ["chrom", "pos", "id", "ref", "alt", "qual", "filter"]

                if include_info:
                    # Get INFO field names from header
                    info_fields = [rec.name for rec in vcf_in.header.info.values()]
                    header.extend(info_fields)

                if include_samples:
                    # Get sample names
                    samples = list(vcf_in.header.samples)
                    for sample in samples:
                        header.extend([f"{sample}_GT", f"{sample}_DP", f"{sample}_GQ"])

                tsv_out.write("\t".join(header) + "\n")

                # Write variants
                count = 0
                for record in vcf_in:
                    if max_variants and count >= max_variants:
                        break

                    row = [
                        record.chrom,
                        str(record.pos),
                        record.id if record.id else ".",
                        record.ref,
                        (
                            ",".join([str(alt) for alt in record.alts])
                            if record.alts
                            else "."
                        ),
                        str(record.qual) if record.qual is not None else ".",
                        ",".join(record.filter) if record.filter else "PASS",
                    ]

                    if include_info:
                        for field in info_fields:
                            value = record.info.get(field, ".")
                            if isinstance(value, (list, tuple)):
                                value = ",".join(map(str, value))
                            row.append(str(value))

                    if include_samples:
                        for sample in samples:
                            sample_data = record.samples[sample]
                            gt = sample_data.get("GT", ".")
                            dp = sample_data.get("DP", ".")
                            gq = sample_data.get("GQ", ".")

                            # Format genotype
                            if gt and gt != ".":
                                gt_str = "/".join(map(str, gt))
                            else:
                                gt_str = "."

                            row.extend([gt_str, str(dp), str(gq)])

                    tsv_out.write("\t".join(row) + "\n")
                    count += 1

    def convert_vcf_to_bed(
        self,
        source: str,
        target: str,
        extend: int = 0,
        variant_type: Optional[str] = None,
    ):
        """
        Convert VCF to BED format (variant positions as intervals).

        Args:
            source: Input VCF file path
            target: Output BED file path
            extend: Extend intervals by N bp on each side
            variant_type: Filter by variant type (SNP, INDEL, etc.)
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for VCF conversion")

        source_path = Path(source)
        target_path = Path(target)

        with self.pysam.VariantFile(source_path) as vcf_in:
            with open(target_path, "w") as bed_out:
                for record in vcf_in:
                    # Filter by variant type if specified
                    if variant_type:
                        if variant_type == "SNP" and not self._is_snp(record):
                            continue
                        elif variant_type == "INDEL" and not self._is_indel(record):
                            continue

                    # Calculate interval
                    start = max(
                        0, record.pos - 1 - extend
                    )  # VCF is 1-based, BED is 0-based
                    end = record.pos + len(record.ref) - 1 + extend

                    # Build BED line
                    name = (
                        record.id if record.id else f"var_{record.chrom}_{record.pos}"
                    )
                    score = int(record.qual) if record.qual is not None else 0

                    bed_line = f"{record.chrom}\t{start}\t{end}\t{name}\t{score}\n"
                    bed_out.write(bed_line)

    def _is_snp(self, record) -> bool:
        """Check if variant is a SNP."""
        if not record.alts:
            return False
        return len(record.ref) == 1 and all(len(str(alt)) == 1 for alt in record.alts)

    def _is_indel(self, record) -> bool:
        """Check if variant is an INDEL."""
        if not record.alts:
            return False
        return len(record.ref) != 1 or any(len(str(alt)) != 1 for alt in record.alts)

    # ============================================================
    # VCF Operations
    # ============================================================

    def filter_vcf(
        self,
        source: str,
        target: str,
        min_qual: float = 0,
        pass_only: bool = False,
        regions: Optional[List[str]] = None,
        variant_types: Optional[Set[str]] = None,
    ):
        """
        Filter VCF file by quality, filter status, and regions.

        Args:
            source: Input VCF file path
            target: Output filtered VCF file path
            min_qual: Minimum variant quality
            pass_only: Only include PASS variants
            regions: List of regions (chr:start-end)
            variant_types: Set of variant types to include (SNP, INDEL)
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for VCF filtering")

        with self.pysam.VariantFile(source) as vcf_in:
            with self.pysam.VariantFile(target, "w", header=vcf_in.header) as vcf_out:
                for record in vcf_in:
                    # Filter by quality
                    if record.qual is not None and record.qual < min_qual:
                        continue

                    # Filter by PASS status
                    if pass_only and record.filter and len(record.filter) > 0:
                        if "PASS" not in record.filter:
                            continue

                    # Filter by variant type
                    if variant_types:
                        if "SNP" in variant_types and not self._is_snp(record):
                            continue
                        if "INDEL" in variant_types and not self._is_indel(record):
                            continue

                    # Filter by regions (simplified - full implementation would use pysam.TabixFile)
                    if regions:
                        # This is a basic implementation
                        # For production, use indexed VCF and fetch specific regions
                        pass

                    vcf_out.write(record)

    def index_vcf(self, vcf_file: str):
        """
        Create VCF index file (.tbi or .csi).

        Args:
            vcf_file: VCF file to index (must be bgzipped)
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for VCF indexing")

        self.pysam.tabix_index(str(vcf_file), preset="vcf", force=True)

    # ============================================================
    # Statistics
    # ============================================================

    def get_variant_stats(self, file_path: str) -> Dict:
        """
        Get variant statistics from VCF/BCF file.

        Args:
            file_path: VCF or BCF file path

        Returns:
            Dictionary with statistics
        """
        if not self.has_pysam:
            raise ImportError("pysam is required for variant statistics")

        stats = {
            "total_variants": 0,
            "snps": 0,
            "indels": 0,
            "multiallelic": 0,
            "pass_variants": 0,
            "filtered_variants": 0,
            "chromosomes": set(),
            "variant_types": {},
            "qual_distribution": [],
        }

        with self.pysam.VariantFile(file_path) as vcf_in:
            # Get sample count
            stats["num_samples"] = len(list(vcf_in.header.samples))

            for record in vcf_in:
                stats["total_variants"] += 1
                stats["chromosomes"].add(record.chrom)

                # Variant types
                if self._is_snp(record):
                    stats["snps"] += 1
                elif self._is_indel(record):
                    stats["indels"] += 1

                # Multiallelic
                if record.alts and len(record.alts) > 1:
                    stats["multiallelic"] += 1

                # Filter status
                if (
                    not record.filter
                    or len(record.filter) == 0
                    or "PASS" in record.filter
                ):
                    stats["pass_variants"] += 1
                else:
                    stats["filtered_variants"] += 1

                # Quality
                if record.qual is not None:
                    stats["qual_distribution"].append(record.qual)

        # Convert chromosome set to list
        stats["chromosomes"] = sorted(list(stats["chromosomes"]))
        stats["num_chromosomes"] = len(stats["chromosomes"])

        # Calculate quality statistics
        if stats["qual_distribution"]:
            stats["avg_qual"] = sum(stats["qual_distribution"]) / len(
                stats["qual_distribution"]
            )
            stats["min_qual"] = min(stats["qual_distribution"])
            stats["max_qual"] = max(stats["qual_distribution"])

        # Calculate percentages
        if stats["total_variants"] > 0:
            stats["snp_percent"] = (stats["snps"] / stats["total_variants"]) * 100
            stats["indel_percent"] = (stats["indels"] / stats["total_variants"]) * 100
            stats["pass_percent"] = (
                stats["pass_variants"] / stats["total_variants"]
            ) * 100

        return stats


def get_variant_converter() -> VariantConverter:
    """Get singleton instance of VariantConverter."""
    return VariantConverter()
