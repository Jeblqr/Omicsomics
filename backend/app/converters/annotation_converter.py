"""
Annotation format converter for GTF/GFF3 files.

Supports conversions between:
- GTF (Gene Transfer Format)
- GFF3 (General Feature Format version 3)
- TSV (Tabular export)
- BED (Genomic intervals)

Dependencies:
- pandas: Data manipulation
"""

from pathlib import Path
from typing import Dict, Optional, List, Set
import csv
import re


class AnnotationConverter:
    """
    Converter for annotation file formats (GTF, GFF3).

    Singleton pattern for efficiency.
    """

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize annotation converter."""
        if not hasattr(self, "_initialized"):
            self._initialized = True

    # ============================================================
    # GTF/GFF3 Parsing
    # ============================================================

    def _parse_gtf_attributes(self, attr_string: str) -> Dict:
        """Parse GTF attribute string into dictionary."""
        attributes = {}

        # GTF format: key "value"; key "value";
        pattern = r'(\w+)\s+"([^"]*)";?'
        matches = re.findall(pattern, attr_string)

        for key, value in matches:
            attributes[key] = value

        return attributes

    def _parse_gff3_attributes(self, attr_string: str) -> Dict:
        """Parse GFF3 attribute string into dictionary."""
        attributes = {}

        # GFF3 format: key=value;key=value
        for item in attr_string.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attributes[key.strip()] = value.strip()

        return attributes

    def _format_gtf_attributes(self, attributes: Dict) -> str:
        """Format dictionary into GTF attribute string."""
        items = []
        for key, value in attributes.items():
            items.append(f'{key} "{value}"')
        return "; ".join(items) + ";"

    def _format_gff3_attributes(self, attributes: Dict) -> str:
        """Format dictionary into GFF3 attribute string."""
        items = []
        for key, value in attributes.items():
            items.append(f"{key}={value}")
        return ";".join(items)

    # ============================================================
    # GTF/GFF3 Conversions
    # ============================================================

    def convert_gtf_to_gff3(self, source: str, target: str):
        """
        Convert GTF to GFF3 format.

        Args:
            source: Input GTF file path
            target: Output GFF3 file path
        """
        source_path = Path(source)
        target_path = Path(target)

        with open(source_path, "r") as gtf_in, open(target_path, "w") as gff3_out:
            # Write GFF3 header
            gff3_out.write("##gff-version 3\n")

            for line in gtf_in:
                # Skip comments and empty lines
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                # Parse GTF attributes
                gtf_attrs = self._parse_gtf_attributes(fields[8])

                # Convert to GFF3 attributes
                gff3_attrs = {}
                if "gene_id" in gtf_attrs:
                    gff3_attrs["gene_id"] = gtf_attrs["gene_id"]
                if "transcript_id" in gtf_attrs:
                    gff3_attrs["transcript_id"] = gtf_attrs["transcript_id"]
                if "gene_name" in gtf_attrs:
                    gff3_attrs["Name"] = gtf_attrs["gene_name"]

                # Add ID based on feature type
                feature_type = fields[2]
                if feature_type == "gene" and "gene_id" in gtf_attrs:
                    gff3_attrs["ID"] = gtf_attrs["gene_id"]
                elif feature_type == "transcript" and "transcript_id" in gtf_attrs:
                    gff3_attrs["ID"] = gtf_attrs["transcript_id"]
                    if "gene_id" in gtf_attrs:
                        gff3_attrs["Parent"] = gtf_attrs["gene_id"]
                elif "transcript_id" in gtf_attrs:
                    gff3_attrs["Parent"] = gtf_attrs["transcript_id"]

                # Copy other attributes
                for key, value in gtf_attrs.items():
                    if key not in ["gene_id", "transcript_id", "gene_name"]:
                        gff3_attrs[key] = value

                # Write GFF3 line
                fields[8] = self._format_gff3_attributes(gff3_attrs)
                gff3_out.write("\t".join(fields) + "\n")

    def convert_gff3_to_gtf(self, source: str, target: str):
        """
        Convert GFF3 to GTF format.

        Args:
            source: Input GFF3 file path
            target: Output GTF file path
        """
        source_path = Path(source)
        target_path = Path(target)

        with open(source_path, "r") as gff3_in, open(target_path, "w") as gtf_out:
            for line in gff3_in:
                # Skip comments and empty lines
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                # Parse GFF3 attributes
                gff3_attrs = self._parse_gff3_attributes(fields[8])

                # Convert to GTF attributes
                gtf_attrs = {}
                if "gene_id" in gff3_attrs:
                    gtf_attrs["gene_id"] = gff3_attrs["gene_id"]
                if "transcript_id" in gff3_attrs:
                    gtf_attrs["transcript_id"] = gff3_attrs["transcript_id"]
                if "Name" in gff3_attrs:
                    gtf_attrs["gene_name"] = gff3_attrs["Name"]

                # Handle ID/Parent for gene/transcript hierarchy
                if "ID" in gff3_attrs and "gene_id" not in gtf_attrs:
                    gtf_attrs["gene_id"] = gff3_attrs["ID"]
                if "Parent" in gff3_attrs and "transcript_id" not in gtf_attrs:
                    gtf_attrs["transcript_id"] = gff3_attrs["Parent"]

                # Copy other attributes
                for key, value in gff3_attrs.items():
                    if key not in ["ID", "Parent", "Name", "gene_id", "transcript_id"]:
                        gtf_attrs[key] = value

                # Write GTF line
                fields[8] = self._format_gtf_attributes(gtf_attrs)
                gtf_out.write("\t".join(fields) + "\n")

    def convert_gtf_to_tsv(
        self,
        source: str,
        target: str,
        max_features: Optional[int] = None,
        feature_types: Optional[Set[str]] = None,
    ):
        """
        Convert GTF to TSV format for analysis.

        Args:
            source: Input GTF file path
            target: Output TSV file path
            max_features: Maximum features to export (None for all)
            feature_types: Set of feature types to include (gene, exon, etc.)
        """
        source_path = Path(source)
        target_path = Path(target)

        with open(source_path, "r") as gtf_in, open(target_path, "w") as tsv_out:
            # Write header
            header = [
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "gene_id",
                "transcript_id",
                "gene_name",
                "gene_type",
                "transcript_name",
                "transcript_type",
            ]
            tsv_out.write("\t".join(header) + "\n")

            count = 0
            for line in gtf_in:
                if line.startswith("#") or not line.strip():
                    continue

                if max_features and count >= max_features:
                    break

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                # Filter by feature type
                if feature_types and fields[2] not in feature_types:
                    continue

                # Parse attributes
                attrs = self._parse_gtf_attributes(fields[8])

                # Build row
                row = fields[:8]  # seqname to frame
                row.extend(
                    [
                        attrs.get("gene_id", ""),
                        attrs.get("transcript_id", ""),
                        attrs.get("gene_name", ""),
                        attrs.get("gene_type", attrs.get("gene_biotype", "")),
                        attrs.get("transcript_name", ""),
                        attrs.get(
                            "transcript_type", attrs.get("transcript_biotype", "")
                        ),
                    ]
                )

                tsv_out.write("\t".join(row) + "\n")
                count += 1

    def convert_gff3_to_tsv(
        self,
        source: str,
        target: str,
        max_features: Optional[int] = None,
        feature_types: Optional[Set[str]] = None,
    ):
        """
        Convert GFF3 to TSV format for analysis.

        Args:
            source: Input GFF3 file path
            target: Output TSV file path
            max_features: Maximum features to export (None for all)
            feature_types: Set of feature types to include
        """
        source_path = Path(source)
        target_path = Path(target)

        with open(source_path, "r") as gff3_in, open(target_path, "w") as tsv_out:
            # Write header
            header = [
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "ID",
                "Parent",
                "Name",
                "gene_id",
                "transcript_id",
            ]
            tsv_out.write("\t".join(header) + "\n")

            count = 0
            for line in gff3_in:
                if line.startswith("#") or not line.strip():
                    continue

                if max_features and count >= max_features:
                    break

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                # Filter by feature type
                if feature_types and fields[2] not in feature_types:
                    continue

                # Parse attributes
                attrs = self._parse_gff3_attributes(fields[8])

                # Build row
                row = fields[:8]  # seqname to frame
                row.extend(
                    [
                        attrs.get("ID", ""),
                        attrs.get("Parent", ""),
                        attrs.get("Name", ""),
                        attrs.get("gene_id", ""),
                        attrs.get("transcript_id", ""),
                    ]
                )

                tsv_out.write("\t".join(row) + "\n")
                count += 1

    def convert_gtf_to_bed(
        self,
        source: str,
        target: str,
        feature_type: str = "exon",
        name_field: str = "gene_name",
    ):
        """
        Convert GTF to BED format.

        Args:
            source: Input GTF file path
            target: Output BED file path
            feature_type: Feature type to extract (gene, exon, etc.)
            name_field: Attribute to use as BED name
        """
        source_path = Path(source)
        target_path = Path(target)

        with open(source_path, "r") as gtf_in, open(target_path, "w") as bed_out:
            for line in gtf_in:
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                # Filter by feature type
                if fields[2] != feature_type:
                    continue

                # Parse attributes
                attrs = self._parse_gtf_attributes(fields[8])

                # Build BED line (BED6)
                chrom = fields[0]
                start = int(fields[3]) - 1  # GTF is 1-based, BED is 0-based
                end = int(fields[4])
                name = attrs.get(name_field, attrs.get("gene_id", "."))
                score = fields[5] if fields[5] != "." else "0"
                strand = fields[6]

                bed_line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"
                bed_out.write(bed_line)

    def convert_gff3_to_bed(
        self,
        source: str,
        target: str,
        feature_type: str = "exon",
        name_field: str = "Name",
    ):
        """
        Convert GFF3 to BED format.

        Args:
            source: Input GFF3 file path
            target: Output BED file path
            feature_type: Feature type to extract
            name_field: Attribute to use as BED name
        """
        source_path = Path(source)
        target_path = Path(target)

        with open(source_path, "r") as gff3_in, open(target_path, "w") as bed_out:
            for line in gff3_in:
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                # Filter by feature type
                if fields[2] != feature_type:
                    continue

                # Parse attributes
                attrs = self._parse_gff3_attributes(fields[8])

                # Build BED line (BED6)
                chrom = fields[0]
                start = int(fields[3]) - 1  # GFF3 is 1-based, BED is 0-based
                end = int(fields[4])
                name = attrs.get(name_field, attrs.get("ID", "."))
                score = fields[5] if fields[5] != "." else "0"
                strand = fields[6]

                bed_line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"
                bed_out.write(bed_line)

    # ============================================================
    # Statistics
    # ============================================================

    def get_annotation_stats(self, file_path: str, file_format: str = "gtf") -> Dict:
        """
        Get annotation statistics from GTF/GFF3 file.

        Args:
            file_path: GTF or GFF3 file path
            file_format: 'gtf' or 'gff3'

        Returns:
            Dictionary with statistics
        """
        stats = {
            "total_features": 0,
            "feature_types": {},
            "chromosomes": set(),
            "genes": set(),
            "transcripts": set(),
            "strands": {"+": 0, "-": 0, ".": 0},
        }

        parse_attrs = (
            self._parse_gtf_attributes
            if file_format == "gtf"
            else self._parse_gff3_attributes
        )

        with open(file_path, "r") as ann_file:
            for line in ann_file:
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                stats["total_features"] += 1

                # Feature type
                feature_type = fields[2]
                stats["feature_types"][feature_type] = (
                    stats["feature_types"].get(feature_type, 0) + 1
                )

                # Chromosome
                stats["chromosomes"].add(fields[0])

                # Strand
                strand = fields[6]
                if strand in stats["strands"]:
                    stats["strands"][strand] += 1

                # Parse attributes
                attrs = parse_attrs(fields[8])

                if file_format == "gtf":
                    if "gene_id" in attrs:
                        stats["genes"].add(attrs["gene_id"])
                    if "transcript_id" in attrs:
                        stats["transcripts"].add(attrs["transcript_id"])
                else:  # gff3
                    if "gene_id" in attrs:
                        stats["genes"].add(attrs["gene_id"])
                    elif fields[2] == "gene" and "ID" in attrs:
                        stats["genes"].add(attrs["ID"])

                    if "transcript_id" in attrs:
                        stats["transcripts"].add(attrs["transcript_id"])
                    elif fields[2] in ["mRNA", "transcript"] and "ID" in attrs:
                        stats["transcripts"].add(attrs["ID"])

        # Convert sets to counts
        stats["num_genes"] = len(stats["genes"])
        stats["num_transcripts"] = len(stats["transcripts"])
        stats["num_chromosomes"] = len(stats["chromosomes"])
        stats["chromosomes"] = sorted(list(stats["chromosomes"]))

        # Remove large sets to keep stats dict manageable
        del stats["genes"]
        del stats["transcripts"]

        return stats


def get_annotation_converter() -> AnnotationConverter:
    """Get singleton instance of AnnotationConverter."""
    return AnnotationConverter()
