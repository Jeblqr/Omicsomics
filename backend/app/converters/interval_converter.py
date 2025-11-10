"""
Genomic interval format converter.

Supports conversions between:
- BED (BED3, BED6, BED12)
- bedGraph
- BigWig
- TSV/CSV
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import pandas as pd

logger = logging.getLogger(__name__)


class IntervalConverter:
    """
    Handle genomic interval format conversions.
    """

    def __init__(self):
        # Try to import optional dependencies
        try:
            import pyBigWig

            self.pyBigWig = pyBigWig
            self.pybigwig_available = True
        except ImportError:
            logger.warning(
                "pyBigWig not available. BigWig conversion will be disabled."
            )
            self.pybigwig_available = False

        try:
            import pybedtools

            self.pybedtools = pybedtools
            self.pybedtools_available = True
        except ImportError:
            logger.warning(
                "pybedtools not available. Some BED conversions will be limited."
            )
            self.pybedtools_available = False

    def convert_bed_to_csv(
        self, source_path: str, target_path: str, bed_type: str = "bed6"
    ) -> Dict:
        """
        Convert BED to CSV.

        Args:
            source_path: Path to input BED file
            target_path: Path to output CSV file
            bed_type: BED format ('bed3', 'bed6', or 'bed12')

        Returns:
            Conversion result info
        """
        start_time = datetime.now()

        # Define column names based on BED type
        if bed_type == "bed3":
            columns = ["chrom", "start", "end"]
        elif bed_type == "bed6":
            columns = ["chrom", "start", "end", "name", "score", "strand"]
        elif bed_type == "bed12":
            columns = [
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts",
            ]
        else:
            raise ValueError(f"Unknown BED type: {bed_type}")

        # Read BED file
        df = pd.read_csv(source_path, sep="\t", header=None, names=columns, comment="#")

        # Write CSV
        df.to_csv(target_path, index=False)

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        return {
            "success": True,
            "intervals_converted": len(df),
            "duration_seconds": duration,
            "source_format": bed_type,
            "target_format": "csv",
            "columns": columns,
        }

    def convert_csv_to_bed(
        self, source_path: str, target_path: str, bed_type: str = "bed6"
    ) -> Dict:
        """
        Convert CSV to BED.

        Args:
            source_path: Path to input CSV file
            target_path: Path to output BED file
            bed_type: Target BED format ('bed3', 'bed6', or 'bed12')

        Returns:
            Conversion result info
        """
        start_time = datetime.now()

        # Read CSV
        df = pd.read_csv(source_path)

        # Select columns based on BED type
        if bed_type == "bed3":
            required_cols = ["chrom", "start", "end"]
        elif bed_type == "bed6":
            required_cols = ["chrom", "start", "end", "name", "score", "strand"]
        elif bed_type == "bed12":
            required_cols = [
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts",
            ]
        else:
            raise ValueError(f"Unknown BED type: {bed_type}")

        # Check if required columns exist
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Select and write BED
        bed_df = df[required_cols]
        bed_df.to_csv(target_path, sep="\t", header=False, index=False)

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        return {
            "success": True,
            "intervals_converted": len(bed_df),
            "duration_seconds": duration,
            "source_format": "csv",
            "target_format": bed_type,
        }

    def convert_bedgraph_to_bed(self, source_path: str, target_path: str) -> Dict:
        """
        Convert bedGraph to BED4.

        bedGraph format: chrom start end value
        BED4 format: chrom start end name (value as name)

        Args:
            source_path: Path to input bedGraph file
            target_path: Path to output BED file

        Returns:
            Conversion result info
        """
        start_time = datetime.now()

        # Read bedGraph
        df = pd.read_csv(
            source_path,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "value"],
            comment="#",
        )

        # Convert value to string for BED name field
        df["name"] = df["value"].astype(str)

        # Write BED4
        bed_df = df[["chrom", "start", "end", "name"]]
        bed_df.to_csv(target_path, sep="\t", header=False, index=False)

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        return {
            "success": True,
            "intervals_converted": len(df),
            "duration_seconds": duration,
            "source_format": "bedgraph",
            "target_format": "bed4",
        }

    def convert_bed_to_bedgraph(self, source_path: str, target_path: str) -> Dict:
        """
        Convert BED to bedGraph.

        Uses the score column (5th column) as the bedGraph value.

        Args:
            source_path: Path to input BED file
            target_path: Path to output bedGraph file

        Returns:
            Conversion result info
        """
        start_time = datetime.now()

        # Read BED file
        df = pd.read_csv(source_path, sep="\t", header=None, comment="#")

        # Check if score column exists (column 4, 0-indexed)
        if len(df.columns) < 5:
            raise ValueError(
                "BED file must have at least 5 columns for bedGraph conversion"
            )

        # Extract bedGraph columns: chrom, start, end, score
        bedgraph_df = df.iloc[:, [0, 1, 2, 4]]
        bedgraph_df.columns = ["chrom", "start", "end", "value"]

        # Write bedGraph
        with open(target_path, "w") as f:
            # Write track line
            f.write("track type=bedGraph\n")
            bedgraph_df.to_csv(f, sep="\t", header=False, index=False)

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        return {
            "success": True,
            "intervals_converted": len(bedgraph_df),
            "duration_seconds": duration,
            "source_format": "bed",
            "target_format": "bedgraph",
        }

    def convert_bigwig_to_bedgraph(
        self, source_path: str, target_path: str, chrom: Optional[str] = None
    ) -> Dict:
        """
        Convert BigWig to bedGraph.

        Args:
            source_path: Path to input BigWig file
            target_path: Path to output bedGraph file
            chrom: Optional chromosome to extract (None = all chromosomes)

        Returns:
            Conversion result info
        """
        if not self.pybigwig_available:
            raise RuntimeError("pyBigWig is required for BigWig conversion")

        start_time = datetime.now()

        bw = self.pyBigWig.open(source_path)

        try:
            total_intervals = 0

            with open(target_path, "w") as f:
                # Write track line
                f.write("track type=bedGraph\n")

                # Get chromosomes to process
                if chrom:
                    chroms = [chrom]
                else:
                    chroms = bw.chroms().keys()

                # Extract intervals for each chromosome
                for chr_name in chroms:
                    try:
                        intervals = bw.intervals(chr_name)
                        if intervals:
                            for start, end, value in intervals:
                                f.write(f"{chr_name}\t{start}\t{end}\t{value}\n")
                                total_intervals += 1
                    except RuntimeError:
                        # Chromosome not in BigWig file
                        logger.warning(f"Chromosome {chr_name} not found in BigWig")
                        continue

            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()

            return {
                "success": True,
                "intervals_converted": total_intervals,
                "chromosomes_processed": len(chroms) if not chrom else 1,
                "duration_seconds": duration,
                "source_format": "bigwig",
                "target_format": "bedgraph",
            }
        finally:
            bw.close()

    def convert_bedgraph_to_bigwig(
        self, source_path: str, target_path: str, chrom_sizes: Dict[str, int]
    ) -> Dict:
        """
        Convert bedGraph to BigWig.

        Args:
            source_path: Path to input bedGraph file
            target_path: Path to output BigWig file
            chrom_sizes: Dictionary mapping chromosome names to sizes

        Returns:
            Conversion result info
        """
        if not self.pybigwig_available:
            raise RuntimeError("pyBigWig is required for BigWig conversion")

        start_time = datetime.now()

        # Read bedGraph
        df = pd.read_csv(
            source_path,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "value"],
            comment="track",
        )

        # Open BigWig for writing
        bw = self.pyBigWig.open(target_path, "w")

        try:
            # Add header with chromosome sizes
            bw.addHeader(list(chrom_sizes.items()))

            # Group by chromosome and add entries
            for chrom in df["chrom"].unique():
                chrom_df = df[df["chrom"] == chrom]

                chroms = [chrom] * len(chrom_df)
                starts = chrom_df["start"].tolist()
                ends = chrom_df["end"].tolist()
                values = chrom_df["value"].tolist()

                bw.addEntries(chroms, starts, ends=ends, values=values)

            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()

            return {
                "success": True,
                "intervals_converted": len(df),
                "chromosomes_processed": df["chrom"].nunique(),
                "duration_seconds": duration,
                "source_format": "bedgraph",
                "target_format": "bigwig",
            }
        finally:
            bw.close()

    def get_bed_stats(self, file_path: str) -> Dict:
        """
        Get statistics about a BED file.

        Args:
            file_path: Path to BED file

        Returns:
            Statistics dictionary
        """
        df = pd.read_csv(file_path, sep="\t", header=None, comment="#")

        # Calculate interval lengths
        if len(df.columns) >= 3:
            df["length"] = df[2] - df[1]  # end - start

            stats = {
                "num_intervals": len(df),
                "num_chromosomes": df[0].nunique(),
                "total_bases": df["length"].sum(),
                "avg_length": df["length"].mean(),
                "min_length": df["length"].min(),
                "max_length": df["length"].max(),
                "chromosomes": df[0].unique().tolist(),
            }

            return stats
        else:
            return {
                "num_intervals": len(df),
                "error": "Insufficient columns for BED format",
            }


# Singleton instance
_interval_converter = None


def get_interval_converter() -> IntervalConverter:
    """Get singleton interval converter instance."""
    global _interval_converter
    if _interval_converter is None:
        _interval_converter = IntervalConverter()
    return _interval_converter
