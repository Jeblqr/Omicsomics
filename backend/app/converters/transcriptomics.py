"""Transcriptomics data format converters."""

from pathlib import Path
from typing import Dict, List, Any
import csv
import pandas as pd

from app.converters.base import BaseConverter, ConverterFactory
from app.schemas.unified_format import (
    UnifiedData,
    UnifiedDataRecord,
    OmicsType,
    TranscriptomicsData,
)


@ConverterFactory.register(OmicsType.TRANSCRIPTOMICS)
class TranscriptomicsConverter(BaseConverter):
    """Converter for transcriptomics data formats (RNA-seq counts, expression matrices)."""

    async def to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> UnifiedData:
        """Convert transcriptomics file to unified format."""

        source_format_lower = source_format.lower()

        if source_format_lower in ["csv"]:
            return await self._csv_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif source_format_lower in ["tsv", "txt", "counts"]:
            return await self._tsv_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif source_format_lower in ["xlsx", "xls"]:
            return await self._excel_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        else:
            raise ValueError(f"Unsupported source format: {source_format}")

    async def _csv_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> TranscriptomicsData:
        """Convert CSV expression matrix to unified format."""

        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )

        records = []
        headers = []

        with open(file_path, "r") as f:
            reader = csv.DictReader(f)
            headers = reader.fieldnames or []

            for idx, row in enumerate(reader):
                # First column is usually gene/feature ID
                gene_id = (
                    row.get(headers[0], f"feature_{idx}")
                    if headers
                    else f"feature_{idx}"
                )

                # Create values dict from row
                values = dict(row)

                record = UnifiedDataRecord(
                    id=gene_id,
                    values=values,
                )
                records.append(record)

        statistics = {
            "total_features": len(records),
            "total_samples": len(headers) - 1 if headers else 0,  # Exclude ID column
            "file_size_bytes": file_path.stat().st_size,
        }

        return TranscriptomicsData(
            metadata=metadata,
            headers=list(headers),
            records=records,
            statistics=statistics,
        )

    async def _tsv_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> TranscriptomicsData:
        """Convert TSV/TXT expression matrix to unified format."""

        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )

        records = []
        headers = []

        with open(file_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            headers = reader.fieldnames or []

            for idx, row in enumerate(reader):
                # First column is usually gene/feature ID
                gene_id = (
                    row.get(headers[0], f"feature_{idx}")
                    if headers
                    else f"feature_{idx}"
                )

                # Create values dict, convert numeric values
                values = {}
                for key, value in row.items():
                    try:
                        # Try to convert to float
                        values[key] = float(value)
                    except (ValueError, TypeError):
                        # Keep as string if not numeric
                        values[key] = value

                record = UnifiedDataRecord(
                    id=gene_id,
                    values=values,
                )
                records.append(record)

        statistics = {
            "total_features": len(records),
            "total_samples": len(headers) - 1 if headers else 0,
            "file_size_bytes": file_path.stat().st_size,
        }

        return TranscriptomicsData(
            metadata=metadata,
            headers=list(headers),
            records=records,
            statistics=statistics,
        )

    async def _excel_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> TranscriptomicsData:
        """Convert Excel expression matrix to unified format."""

        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )

        try:
            # Read Excel file
            df = pd.read_excel(file_path)
            headers = list(df.columns)

            records = []
            for idx, row in df.iterrows():
                gene_id = str(row.iloc[0]) if len(row) > 0 else f"feature_{idx}"

                values = {}
                for col in df.columns:
                    values[col] = row[col]

                record = UnifiedDataRecord(
                    id=gene_id,
                    values=values,
                )
                records.append(record)

            statistics = {
                "total_features": len(records),
                "total_samples": len(headers) - 1 if headers else 0,
                "file_size_bytes": file_path.stat().st_size,
            }

            return TranscriptomicsData(
                metadata=metadata,
                headers=headers,
                records=records,
                statistics=statistics,
            )
        except Exception as e:
            raise ValueError(f"Failed to parse Excel file: {str(e)}")

    async def from_unified(
        self, unified_data: UnifiedData, target_format: str, output_path: Path, **kwargs
    ) -> Path:
        """Convert unified format to target transcriptomics format."""

        target_format_lower = target_format.lower()

        if target_format_lower == "csv":
            return await self._unified_to_csv(unified_data, output_path, **kwargs)
        elif target_format_lower in ["tsv", "txt"]:
            return await self._unified_to_tsv(unified_data, output_path, **kwargs)
        elif target_format_lower == "deseq2":
            return await self._unified_to_deseq2(unified_data, output_path, **kwargs)
        else:
            raise ValueError(f"Unsupported target format: {target_format}")

    async def _unified_to_csv(
        self, unified_data: UnifiedData, output_path: Path, **kwargs
    ) -> Path:
        """Convert unified format to CSV."""

        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", newline="") as f:
            if unified_data.headers:
                writer = csv.DictWriter(f, fieldnames=unified_data.headers)
                writer.writeheader()

                for record in unified_data.records:
                    writer.writerow(record.values)
            else:
                # If no headers, write simple format
                writer = csv.writer(f)
                for record in unified_data.records:
                    writer.writerow([record.id] + list(record.values.values()))

        return output_path

    async def _unified_to_tsv(
        self, unified_data: UnifiedData, output_path: Path, **kwargs
    ) -> Path:
        """Convert unified format to TSV."""

        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", newline="") as f:
            if unified_data.headers:
                writer = csv.DictWriter(
                    f, fieldnames=unified_data.headers, delimiter="\t"
                )
                writer.writeheader()

                for record in unified_data.records:
                    writer.writerow(record.values)
            else:
                writer = csv.writer(f, delimiter="\t")
                for record in unified_data.records:
                    writer.writerow([record.id] + list(record.values.values()))

        return output_path

    async def _unified_to_deseq2(
        self, unified_data: UnifiedData, output_path: Path, **kwargs
    ) -> Path:
        """Convert unified format to DESeq2 compatible format."""

        # DESeq2 expects: first column = gene names, subsequent columns = counts
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")

            # Write header
            if unified_data.headers:
                writer.writerow(unified_data.headers)

            # Write data
            for record in unified_data.records:
                row = [record.id] + [
                    record.values.get(h, 0) for h in unified_data.headers[1:]
                ]
                writer.writerow(row)

        return output_path
