"""
Proteomics data converter.

Handles conversion between proteomics file formats and unified format.
Supports: mzML, mzXML, MGF (Mascot Generic Format), CSV, TSV
"""

import csv
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Any

from app.converters.base import BaseConverter, ConverterFactory
from app.schemas.unified_format import (
    UnifiedData,
    UnifiedDataRecord,
    OmicsType,
    ProteomicsData,
)


@ConverterFactory.register(OmicsType.PROTEOMICS)
class ProteomicsConverter(BaseConverter):
    """Converter for proteomics data formats."""

    async def to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> UnifiedData:
        """
        Convert proteomics data to unified format.

        Args:
            file_path: Path to the proteomics data file
            sample_id: Sample identifier
            source_format: Format of the input file (mzML, mzXML, MGF, CSV, TSV)
            **kwargs: Additional parameters (organism, reference_genome, etc.)

        Returns:
            ProteomicsData object
        """
        format_lower = source_format.lower()

        if format_lower in ["mzml"]:
            return await self._mzml_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif format_lower in ["mzxml"]:
            return await self._mzxml_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif format_lower in ["mgf"]:
            return await self._mgf_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif format_lower in ["csv", "tsv", "txt"]:
            return await self._tabular_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        else:
            raise ValueError(f"Unsupported proteomics format: {source_format}")

    async def _mzml_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> ProteomicsData:
        """
        Convert mzML to unified format.
        Extracts spectrum information (m/z, intensity, retention time).
        """
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )

        tree = ET.parse(file_path)
        root = tree.getroot()

        # Extract namespace
        ns = {"mzML": "http://psi.hupo.org/ms/mzml"}
        if root.tag.startswith("{"):
            ns["mzML"] = root.tag.split("}")[0][1:]

        headers = [
            "spectrum_id",
            "ms_level",
            "retention_time",
            "base_peak_mz",
            "base_peak_intensity",
            "total_ion_current",
        ]
        records = []

        # Extract spectra
        for spectrum in root.findall(".//mzML:spectrum", ns):
            spectrum_id = spectrum.get("id", "")
            ms_level = spectrum.get("msLevel", "1")

            # Extract retention time and other CV params
            rt = None
            base_peak_mz = None
            base_peak_intensity = None
            tic = None

            for cv_param in spectrum.findall(".//mzML:cvParam", ns):
                name = cv_param.get("name", "")
                value = cv_param.get("value", "")

                if "scan start time" in name.lower():
                    rt = value
                elif "base peak m/z" in name.lower():
                    base_peak_mz = value
                elif "base peak intensity" in name.lower():
                    base_peak_intensity = value
                elif "total ion current" in name.lower():
                    tic = value

            record = UnifiedDataRecord(
                id=spectrum_id or f"spectrum_{len(records)}",
                values={
                    "spectrum_id": spectrum_id,
                    "ms_level": ms_level,
                    "retention_time": rt or "",
                    "base_peak_mz": base_peak_mz or "",
                    "base_peak_intensity": base_peak_intensity or "",
                    "total_ion_current": tic or "",
                },
            )
            records.append(record)

        metadata.custom_fields["total_spectra"] = len(records)

        return ProteomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={"total_spectra": len(records)},
        )

    async def _mzxml_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> ProteomicsData:
        """
        Convert mzXML to unified format.
        Similar to mzML but with different XML schema.
        """
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )

        tree = ET.parse(file_path)
        root = tree.getroot()

        # Extract namespace
        ns = {}
        if root.tag.startswith("{"):
            ns["mzXML"] = root.tag.split("}")[0][1:]

        headers = [
            "scan_num",
            "ms_level",
            "retention_time",
            "base_peak_mz",
            "base_peak_intensity",
            "total_ion_current",
        ]
        records = []

        # Extract scans
        for scan in (
            root.iter("scan") if not ns else root.iter(f'{{{ns["mzXML"]}}}scan')
        ):
            scan_num = scan.get("num", "")
            ms_level = scan.get("msLevel", "1")
            rt = scan.get("retentionTime", "")
            base_peak_mz = scan.get("basePeakMz", "")
            base_peak_intensity = scan.get("basePeakIntensity", "")
            tic = scan.get("totIonCurrent", "")

            record = UnifiedDataRecord(
                id=f"scan_{scan_num}",
                values={
                    "scan_num": scan_num,
                    "ms_level": ms_level,
                    "retention_time": rt,
                    "base_peak_mz": base_peak_mz,
                    "base_peak_intensity": base_peak_intensity,
                    "total_ion_current": tic,
                },
            )
            records.append(record)

        metadata.custom_fields["total_scans"] = len(records)

        return ProteomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={"total_scans": len(records)},
        )

    async def _mgf_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> ProteomicsData:
        """
        Convert MGF (Mascot Generic Format) to unified format.
        MGF is a text format for MS/MS spectra.
        """
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )

        headers = ["title", "rtinseconds", "pepmass", "charge", "num_peaks"]
        records = []

        with open(file_path, "r") as f:
            current_spectrum = {}
            peak_count = 0
            spectrum_id = 0

            for line in f:
                line = line.strip()

                if line == "BEGIN IONS":
                    current_spectrum = {}
                    peak_count = 0
                    spectrum_id += 1
                elif line == "END IONS":
                    current_spectrum["num_peaks"] = str(peak_count)
                    record = UnifiedDataRecord(
                        id=current_spectrum.get("title", f"spectrum_{spectrum_id}"),
                        values=current_spectrum,
                    )
                    records.append(record)
                elif "=" in line:
                    key, value = line.split("=", 1)
                    key_lower = key.lower()
                    if key_lower in ["title", "rtinseconds", "pepmass", "charge"]:
                        current_spectrum[key_lower] = value
                elif line and not line.startswith("#"):
                    # Peak data (m/z intensity pairs)
                    peak_count += 1

        metadata.custom_fields["total_spectra"] = len(records)

        return ProteomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={"total_spectra": len(records)},
        )

    async def _tabular_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> ProteomicsData:
        """
        Convert tabular proteomics data (CSV/TSV) to unified format.
        Common for protein identification results (e.g., from search engines).
        """
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )

        delimiter = "," if source_format.lower() == "csv" else "\t"

        with open(file_path, "r") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            headers = list(reader.fieldnames or [])
            records = []

            for idx, row in enumerate(reader):
                record = UnifiedDataRecord(
                    id=row.get(
                        "id",
                        row.get("protein_id", row.get("peptide_id", f"record_{idx}")),
                    ),
                    values=dict(row),
                )
                records.append(record)

        metadata.custom_fields["total_proteins"] = len(records)

        return ProteomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={"total_records": len(records)},
        )

    async def from_unified(
        self, unified_data: UnifiedData, target_format: str, output_path: Path, **kwargs
    ) -> Path:
        """
        Convert unified format to proteomics data format.

        Args:
            unified_data: UnifiedData object to convert
            target_format: Desired output format (CSV, TSV, MGF)
            output_path: Path where to save the output file
            **kwargs: Additional parameters

        Returns:
            Path to the converted file
        """
        format_lower = target_format.lower()

        if format_lower == "csv":
            await self._unified_to_csv(unified_data, output_path)
        elif format_lower == "tsv":
            await self._unified_to_tsv(unified_data, output_path)
        elif format_lower == "mgf":
            await self._unified_to_mgf(unified_data, output_path)
        else:
            raise ValueError(f"Unsupported output format: {target_format}")

        return output_path

    async def _unified_to_csv(
        self, unified_data: UnifiedData, output_path: Path
    ) -> None:
        """Convert unified format to CSV."""
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=unified_data.headers)
            writer.writeheader()
            for record in unified_data.records:
                writer.writerow(record.values)

    async def _unified_to_tsv(
        self, unified_data: UnifiedData, output_path: Path
    ) -> None:
        """Convert unified format to TSV."""
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=unified_data.headers, delimiter="\t")
            writer.writeheader()
            for record in unified_data.records:
                writer.writerow(record.values)

    async def _unified_to_mgf(
        self, unified_data: UnifiedData, output_path: Path
    ) -> None:
        """
        Convert unified format to MGF.
        Note: This is a simplified conversion - full spectrum data would need to be stored.
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            for record in unified_data.records:
                f.write("BEGIN IONS\n")
                for key, value in record.values.items():
                    if value and key != "num_peaks":
                        f.write(f"{key.upper()}={value}\n")
                f.write("END IONS\n\n")
