"""
Metabolomics data converter.

Handles conversion between metabolomics file formats and unified format.
Supports: mzData, mzML, netCDF, CSV, TSV
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
    MetabolomicsData,
)


@ConverterFactory.register(OmicsType.METABOLOMICS)
class MetabolomicsConverter(BaseConverter):
    """Converter for metabolomics data formats."""

    async def to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> UnifiedData:
        """
        Convert metabolomics data to unified format.

        Args:
            file_path: Path to the metabolomics data file
            sample_id: Sample identifier
            source_format: Format of the input file (mzData, mzML, CDF, CSV, TSV)
            **kwargs: Additional parameters (organism, reference_genome, etc.)

        Returns:
            MetabolomicsData object
        """
        format_lower = source_format.lower()

        if format_lower in ["mzdata"]:
            return await self._mzdata_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif format_lower in ["mzml"]:
            return await self._mzml_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif format_lower in ["cdf", "netcdf"]:
            return await self._cdf_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif format_lower in ["csv", "tsv", "txt"]:
            return await self._tabular_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        else:
            raise ValueError(f"Unsupported metabolomics format: {source_format}")

    async def _mzdata_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> MetabolomicsData:
        """
        Convert mzData to unified format.
        mzData is an older XML format for mass spectrometry data.
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
            ns["mzData"] = root.tag.split("}")[0][1:]

        headers = [
            "spectrum_id",
            "ms_level",
            "retention_time",
            "mz_range_start",
            "mz_range_stop",
            "peak_count",
        ]
        records = []

        # Extract spectra
        spectrum_list = (
            root.find(".//spectrumList", ns)
            if not ns
            else root.find(f'.//{{{ns["mzData"]}}}spectrumList', ns)
        )

        if spectrum_list is not None:
            for spectrum in spectrum_list.findall(
                "spectrum" if not ns else f'{{{ns["mzData"]}}}spectrum'
            ):
                spectrum_id = spectrum.get("id", "")

                # Extract spectrum description
                ms_level = ""
                rt = ""
                mz_start = ""
                mz_stop = ""

                spec_desc = spectrum.find(
                    "spectrumDesc" if not ns else f'{{{ns["mzData"]}}}spectrumDesc'
                )
                if spec_desc is not None:
                    # Extract CV params
                    for cv_param in spec_desc.iter(
                        "cvParam" if not ns else f'{{{ns["mzData"]}}}cvParam'
                    ):
                        name = cv_param.get("name", "")
                        value = cv_param.get("value", "")

                        if "ms level" in name.lower():
                            ms_level = value
                        elif "time" in name.lower() and "retention" in name.lower():
                            rt = value

                    # Extract m/z range
                    mz_range = spec_desc.find(
                        ".//mzRangeStart"
                        if not ns
                        else f'.//{{{ns["mzData"]}}}mzRangeStart'
                    )
                    if mz_range is not None:
                        mz_start = mz_range.get("value", "")

                    mz_range = spec_desc.find(
                        ".//mzRangeStop"
                        if not ns
                        else f'.//{{{ns["mzData"]}}}mzRangeStop'
                    )
                    if mz_range is not None:
                        mz_stop = mz_range.get("value", "")

                # Count peaks
                peak_count = 0
                mz_array_length = spectrum.find(
                    ".//mzArrayBinary"
                    if not ns
                    else f'.//{{{ns["mzData"]}}}mzArrayBinary'
                )
                if mz_array_length is not None:
                    length_attr = mz_array_length.find(
                        ".//data" if not ns else f'{{{ns["mzData"]}}}data'
                    )
                    if length_attr is not None:
                        peak_count = length_attr.get("length", "0")

                record = UnifiedDataRecord(
                    id=spectrum_id or f"spectrum_{len(records)}",
                    values={
                        "spectrum_id": spectrum_id,
                        "ms_level": ms_level,
                        "retention_time": rt,
                        "mz_range_start": mz_start,
                        "mz_range_stop": mz_stop,
                        "peak_count": peak_count,
                    },
                )
                records.append(record)

        metadata.custom_fields["total_spectra"] = len(records)

        return MetabolomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={"total_spectra": len(records)},
        )

    async def _mzml_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> MetabolomicsData:
        """
        Convert mzML to unified format (metabolomics context).
        Similar to proteomics but may have different metadata emphasis.
        """
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )
        metadata.custom_fields["data_type"] = "metabolomics"

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

        return MetabolomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={"total_spectra": len(records)},
        )

    async def _cdf_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> MetabolomicsData:
        """
        Convert netCDF to unified format.
        Note: This requires netCDF4 library. Providing a placeholder implementation.
        """
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )
        metadata.custom_fields["note"] = "Full netCDF parsing requires netCDF4 library"

        headers = ["file_info"]
        records = [
            UnifiedDataRecord(
                id="netcdf_file", values={"file_info": f"netCDF file: {file_path.name}"}
            )
        ]

        return MetabolomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={"requires_netcdf4": True},
        )

    async def _tabular_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> MetabolomicsData:
        """
        Convert tabular metabolomics data (CSV/TSV) to unified format.
        Common for metabolite identification results and quantification tables.
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
                        row.get(
                            "metabolite_id", row.get("compound_id", f"record_{idx}")
                        ),
                    ),
                    values=dict(row),
                )
                records.append(record)

        metadata.custom_fields["total_metabolites"] = len(records)

        return MetabolomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={"total_records": len(records)},
        )

    async def from_unified(
        self, unified_data: UnifiedData, target_format: str, output_path: Path, **kwargs
    ) -> Path:
        """
        Convert unified format to metabolomics data format.

        Args:
            unified_data: UnifiedData object to convert
            target_format: Desired output format (CSV, TSV)
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
