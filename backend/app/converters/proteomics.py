"""
Proteomics data converter.

Handles conversion between proteomics file formats and unified format.
Supports: mzML, mzXML, MGF (Mascot Generic Format), CSV, TSV, Thermo RAW
"""

import csv
import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Any, Optional

from app.converters.base import BaseConverter, ConverterFactory
from app.schemas.unified_format import (
    UnifiedData,
    UnifiedDataRecord,
    OmicsType,
    ProteomicsData,
)

logger = logging.getLogger(__name__)

# Optional pyteomics import for advanced RAW file parsing
try:
    from pyteomics import mzml as pymzml
    from pyteomics import mzxml as pymzxml

    PYTEOMICS_AVAILABLE = True
except ImportError:
    PYTEOMICS_AVAILABLE = False
    logger.warning(
        "pyteomics not available. Advanced binary format parsing will be limited. "
        "Install with: pip install pyteomics"
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
            source_format: Format of the input file (mzML, mzXML, MGF, CSV, TSV, RAW)
            **kwargs: Additional parameters (organism, reference_genome, etc.)

        Returns:
            ProteomicsData object
        """
        format_lower = source_format.lower()

        if format_lower in ["mzml"]:
            # Try pyteomics first for better parsing
            if PYTEOMICS_AVAILABLE:
                return await self._mzml_to_unified_pyteomics(
                    file_path, sample_id, source_format, **kwargs
                )
            else:
                return await self._mzml_to_unified(
                    file_path, sample_id, source_format, **kwargs
                )
        elif format_lower in ["mzxml"]:
            if PYTEOMICS_AVAILABLE:
                return await self._mzxml_to_unified_pyteomics(
                    file_path, sample_id, source_format, **kwargs
                )
            else:
                return await self._mzxml_to_unified(
                    file_path, sample_id, source_format, **kwargs
                )
        elif format_lower in ["mgf"]:
            return await self._mgf_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif format_lower in ["raw"]:
            # Thermo RAW files require ThermoRawFileParser conversion first
            return await self._raw_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        elif format_lower in ["csv", "tsv", "txt"]:
            return await self._tabular_to_unified(
                file_path, sample_id, source_format, **kwargs
            )
        else:
            raise ValueError(f"Unsupported proteomics format: {source_format}")

    async def _raw_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> ProteomicsData:
        """
        Convert Thermo RAW to unified format.
        
        This method provides metadata extraction from RAW files.
        For full spectrum data, ThermoRawFileParser should be used to convert
        RAW -> mzML first, then process with _mzml_to_unified.
        
        Note: Direct RAW file parsing requires platform-specific libraries.
        For production use, convert RAW to mzML using ThermoRawFileParser:
            ThermoRawFileParser -i input.raw -o output_dir -f 2
        """
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )
        
        # Extract basic file metadata
        file_stats = file_path.stat()
        metadata.custom_fields.update({
            "file_size_mb": round(file_stats.st_size / (1024 * 1024), 2),
            "file_name": file_path.name,
            "parsing_note": (
                "RAW file metadata only. For full spectrum data, "
                "convert to mzML using ThermoRawFileParser first."
            ),
            "conversion_command": (
                f"ThermoRawFileParser -i {file_path.name} -o output_dir -f 2"
            ),
        })
        
        headers = ["property", "value"]
        records = [
            UnifiedDataRecord(
                id="file_info",
                values={
                    "property": "File Name",
                    "value": file_path.name,
                },
            ),
            UnifiedDataRecord(
                id="file_size",
                values={
                    "property": "File Size (MB)",
                    "value": str(metadata.custom_fields["file_size_mb"]),
                },
            ),
            UnifiedDataRecord(
                id="format",
                values={
                    "property": "Format",
                    "value": "Thermo RAW",
                },
            ),
            UnifiedDataRecord(
                id="conversion_required",
                values={
                    "property": "Conversion Required",
                    "value": "Yes - Use ThermoRawFileParser to convert to mzML",
                },
            ),
        ]
        
        return ProteomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={"metadata_records": len(records)},
        )

    async def _mzml_to_unified_pyteomics(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> ProteomicsData:
        """
        Convert mzML to unified format using pyteomics for better parsing.
        
        Pyteomics provides more robust parsing of complex mzML files with
        better handling of different vendor formats and edge cases.
        """
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )
        
        headers = [
            "spectrum_id",
            "ms_level",
            "scan_number",
            "retention_time_seconds",
            "precursor_mz",
            "precursor_charge",
            "num_peaks",
            "base_peak_mz",
            "base_peak_intensity",
            "total_ion_current",
        ]
        records = []
        
        try:
            # Use pyteomics for robust parsing
            with pymzml.read(str(file_path)) as reader:
                for spectrum in reader:
                    spectrum_id = spectrum.get("id", f"scan_{spectrum.get('index', len(records))}")
                    ms_level = spectrum.get("ms level", 1)
                    scan_num = spectrum.get("index", len(records))
                    
                    # Extract retention time (in seconds)
                    rt = spectrum.get("scanList", {}).get("scan", [{}])[0].get(
                        "scan start time", None
                    )
                    
                    # Extract precursor info for MS2/MS3
                    precursor_mz = None
                    precursor_charge = None
                    if ms_level > 1:
                        precursors = spectrum.get("precursorList", {}).get("precursor", [])
                        if precursors:
                            selected_ion = precursors[0].get("selectedIonList", {}).get(
                                "selectedIon", [{}]
                            )[0]
                            precursor_mz = selected_ion.get("selected ion m/z")
                            precursor_charge = selected_ion.get("charge state")
                    
                    # Extract peak data
                    mz_array = spectrum.get("m/z array", [])
                    intensity_array = spectrum.get("intensity array", [])
                    num_peaks = len(mz_array)
                    
                    # Calculate base peak
                    base_peak_mz = None
                    base_peak_intensity = None
                    if num_peaks > 0:
                        max_idx = intensity_array.argmax() if len(intensity_array) > 0 else 0
                        base_peak_mz = float(mz_array[max_idx]) if max_idx < len(mz_array) else None
                        base_peak_intensity = float(intensity_array[max_idx])
                    
                    # Total ion current
                    tic = spectrum.get("total ion current", 
                                      intensity_array.sum() if len(intensity_array) > 0 else 0)
                    
                    record = UnifiedDataRecord(
                        id=spectrum_id,
                        values={
                            "spectrum_id": spectrum_id,
                            "ms_level": str(ms_level),
                            "scan_number": str(scan_num),
                            "retention_time_seconds": str(rt) if rt else "",
                            "precursor_mz": str(precursor_mz) if precursor_mz else "",
                            "precursor_charge": str(precursor_charge) if precursor_charge else "",
                            "num_peaks": str(num_peaks),
                            "base_peak_mz": str(base_peak_mz) if base_peak_mz else "",
                            "base_peak_intensity": str(base_peak_intensity) if base_peak_intensity else "",
                            "total_ion_current": str(tic),
                        },
                    )
                    records.append(record)
                    
        except Exception as e:
            logger.error(f"Error parsing mzML with pyteomics: {e}")
            # Fall back to basic XML parsing
            return await self._mzml_to_unified(file_path, sample_id, source_format, **kwargs)
        
        metadata.custom_fields.update({
            "total_spectra": len(records),
            "parser": "pyteomics",
        })
        
        return ProteomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={
                "total_spectra": len(records),
                "ms1_spectra": sum(1 for r in records if r.values.get("ms_level") == "1"),
                "ms2_spectra": sum(1 for r in records if r.values.get("ms_level") == "2"),
            },
        )

    async def _mzxml_to_unified_pyteomics(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> ProteomicsData:
        """
        Convert mzXML to unified format using pyteomics.
        
        Similar to mzML parsing but for mzXML format.
        """
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get("organism"),
            reference_genome=kwargs.get("reference_genome"),
        )
        
        headers = [
            "scan_number",
            "ms_level",
            "retention_time_seconds",
            "precursor_mz",
            "num_peaks",
            "base_peak_mz",
            "base_peak_intensity",
            "total_ion_current",
        ]
        records = []
        
        try:
            with pymzxml.read(str(file_path)) as reader:
                for scan in reader:
                    scan_num = scan.get("num", len(records))
                    ms_level = scan.get("msLevel", 1)
                    rt = scan.get("retentionTime", "")  # in seconds
                    
                    # Precursor info
                    precursor_mz = None
                    if "precursorMz" in scan:
                        precursors = scan["precursorMz"]
                        if precursors and len(precursors) > 0:
                            precursor_mz = precursors[0].get("precursorMz")
                    
                    # Peak data
                    mz_array = scan.get("m/z array", [])
                    intensity_array = scan.get("intensity array", [])
                    num_peaks = len(mz_array)
                    
                    # Base peak
                    base_peak_mz = None
                    base_peak_intensity = None
                    if num_peaks > 0:
                        max_idx = intensity_array.argmax() if len(intensity_array) > 0 else 0
                        base_peak_mz = float(mz_array[max_idx]) if max_idx < len(mz_array) else None
                        base_peak_intensity = float(intensity_array[max_idx])
                    
                    # TIC
                    tic = scan.get("totIonCurrent", 
                                  intensity_array.sum() if len(intensity_array) > 0 else 0)
                    
                    record = UnifiedDataRecord(
                        id=f"scan_{scan_num}",
                        values={
                            "scan_number": str(scan_num),
                            "ms_level": str(ms_level),
                            "retention_time_seconds": str(rt) if rt else "",
                            "precursor_mz": str(precursor_mz) if precursor_mz else "",
                            "num_peaks": str(num_peaks),
                            "base_peak_mz": str(base_peak_mz) if base_peak_mz else "",
                            "base_peak_intensity": str(base_peak_intensity) if base_peak_intensity else "",
                            "total_ion_current": str(tic),
                        },
                    )
                    records.append(record)
                    
        except Exception as e:
            logger.error(f"Error parsing mzXML with pyteomics: {e}")
            # Fall back to basic XML parsing
            return await self._mzxml_to_unified(file_path, sample_id, source_format, **kwargs)
        
        metadata.custom_fields.update({
            "total_scans": len(records),
            "parser": "pyteomics",
        })
        
        return ProteomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics={
                "total_scans": len(records),
                "ms1_scans": sum(1 for r in records if r.values.get("ms_level") == "1"),
                "ms2_scans": sum(1 for r in records if r.values.get("ms_level") == "2"),
            },
        )

    async def _mzml_to_unified(
        self, file_path: Path, sample_id: str, source_format: str, **kwargs
    ) -> ProteomicsData:
        """
        Convert mzML to unified format using basic XML parsing.
        Extracts spectrum information (m/z, intensity, retention time).
        
        This is a fallback method when pyteomics is not available.
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
        metadata.custom_fields["parser"] = "xml_fallback"

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
