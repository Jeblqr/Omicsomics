"""Genomics data format converters."""

from pathlib import Path
from typing import Dict, List, Any
import gzip
import csv

from app.converters.base import BaseConverter, ConverterFactory
from app.schemas.unified_format import (
    UnifiedData, UnifiedDataRecord, OmicsType, GenomicsData
)


@ConverterFactory.register(OmicsType.GENOMICS)
class GenomicsConverter(BaseConverter):
    """Converter for genomics data formats."""
    
    async def to_unified(
        self,
        file_path: Path,
        sample_id: str,
        source_format: str,
        **kwargs
    ) -> UnifiedData:
        """Convert genomics file to unified format."""
        
        source_format_lower = source_format.lower()
        
        if source_format_lower in ['vcf', 'vcf.gz']:
            return await self._vcf_to_unified(file_path, sample_id, source_format, **kwargs)
        elif source_format_lower in ['bed']:
            return await self._bed_to_unified(file_path, sample_id, source_format, **kwargs)
        elif source_format_lower in ['gtf', 'gff', 'gff3']:
            return await self._gtf_to_unified(file_path, sample_id, source_format, **kwargs)
        else:
            raise ValueError(f"Unsupported source format: {source_format}")
    
    async def _vcf_to_unified(
        self,
        file_path: Path,
        sample_id: str,
        source_format: str,
        **kwargs
    ) -> GenomicsData:
        """Convert VCF to unified format."""
        
        # Open file (handle gzipped)
        open_func = gzip.open if file_path.suffix == '.gz' else open
        
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get('organism'),
            reference_genome=kwargs.get('reference_genome'),
            vcf_version=None
        )
        
        headers = [
            "CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
            "FILTER", "INFO", "FORMAT", "SAMPLE"
        ]
        
        records = []
        with open_func(file_path, 'rt') as f:
            for line in f:
                if line.startswith('##'):
                    # Parse header lines for metadata
                    if line.startswith('##fileformat='):
                        metadata.custom_fields['vcf_version'] = line.strip().split('=')[1]
                    continue
                
                if line.startswith('#CHROM'):
                    # Column header line
                    continue
                
                # Parse variant line
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    record = UnifiedDataRecord(
                        id=f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}",
                        values={
                            "CHROM": fields[0],
                            "POS": int(fields[1]),
                            "ID": fields[2] if fields[2] != '.' else None,
                            "REF": fields[3],
                            "ALT": fields[4],
                            "QUAL": float(fields[5]) if fields[5] != '.' else None,
                            "FILTER": fields[6],
                            "INFO": fields[7] if len(fields) > 7 else None,
                            "FORMAT": fields[8] if len(fields) > 8 else None,
                            "SAMPLE": fields[9] if len(fields) > 9 else None
                        }
                    )
                    records.append(record)
        
        statistics = {
            "total_variants": len(records),
            "file_size_bytes": file_path.stat().st_size
        }
        
        return GenomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics=statistics
        )
    
    async def _bed_to_unified(
        self,
        file_path: Path,
        sample_id: str,
        source_format: str,
        **kwargs
    ) -> GenomicsData:
        """Convert BED to unified format."""
        
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get('organism'),
            reference_genome=kwargs.get('reference_genome')
        )
        
        headers = ["CHROM", "START", "END", "NAME", "SCORE", "STRAND"]
        records = []
        
        with open(file_path, 'r') as f:
            for idx, line in enumerate(f):
                if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    record = UnifiedDataRecord(
                        id=f"{fields[0]}:{fields[1]}-{fields[2]}",
                        values={
                            "CHROM": fields[0],
                            "START": int(fields[1]),
                            "END": int(fields[2]),
                            "NAME": fields[3] if len(fields) > 3 else f"region_{idx}",
                            "SCORE": float(fields[4]) if len(fields) > 4 else None,
                            "STRAND": fields[5] if len(fields) > 5 else None
                        }
                    )
                    records.append(record)
        
        statistics = {
            "total_regions": len(records),
            "file_size_bytes": file_path.stat().st_size
        }
        
        return GenomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics=statistics
        )
    
    async def _gtf_to_unified(
        self,
        file_path: Path,
        sample_id: str,
        source_format: str,
        **kwargs
    ) -> GenomicsData:
        """Convert GTF/GFF to unified format."""
        
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            organism=kwargs.get('organism'),
            reference_genome=kwargs.get('reference_genome')
        )
        
        headers = ["SEQNAME", "SOURCE", "FEATURE", "START", "END", "SCORE", "STRAND", "FRAME", "ATTRIBUTES"]
        records = []
        
        with open(file_path, 'r') as f:
            for idx, line in enumerate(f):
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    record = UnifiedDataRecord(
                        id=f"feature_{idx}",
                        values={
                            "SEQNAME": fields[0],
                            "SOURCE": fields[1],
                            "FEATURE": fields[2],
                            "START": int(fields[3]),
                            "END": int(fields[4]),
                            "SCORE": float(fields[5]) if fields[5] != '.' else None,
                            "STRAND": fields[6],
                            "FRAME": fields[7],
                            "ATTRIBUTES": fields[8]
                        }
                    )
                    records.append(record)
        
        statistics = {
            "total_features": len(records),
            "file_size_bytes": file_path.stat().st_size
        }
        
        return GenomicsData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics=statistics
        )
    
    async def from_unified(
        self,
        unified_data: UnifiedData,
        target_format: str,
        output_path: Path,
        **kwargs
    ) -> Path:
        """Convert unified format to target genomics format."""
        
        target_format_lower = target_format.lower()
        
        if target_format_lower in ['vcf', 'vcf.gz']:
            return await self._unified_to_vcf(unified_data, output_path, **kwargs)
        elif target_format_lower == 'bed':
            return await self._unified_to_bed(unified_data, output_path, **kwargs)
        elif target_format_lower in ['plink-bed']:  # PLINK binary format
            return await self._unified_to_plink(unified_data, output_path, **kwargs)
        else:
            raise ValueError(f"Unsupported target format: {target_format}")
    
    async def _unified_to_vcf(
        self,
        unified_data: UnifiedData,
        output_path: Path,
        **kwargs
    ) -> Path:
        """Convert unified format to VCF."""
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        open_func = gzip.open if output_path.suffix == '.gz' else open
        mode = 'wt' if output_path.suffix == '.gz' else 'w'
        
        with open_func(output_path, mode) as f:
            # Write header
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##source={unified_data.metadata.source_format}\n")
            f.write(f"##reference={unified_data.metadata.reference_genome or 'unknown'}\n")
            
            # Write column headers
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            # Write records
            for record in unified_data.records:
                values = record.values
                f.write("\t".join([
                    str(values.get("CHROM", ".")),
                    str(values.get("POS", ".")),
                    str(values.get("ID") or "."),
                    str(values.get("REF", ".")),
                    str(values.get("ALT", ".")),
                    str(values.get("QUAL") or "."),
                    str(values.get("FILTER", ".")),
                    str(values.get("INFO", "."))
                ]) + "\n")
        
        return output_path
    
    async def _unified_to_bed(
        self,
        unified_data: UnifiedData,
        output_path: Path,
        **kwargs
    ) -> Path:
        """Convert unified format to BED."""
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            for record in unified_data.records:
                values = record.values
                fields = [
                    str(values.get("CHROM", "chr1")),
                    str(values.get("START", 0)),
                    str(values.get("END", 0))
                ]
                
                if "NAME" in values:
                    fields.append(str(values["NAME"]))
                if "SCORE" in values:
                    fields.append(str(values["SCORE"]))
                if "STRAND" in values:
                    fields.append(str(values["STRAND"]))
                
                f.write("\t".join(fields) + "\n")
        
        return output_path
    
    async def _unified_to_plink(
        self,
        unified_data: UnifiedData,
        output_path: Path,
        **kwargs
    ) -> Path:
        """Convert unified format to PLINK format."""
        # This would require more complex implementation
        # For now, raise not implemented
        raise NotImplementedError("PLINK conversion not yet implemented")
