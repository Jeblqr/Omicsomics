"""Base converter class for data format conversion."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Optional
import json

from app.schemas.unified_format import (
    UnifiedData, UnifiedMetadata, UnifiedDataRecord,
    OmicsType, ProcessingStep, SUPPORTED_FORMATS, EXPORT_FORMATS
)


class BaseConverter(ABC):
    """Base class for all format converters."""
    
    def __init__(self, omics_type: OmicsType):
        """Initialize converter with omics type."""
        self.omics_type = omics_type
        self.supported_inputs = SUPPORTED_FORMATS.get(omics_type, [])
    
    def is_supported_format(self, file_format: str) -> bool:
        """Check if file format is supported."""
        return file_format.lower() in self.supported_inputs
    
    @abstractmethod
    async def to_unified(
        self, 
        file_path: Path,
        sample_id: str,
        source_format: str,
        **kwargs
    ) -> UnifiedData:
        """Convert from source format to unified format."""
        pass
    
    @abstractmethod
    async def from_unified(
        self,
        unified_data: UnifiedData,
        target_format: str,
        output_path: Path,
        **kwargs
    ) -> Path:
        """Convert from unified format to target format."""
        pass
    
    def create_metadata(
        self,
        sample_id: str,
        source_format: str,
        organism: Optional[str] = None,
        reference_genome: Optional[str] = None,
        **custom_fields
    ) -> UnifiedMetadata:
        """Create metadata object."""
        return UnifiedMetadata(
            omics_type=self.omics_type,
            source_format=source_format,
            sample_id=sample_id,
            organism=organism,
            reference_genome=reference_genome,
            custom_fields=custom_fields
        )
    
    def add_processing_step(
        self,
        unified_data: UnifiedData,
        step_name: str,
        tool: str,
        version: Optional[str] = None,
        **parameters
    ) -> UnifiedData:
        """Add a processing step to the metadata."""
        step = ProcessingStep(
            step_name=step_name,
            tool=tool,
            version=version,
            parameters=parameters
        )
        unified_data.metadata.processing_steps.append(step)
        return unified_data
    
    async def save_unified(
        self,
        unified_data: UnifiedData,
        output_path: Path
    ) -> Path:
        """Save unified data to JSON file."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(
                unified_data.model_dump(mode='json'),
                f,
                indent=2,
                default=str
            )
        
        return output_path
    
    async def load_unified(self, file_path: Path) -> UnifiedData:
        """Load unified data from JSON file."""
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        return UnifiedData(**data)


class ConverterFactory:
    """Factory for creating appropriate converters."""
    
    _converters: Dict[OmicsType, type] = {}
    
    @classmethod
    def register(cls, omics_type: OmicsType):
        """Decorator to register a converter class."""
        def decorator(converter_class):
            cls._converters[omics_type] = converter_class
            return converter_class
        return decorator
    
    @classmethod
    def create(cls, omics_type: OmicsType) -> BaseConverter:
        """Create a converter instance for the given omics type."""
        converter_class = cls._converters.get(omics_type)
        if not converter_class:
            raise ValueError(f"No converter registered for {omics_type}")
        return converter_class(omics_type)
    
    @classmethod
    def get_supported_formats(cls, omics_type: OmicsType) -> list:
        """Get supported input formats for an omics type."""
        return SUPPORTED_FORMATS.get(omics_type, [])
    
    @classmethod
    def get_export_formats(cls, tool: str) -> list:
        """Get supported export formats for a tool."""
        return EXPORT_FORMATS.get(tool, [])
