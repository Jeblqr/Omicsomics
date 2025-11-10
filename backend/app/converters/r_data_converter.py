"""
R data format converter for RData files.

Supports conversions:
- RData -> CSV (data frames, matrices)
- RData -> JSON (lists, vectors)

Dependencies:
- rpy2: R interface for Python
"""
from pathlib import Path
from typing import Dict, Optional, List
import json


class RDataConverter:
    """
    Converter for R data formats (.RData, .rda, .rds).
    
    Singleton pattern for efficiency.
    """
    
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def __init__(self):
        """Initialize R data converter."""
        if not hasattr(self, '_initialized'):
            self._initialized = True
            self._check_dependencies()
    
    def _check_dependencies(self):
        """Check if required tools are available."""
        try:
            import rpy2.robjects as robjects
            from rpy2.robjects import pandas2ri
            self.robjects = robjects
            self.pandas2ri = pandas2ri
            self.has_rpy2 = True
        except ImportError:
            self.has_rpy2 = False
        
        try:
            import pandas as pd
            self.pandas = pd
            self.has_pandas = True
        except ImportError:
            self.has_pandas = False
    
    # ============================================================
    # RData Conversions
    # ============================================================
    
    def convert_rdata_to_csv(self, source: str, target: str,
                            object_name: Optional[str] = None):
        """
        Convert RData file to CSV format.
        
        Only works for data frames and matrices.
        
        Args:
            source: Input RData file path
            target: Output CSV file path
            object_name: Name of R object to extract (first suitable if not specified)
        """
        if not self.has_rpy2:
            raise ImportError("rpy2 is required for RData conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for RData conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Activate pandas conversion
        self.pandas2ri.activate()
        
        try:
            # Load RData file
            self.robjects.r['load'](str(source_path))
            
            # Get object names
            object_names = list(self.robjects.globalenv.keys())
            
            if not object_names:
                raise ValueError("No objects found in RData file")
            
            # Select object
            if object_name:
                if object_name not in object_names:
                    raise ValueError(f"Object '{object_name}' not found. Available: {object_names}")
                r_obj = self.robjects.globalenv[object_name]
            else:
                # Use first data frame or matrix
                r_obj = None
                for name in object_names:
                    obj = self.robjects.globalenv[name]
                    if self._is_dataframe(obj) or self._is_matrix(obj):
                        r_obj = obj
                        break
                
                if r_obj is None:
                    raise ValueError("No data frame or matrix found in RData file")
            
            # Convert to pandas DataFrame
            if self._is_dataframe(r_obj):
                df = self.pandas2ri.rpy2py(r_obj)
            elif self._is_matrix(r_obj):
                # Convert matrix to DataFrame
                matrix_data = self.pandas2ri.rpy2py(r_obj)
                df = self.pandas.DataFrame(matrix_data)
            else:
                raise ValueError("Object is not a data frame or matrix")
            
            # Write CSV
            df.to_csv(target_path, index=False)
            
        finally:
            self.pandas2ri.deactivate()
    
    def convert_rdata_to_json(self, source: str, target: str,
                             object_name: Optional[str] = None):
        """
        Convert RData file to JSON format.
        
        Works for vectors, lists, and simple objects.
        
        Args:
            source: Input RData file path
            target: Output JSON file path
            object_name: Name of R object to extract (first if not specified)
        """
        if not self.has_rpy2:
            raise ImportError("rpy2 is required for RData conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Load RData file
        self.robjects.r['load'](str(source_path))
        
        # Get object names
        object_names = list(self.robjects.globalenv.keys())
        
        if not object_names:
            raise ValueError("No objects found in RData file")
        
        # Select object
        if object_name:
            if object_name not in object_names:
                raise ValueError(f"Object '{object_name}' not found. Available: {object_names}")
            r_obj = self.robjects.globalenv[object_name]
        else:
            r_obj = self.robjects.globalenv[object_names[0]]
        
        # Convert to Python object
        py_obj = self._r_to_python(r_obj)
        
        # Write JSON
        with open(target_path, 'w') as f:
            json.dump(py_obj, f, indent=2)
    
    def convert_rds_to_csv(self, source: str, target: str):
        """
        Convert RDS file to CSV format.
        
        Only works for data frames and matrices.
        
        Args:
            source: Input RDS file path
            target: Output CSV file path
        """
        if not self.has_rpy2:
            raise ImportError("rpy2 is required for RDS conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for RDS conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Activate pandas conversion
        self.pandas2ri.activate()
        
        try:
            # Load RDS file
            r_obj = self.robjects.r['readRDS'](str(source_path))
            
            # Convert to pandas DataFrame
            if self._is_dataframe(r_obj):
                df = self.pandas2ri.rpy2py(r_obj)
            elif self._is_matrix(r_obj):
                matrix_data = self.pandas2ri.rpy2py(r_obj)
                df = self.pandas.DataFrame(matrix_data)
            else:
                raise ValueError("RDS object is not a data frame or matrix")
            
            # Write CSV
            df.to_csv(target_path, index=False)
            
        finally:
            self.pandas2ri.deactivate()
    
    # ============================================================
    # Helper Methods
    # ============================================================
    
    def _is_dataframe(self, r_obj) -> bool:
        """Check if R object is a data frame."""
        if not self.has_rpy2:
            return False
        
        return 'data.frame' in self.robjects.r['class'](r_obj)
    
    def _is_matrix(self, r_obj) -> bool:
        """Check if R object is a matrix."""
        if not self.has_rpy2:
            return False
        
        return 'matrix' in self.robjects.r['class'](r_obj)
    
    def _r_to_python(self, r_obj):
        """Convert R object to Python object."""
        if not self.has_rpy2:
            return None
        
        r_class = self.robjects.r['class'](r_obj)[0]
        
        if r_class in ['numeric', 'integer']:
            # Vector
            return list(r_obj)
        elif r_class == 'character':
            # Character vector
            return list(r_obj)
        elif r_class == 'logical':
            # Logical vector
            return [bool(x) for x in r_obj]
        elif r_class == 'list':
            # List
            result = {}
            names = self.robjects.r['names'](r_obj)
            if names:
                for i, name in enumerate(names):
                    result[name] = self._r_to_python(r_obj[i])
            else:
                result = [self._r_to_python(item) for item in r_obj]
            return result
        else:
            # Try to convert to string
            return str(r_obj)
    
    # ============================================================
    # Statistics
    # ============================================================
    
    def get_rdata_stats(self, file_path: str) -> Dict:
        """
        Get statistics from RData file.
        
        Args:
            file_path: RData file path
            
        Returns:
            Dictionary with statistics
        """
        if not self.has_rpy2:
            raise ImportError("rpy2 is required for RData statistics")
        
        source_path = Path(file_path)
        
        # Load RData file
        self.robjects.r['load'](str(source_path))
        
        stats = {
            'num_objects': 0,
            'objects': []
        }
        
        # Get object names
        object_names = list(self.robjects.globalenv.keys())
        stats['num_objects'] = len(object_names)
        
        # Get info for each object
        for name in object_names:
            r_obj = self.robjects.globalenv[name]
            r_class = list(self.robjects.r['class'](r_obj))
            
            obj_info = {
                'name': name,
                'class': r_class
            }
            
            # Add dimensions if applicable
            if 'data.frame' in r_class or 'matrix' in r_class:
                dims = self.robjects.r['dim'](r_obj)
                obj_info['dimensions'] = {
                    'rows': dims[0],
                    'columns': dims[1]
                }
            
            # Add length for vectors
            if any(c in r_class for c in ['numeric', 'integer', 'character', 'logical']):
                obj_info['length'] = self.robjects.r['length'](r_obj)[0]
            
            stats['objects'].append(obj_info)
        
        return stats
    
    def get_rds_stats(self, file_path: str) -> Dict:
        """
        Get statistics from RDS file.
        
        Args:
            file_path: RDS file path
            
        Returns:
            Dictionary with statistics
        """
        if not self.has_rpy2:
            raise ImportError("rpy2 is required for RDS statistics")
        
        source_path = Path(file_path)
        
        # Load RDS file
        r_obj = self.robjects.r['readRDS'](str(source_path))
        
        r_class = list(self.robjects.r['class'](r_obj))
        
        stats = {
            'class': r_class
        }
        
        # Add dimensions if applicable
        if 'data.frame' in r_class or 'matrix' in r_class:
            dims = self.robjects.r['dim'](r_obj)
            stats['dimensions'] = {
                'rows': dims[0],
                'columns': dims[1]
            }
        
        # Add length for vectors
        if any(c in r_class for c in ['numeric', 'integer', 'character', 'logical']):
            stats['length'] = self.robjects.r['length'](r_obj)[0]
        
        return stats
    
    def list_rdata_objects(self, file_path: str) -> List[str]:
        """
        List all objects in RData file.
        
        Args:
            file_path: RData file path
            
        Returns:
            List of object names
        """
        if not self.has_rpy2:
            raise ImportError("rpy2 is required for RData operations")
        
        source_path = Path(file_path)
        
        # Load RData file
        self.robjects.r['load'](str(source_path))
        
        # Get object names
        return list(self.robjects.globalenv.keys())


def get_r_data_converter() -> RDataConverter:
    """Get singleton instance of RDataConverter."""
    return RDataConverter()
