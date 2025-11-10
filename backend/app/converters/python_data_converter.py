"""
Python data format converter for HDF5, NumPy arrays.

Supports conversions between:
- h5 (HDF5 generic format)
- npy (NumPy array)
- npz (NumPy compressed archive)
- CSV (tabular data)

Dependencies:
- h5py: HDF5 file handling
- numpy: Array operations
- pandas: CSV handling
"""
from pathlib import Path
from typing import Dict, Optional, List
import csv


class PythonDataConverter:
    """
    Converter for Python data formats (HDF5, NumPy).
    
    Singleton pattern for efficiency.
    """
    
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def __init__(self):
        """Initialize Python data converter."""
        if not hasattr(self, '_initialized'):
            self._initialized = True
            self._check_dependencies()
    
    def _check_dependencies(self):
        """Check if required tools are available."""
        try:
            import h5py
            self.h5py = h5py
            self.has_h5py = True
        except ImportError:
            self.has_h5py = False
        
        try:
            import numpy as np
            self.numpy = np
            self.has_numpy = True
        except ImportError:
            self.has_numpy = False
        
        try:
            import pandas as pd
            self.pandas = pd
            self.has_pandas = True
        except ImportError:
            self.has_pandas = False
    
    # ============================================================
    # HDF5 (h5) Conversions
    # ============================================================
    
    def convert_h5_to_csv(self, source: str, target: str,
                         dataset_path: str = 'data'):
        """
        Convert HDF5 dataset to CSV format.
        
        Args:
            source: Input HDF5 file path
            target: Output CSV file path
            dataset_path: Path to dataset within HDF5 file
        """
        if not self.has_h5py:
            raise ImportError("h5py is required for HDF5 conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for HDF5 conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        with self.h5py.File(source_path, 'r') as f:
            # Get dataset
            if dataset_path not in f:
                raise ValueError(f"Dataset '{dataset_path}' not found in HDF5 file")
            
            dataset = f[dataset_path]
            data = dataset[:]
            
            # Convert to DataFrame
            if len(data.shape) == 1:
                # 1D array
                df = self.pandas.DataFrame(data, columns=['value'])
            elif len(data.shape) == 2:
                # 2D array
                df = self.pandas.DataFrame(data)
            else:
                raise ValueError(f"Cannot convert {len(data.shape)}D array to CSV")
            
            # Write CSV
            df.to_csv(target_path, index=False)
    
    def convert_csv_to_h5(self, source: str, target: str,
                         dataset_path: str = 'data'):
        """
        Convert CSV to HDF5 format.
        
        Args:
            source: Input CSV file path
            target: Output HDF5 file path
            dataset_path: Path for dataset within HDF5 file
        """
        if not self.has_h5py:
            raise ImportError("h5py is required for HDF5 conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for HDF5 conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Read CSV
        df = self.pandas.read_csv(source_path)
        
        # Write HDF5
        with self.h5py.File(target_path, 'w') as f:
            f.create_dataset(dataset_path, data=df.values)
            
            # Store column names as attributes
            f[dataset_path].attrs['columns'] = [str(col) for col in df.columns]
    
    def convert_h5_to_npy(self, source: str, target: str,
                         dataset_path: str = 'data'):
        """
        Convert HDF5 dataset to NumPy array.
        
        Args:
            source: Input HDF5 file path
            target: Output NPY file path
            dataset_path: Path to dataset within HDF5 file
        """
        if not self.has_h5py:
            raise ImportError("h5py is required for HDF5 conversion")
        
        if not self.has_numpy:
            raise ImportError("numpy is required for HDF5 conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        with self.h5py.File(source_path, 'r') as f:
            if dataset_path not in f:
                raise ValueError(f"Dataset '{dataset_path}' not found in HDF5 file")
            
            data = f[dataset_path][:]
            
        # Save as NPY
        self.numpy.save(target_path, data)
    
    # ============================================================
    # NumPy (npy/npz) Conversions
    # ============================================================
    
    def convert_npy_to_csv(self, source: str, target: str):
        """
        Convert NumPy array to CSV format.
        
        Args:
            source: Input NPY file path
            target: Output CSV file path
        """
        if not self.has_numpy:
            raise ImportError("numpy is required for NPY conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for NPY conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Load NPY
        data = self.numpy.load(source_path)
        
        # Convert to DataFrame
        if len(data.shape) == 1:
            df = self.pandas.DataFrame(data, columns=['value'])
        elif len(data.shape) == 2:
            df = self.pandas.DataFrame(data)
        else:
            raise ValueError(f"Cannot convert {len(data.shape)}D array to CSV")
        
        # Write CSV
        df.to_csv(target_path, index=False)
    
    def convert_csv_to_npy(self, source: str, target: str):
        """
        Convert CSV to NumPy array format.
        
        Args:
            source: Input CSV file path
            target: Output NPY file path
        """
        if not self.has_numpy:
            raise ImportError("numpy is required for NPY conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for NPY conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Read CSV
        df = self.pandas.read_csv(source_path)
        
        # Save as NPY
        self.numpy.save(target_path, df.values)
    
    def convert_npy_to_npz(self, source: str, target: str,
                          array_name: str = 'data'):
        """
        Convert NPY to NPZ (compressed) format.
        
        Args:
            source: Input NPY file path
            target: Output NPZ file path
            array_name: Name for array in NPZ archive
        """
        if not self.has_numpy:
            raise ImportError("numpy is required for NPY/NPZ conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Load NPY
        data = self.numpy.load(source_path)
        
        # Save as NPZ
        self.numpy.savez_compressed(target_path, **{array_name: data})
    
    def convert_npz_to_npy(self, source: str, target: str,
                          array_name: Optional[str] = None):
        """
        Convert NPZ to NPY format.
        
        Args:
            source: Input NPZ file path
            target: Output NPY file path
            array_name: Name of array in NPZ archive (first if not specified)
        """
        if not self.has_numpy:
            raise ImportError("numpy is required for NPY/NPZ conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Load NPZ
        npz_data = self.numpy.load(source_path)
        
        # Get array
        if array_name:
            if array_name not in npz_data:
                raise ValueError(f"Array '{array_name}' not found in NPZ file")
            data = npz_data[array_name]
        else:
            # Use first array
            first_key = list(npz_data.keys())[0]
            data = npz_data[first_key]
        
        # Save as NPY
        self.numpy.save(target_path, data)
    
    def convert_npz_to_csv(self, source: str, target: str,
                          array_name: Optional[str] = None):
        """
        Convert NPZ to CSV format.
        
        Args:
            source: Input NPZ file path
            target: Output CSV file path
            array_name: Name of array in NPZ archive (first if not specified)
        """
        if not self.has_numpy:
            raise ImportError("numpy is required for NPZ conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for NPZ conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Load NPZ
        npz_data = self.numpy.load(source_path)
        
        # Get array
        if array_name:
            if array_name not in npz_data:
                raise ValueError(f"Array '{array_name}' not found in NPZ file")
            data = npz_data[array_name]
        else:
            # Use first array
            first_key = list(npz_data.keys())[0]
            data = npz_data[first_key]
        
        # Convert to DataFrame
        if len(data.shape) == 1:
            df = self.pandas.DataFrame(data, columns=['value'])
        elif len(data.shape) == 2:
            df = self.pandas.DataFrame(data)
        else:
            raise ValueError(f"Cannot convert {len(data.shape)}D array to CSV")
        
        # Write CSV
        df.to_csv(target_path, index=False)
    
    # ============================================================
    # Statistics
    # ============================================================
    
    def get_h5_stats(self, file_path: str) -> Dict:
        """
        Get statistics from HDF5 file.
        
        Args:
            file_path: HDF5 file path
            
        Returns:
            Dictionary with statistics
        """
        if not self.has_h5py:
            raise ImportError("h5py is required for HDF5 statistics")
        
        stats = {
            'datasets': [],
            'groups': [],
            'attributes': {}
        }
        
        def visit_item(name, obj):
            if isinstance(obj, self.h5py.Dataset):
                stats['datasets'].append({
                    'name': name,
                    'shape': obj.shape,
                    'dtype': str(obj.dtype),
                    'size': obj.size
                })
            elif isinstance(obj, self.h5py.Group):
                stats['groups'].append(name)
        
        with self.h5py.File(file_path, 'r') as f:
            # Get attributes
            stats['attributes'] = dict(f.attrs)
            
            # Visit all items
            f.visititems(visit_item)
        
        return stats
    
    def get_npy_stats(self, file_path: str) -> Dict:
        """
        Get statistics from NPY file.
        
        Args:
            file_path: NPY file path
            
        Returns:
            Dictionary with statistics
        """
        if not self.has_numpy:
            raise ImportError("numpy is required for NPY statistics")
        
        # Load NPY
        data = self.numpy.load(file_path)
        
        stats = {
            'shape': data.shape,
            'dtype': str(data.dtype),
            'size': data.size,
            'ndim': data.ndim
        }
        
        # Add numeric statistics if applicable
        if self.numpy.issubdtype(data.dtype, self.numpy.number):
            stats['min'] = float(data.min())
            stats['max'] = float(data.max())
            stats['mean'] = float(data.mean())
            stats['std'] = float(data.std())
        
        return stats
    
    def get_npz_stats(self, file_path: str) -> Dict:
        """
        Get statistics from NPZ file.
        
        Args:
            file_path: NPZ file path
            
        Returns:
            Dictionary with statistics
        """
        if not self.has_numpy:
            raise ImportError("numpy is required for NPZ statistics")
        
        # Load NPZ
        npz_data = self.numpy.load(file_path)
        
        stats = {
            'num_arrays': len(npz_data.files),
            'array_names': npz_data.files,
            'arrays': {}
        }
        
        # Get stats for each array
        for name in npz_data.files:
            data = npz_data[name]
            array_stats = {
                'shape': data.shape,
                'dtype': str(data.dtype),
                'size': data.size
            }
            
            if self.numpy.issubdtype(data.dtype, self.numpy.number):
                array_stats['min'] = float(data.min())
                array_stats['max'] = float(data.max())
                array_stats['mean'] = float(data.mean())
            
            stats['arrays'][name] = array_stats
        
        return stats


def get_python_data_converter() -> PythonDataConverter:
    """Get singleton instance of PythonDataConverter."""
    return PythonDataConverter()
