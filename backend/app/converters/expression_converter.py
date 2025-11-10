"""
Expression matrix format converter.

Supports conversions between:
- MTX (Market Matrix Format - sparse matrices)
- 10X HDF5 (10X Genomics format)
- h5ad (AnnData format)
- CSV (dense matrix)

Dependencies:
- scipy: Sparse matrix operations
- scanpy: Single-cell analysis (for h5ad)
- h5py: HDF5 file handling
"""
from pathlib import Path
from typing import Dict, Optional, List
import csv


class ExpressionConverter:
    """
    Converter for expression matrix formats.
    
    Singleton pattern for efficiency.
    """
    
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def __init__(self):
        """Initialize expression converter."""
        if not hasattr(self, '_initialized'):
            self._initialized = True
            self._check_dependencies()
    
    def _check_dependencies(self):
        """Check if required tools are available."""
        try:
            import scipy.io
            import scipy.sparse
            self.scipy = scipy
            self.has_scipy = True
        except ImportError:
            self.has_scipy = False
        
        try:
            import scanpy as sc
            self.scanpy = sc
            self.has_scanpy = True
        except ImportError:
            self.has_scanpy = False
        
        try:
            import h5py
            self.h5py = h5py
            self.has_h5py = True
        except ImportError:
            self.has_h5py = False
        
        try:
            import pandas as pd
            self.pandas = pd
            self.has_pandas = True
        except ImportError:
            self.has_pandas = False
    
    # ============================================================
    # MTX (Market Matrix) Conversions
    # ============================================================
    
    def convert_mtx_to_csv(self, source: str, target: str,
                          genes_file: Optional[str] = None,
                          barcodes_file: Optional[str] = None):
        """
        Convert MTX (sparse matrix) to CSV (dense matrix).
        
        Args:
            source: Input MTX file path
            target: Output CSV file path
            genes_file: Optional genes/features file (one per line)
            barcodes_file: Optional barcodes/cells file (one per line)
        """
        if not self.has_scipy:
            raise ImportError("scipy is required for MTX conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for MTX conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Read MTX file (sparse matrix)
        matrix = self.scipy.io.mmread(source_path)
        
        # Convert to dense if needed
        if self.scipy.sparse.issparse(matrix):
            matrix = matrix.toarray()
        
        # Read gene names if provided
        if genes_file:
            with open(genes_file, 'r') as f:
                genes = [line.strip().split('\t')[0] for line in f]
        else:
            genes = [f"Gene_{i}" for i in range(matrix.shape[0])]
        
        # Read barcode names if provided
        if barcodes_file:
            with open(barcodes_file, 'r') as f:
                barcodes = [line.strip() for line in f]
        else:
            barcodes = [f"Cell_{i}" for i in range(matrix.shape[1])]
        
        # Create DataFrame
        df = self.pandas.DataFrame(matrix, index=genes, columns=barcodes)
        
        # Write CSV
        df.to_csv(target_path)
    
    def convert_csv_to_mtx(self, source: str, target: str,
                          genes_file: Optional[str] = None,
                          barcodes_file: Optional[str] = None):
        """
        Convert CSV (dense matrix) to MTX (sparse matrix).
        
        Args:
            source: Input CSV file path
            target: Output MTX file path
            genes_file: Optional output genes file
            barcodes_file: Optional output barcodes file
        """
        if not self.has_scipy:
            raise ImportError("scipy is required for MTX conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for MTX conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Read CSV
        df = self.pandas.read_csv(source_path, index_col=0)
        
        # Convert to sparse matrix
        sparse_matrix = self.scipy.sparse.csr_matrix(df.values)
        
        # Write MTX file
        self.scipy.io.mmwrite(target_path, sparse_matrix)
        
        # Write genes file if requested
        if genes_file:
            with open(genes_file, 'w') as f:
                for gene in df.index:
                    f.write(f"{gene}\n")
        
        # Write barcodes file if requested
        if barcodes_file:
            with open(barcodes_file, 'w') as f:
                for barcode in df.columns:
                    f.write(f"{barcode}\n")
    
    # ============================================================
    # 10X HDF5 Conversions
    # ============================================================
    
    def convert_10x_h5_to_h5ad(self, source: str, target: str):
        """
        Convert 10X HDF5 to h5ad (AnnData) format.
        
        Args:
            source: Input 10X HDF5 file path
            target: Output h5ad file path
        """
        if not self.has_scanpy:
            raise ImportError("scanpy is required for 10X HDF5 conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Read 10X HDF5
        adata = self.scanpy.read_10x_h5(source_path)
        
        # Write h5ad
        adata.write_h5ad(target_path)
    
    def convert_10x_h5_to_csv(self, source: str, target: str,
                              layer: Optional[str] = None):
        """
        Convert 10X HDF5 to CSV format.
        
        Args:
            source: Input 10X HDF5 file path
            target: Output CSV file path
            layer: Optional layer to export (default: X)
        """
        if not self.has_scanpy:
            raise ImportError("scanpy is required for 10X HDF5 conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for 10X HDF5 conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Read 10X HDF5
        adata = self.scanpy.read_10x_h5(source_path)
        
        # Get matrix
        if layer and layer in adata.layers:
            matrix = adata.layers[layer]
        else:
            matrix = adata.X
        
        # Convert to dense if sparse
        if self.scipy.sparse.issparse(matrix):
            matrix = matrix.toarray()
        
        # Create DataFrame
        df = self.pandas.DataFrame(
            matrix,
            index=adata.obs_names,
            columns=adata.var_names
        )
        
        # Write CSV
        df.to_csv(target_path)
    
    # ============================================================
    # h5ad (AnnData) Conversions
    # ============================================================
    
    def convert_h5ad_to_csv(self, source: str, target: str,
                           layer: Optional[str] = None):
        """
        Convert h5ad (AnnData) to CSV format.
        
        Args:
            source: Input h5ad file path
            target: Output CSV file path
            layer: Optional layer to export (default: X)
        """
        if not self.has_scanpy:
            raise ImportError("scanpy is required for h5ad conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for h5ad conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Read h5ad
        adata = self.scanpy.read_h5ad(source_path)
        
        # Get matrix
        if layer and layer in adata.layers:
            matrix = adata.layers[layer]
        else:
            matrix = adata.X
        
        # Convert to dense if sparse
        if self.scipy.sparse.issparse(matrix):
            matrix = matrix.toarray()
        
        # Create DataFrame
        df = self.pandas.DataFrame(
            matrix,
            index=adata.obs_names,
            columns=adata.var_names
        )
        
        # Write CSV
        df.to_csv(target_path)
    
    def convert_csv_to_h5ad(self, source: str, target: str):
        """
        Convert CSV to h5ad (AnnData) format.
        
        Args:
            source: Input CSV file path
            target: Output h5ad file path
        """
        if not self.has_scanpy:
            raise ImportError("scanpy is required for h5ad conversion")
        
        if not self.has_pandas:
            raise ImportError("pandas is required for h5ad conversion")
        
        source_path = Path(source)
        target_path = Path(target)
        
        # Read CSV
        df = self.pandas.read_csv(source_path, index_col=0)
        
        # Create AnnData object
        adata = self.scanpy.AnnData(df.T)  # Transpose: genes x cells -> cells x genes
        
        # Write h5ad
        adata.write_h5ad(target_path)
    
    # ============================================================
    # Statistics
    # ============================================================
    
    def get_mtx_stats(self, file_path: str,
                     genes_file: Optional[str] = None,
                     barcodes_file: Optional[str] = None) -> Dict:
        """
        Get statistics from MTX file.
        
        Args:
            file_path: MTX file path
            genes_file: Optional genes file
            barcodes_file: Optional barcodes file
            
        Returns:
            Dictionary with statistics
        """
        if not self.has_scipy:
            raise ImportError("scipy is required for MTX statistics")
        
        # Read MTX file
        matrix = self.scipy.io.mmread(file_path)
        
        stats = {
            'shape': matrix.shape,
            'num_genes': matrix.shape[0],
            'num_cells': matrix.shape[1],
            'total_counts': matrix.sum(),
            'is_sparse': self.scipy.sparse.issparse(matrix)
        }
        
        if self.scipy.sparse.issparse(matrix):
            stats['sparsity'] = 1.0 - (matrix.nnz / (matrix.shape[0] * matrix.shape[1]))
            stats['nnz'] = matrix.nnz
        
        # Count genes and barcodes from files if provided
        if genes_file:
            with open(genes_file, 'r') as f:
                stats['num_genes_file'] = sum(1 for _ in f)
        
        if barcodes_file:
            with open(barcodes_file, 'r') as f:
                stats['num_barcodes_file'] = sum(1 for _ in f)
        
        return stats
    
    def get_h5ad_stats(self, file_path: str) -> Dict:
        """
        Get statistics from h5ad file.
        
        Args:
            file_path: h5ad file path
            
        Returns:
            Dictionary with statistics
        """
        if not self.has_scanpy:
            raise ImportError("scanpy is required for h5ad statistics")
        
        # Read h5ad
        adata = self.scanpy.read_h5ad(file_path)
        
        stats = {
            'n_obs': adata.n_obs,
            'n_vars': adata.n_vars,
            'obs_keys': list(adata.obs.columns),
            'var_keys': list(adata.var.columns),
            'layers': list(adata.layers.keys()) if adata.layers else [],
            'obsm_keys': list(adata.obsm.keys()) if adata.obsm else [],
            'varm_keys': list(adata.varm.keys()) if adata.varm else [],
            'uns_keys': list(adata.uns.keys()) if adata.uns else []
        }
        
        # Check if sparse
        if self.scipy.sparse.issparse(adata.X):
            stats['is_sparse'] = True
            stats['sparsity'] = 1.0 - (adata.X.nnz / (adata.n_obs * adata.n_vars))
        else:
            stats['is_sparse'] = False
        
        return stats
    
    def get_10x_h5_stats(self, file_path: str) -> Dict:
        """
        Get statistics from 10X HDF5 file.
        
        Args:
            file_path: 10X HDF5 file path
            
        Returns:
            Dictionary with statistics
        """
        if not self.has_h5py:
            raise ImportError("h5py is required for 10X HDF5 statistics")
        
        stats = {}
        
        with self.h5py.File(file_path, 'r') as f:
            # Try to get matrix group
            if 'matrix' in f:
                matrix_group = f['matrix']
                stats['shape'] = matrix_group['shape'][:]
                stats['nnz'] = len(matrix_group['data'])
                
                if 'features' in matrix_group:
                    stats['num_features'] = len(matrix_group['features']['name'])
                if 'barcodes' in matrix_group:
                    stats['num_barcodes'] = len(matrix_group['barcodes'])
        
        return stats


def get_expression_converter() -> ExpressionConverter:
    """Get singleton instance of ExpressionConverter."""
    return ExpressionConverter()
