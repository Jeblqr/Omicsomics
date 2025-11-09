"""
Data Quality Control Service.

Generates comprehensive quality reports for omics data files:
- Missing value statistics
- Outlier detection
- Distribution analysis
- Data completeness metrics
"""

import numpy as np
from typing import Dict, List, Any
from app.schemas.unified_format import UnifiedData


class QualityControlService:
    """Service for generating data quality reports."""
    
    @staticmethod
    def generate_report(unified_data: UnifiedData) -> Dict[str, Any]:
        """Generate comprehensive quality control report."""
        records = unified_data.records
        headers = unified_data.headers
        
        if not records:
            return {
                "error": "No data records available",
                "record_count": 0,
            }
        
        # Calculate metrics
        missing_stats = QualityControlService._calculate_missing_values(records, headers)
        outlier_stats = QualityControlService._detect_outliers(records, headers)
        distribution_stats = QualityControlService._analyze_distribution(records, headers)
        completeness_stats = QualityControlService._calculate_completeness(records, headers)
        
        return {
            "record_count": len(records),
            "field_count": len(headers),
            "missing_values": missing_stats,
            "outliers": outlier_stats,
            "distributions": distribution_stats,
            "completeness": completeness_stats,
            "overall_quality_score": QualityControlService._calculate_quality_score(
                missing_stats, outlier_stats, completeness_stats
            ),
        }
    
    @staticmethod
    def _calculate_missing_values(records: List[Any], headers: List[str]) -> Dict[str, Any]:
        """Calculate missing value statistics for each field."""
        total_records = len(records)
        field_stats = {}
        
        for header in headers:
            missing_count = 0
            for record in records:
                # Check both record level and values dict
                value = None
                if hasattr(record, header):
                    value = getattr(record, header)
                elif hasattr(record, 'values') and header in record.values:
                    value = record.values[header]
                
                if value is None or (isinstance(value, str) and value.strip() == ""):
                    missing_count += 1
            
            missing_ratio = missing_count / total_records if total_records > 0 else 0
            field_stats[header] = {
                "missing_count": missing_count,
                "missing_ratio": missing_ratio,
                "status": QualityControlService._get_missing_status(missing_ratio),
            }
        
        # Overall statistics
        total_missing = sum(s["missing_count"] for s in field_stats.values())
        total_values = total_records * len(headers)
        overall_missing_ratio = total_missing / total_values if total_values > 0 else 0
        
        return {
            "by_field": field_stats,
            "overall_missing_count": total_missing,
            "overall_missing_ratio": overall_missing_ratio,
            "fields_with_high_missing": [
                field for field, stats in field_stats.items()
                if stats["missing_ratio"] > 0.5
            ],
        }
    
    @staticmethod
    def _detect_outliers(records: List[Any], headers: List[str]) -> Dict[str, Any]:
        """Detect outliers using IQR method for numeric fields."""
        outlier_stats = {}
        
        for header in headers:
            values = []
            for record in records:
                value = None
                if hasattr(record, header):
                    value = getattr(record, header)
                elif hasattr(record, 'values') and header in record.values:
                    value = record.values[header]
                
                # Try to convert to float
                try:
                    if value is not None and value != "":
                        values.append(float(value))
                except (ValueError, TypeError):
                    continue
            
            if len(values) < 4:  # Need at least 4 values for IQR
                continue
            
            # Calculate IQR
            q1 = np.percentile(values, 25)
            q3 = np.percentile(values, 75)
            iqr = q3 - q1
            lower_bound = q1 - 1.5 * iqr
            upper_bound = q3 + 1.5 * iqr
            
            outliers = [v for v in values if v < lower_bound or v > upper_bound]
            outlier_ratio = len(outliers) / len(values) if len(values) > 0 else 0
            
            outlier_stats[header] = {
                "outlier_count": len(outliers),
                "outlier_ratio": outlier_ratio,
                "lower_bound": lower_bound,
                "upper_bound": upper_bound,
                "status": QualityControlService._get_outlier_status(outlier_ratio),
            }
        
        return {
            "by_field": outlier_stats,
            "fields_with_high_outliers": [
                field for field, stats in outlier_stats.items()
                if stats["outlier_ratio"] > 0.1
            ],
        }
    
    @staticmethod
    def _analyze_distribution(records: List[Any], headers: List[str]) -> Dict[str, Any]:
        """Analyze data distribution for numeric fields."""
        distribution_stats = {}
        
        for header in headers:
            values = []
            for record in records:
                value = None
                if hasattr(record, header):
                    value = getattr(record, header)
                elif hasattr(record, 'values') and header in record.values:
                    value = record.values[header]
                
                try:
                    if value is not None and value != "":
                        values.append(float(value))
                except (ValueError, TypeError):
                    continue
            
            if len(values) < 2:
                continue
            
            distribution_stats[header] = {
                "min": float(np.min(values)),
                "max": float(np.max(values)),
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "std": float(np.std(values)),
                "q25": float(np.percentile(values, 25)),
                "q75": float(np.percentile(values, 75)),
                "value_count": len(values),
            }
        
        return distribution_stats
    
    @staticmethod
    def _calculate_completeness(records: List[Any], headers: List[str]) -> Dict[str, Any]:
        """Calculate data completeness metrics."""
        total_records = len(records)
        complete_records = 0
        
        for record in records:
            is_complete = True
            for header in headers:
                value = None
                if hasattr(record, header):
                    value = getattr(record, header)
                elif hasattr(record, 'values') and header in record.values:
                    value = record.values[header]
                
                if value is None or (isinstance(value, str) and value.strip() == ""):
                    is_complete = False
                    break
            
            if is_complete:
                complete_records += 1
        
        completeness_ratio = complete_records / total_records if total_records > 0 else 0
        
        return {
            "complete_records": complete_records,
            "incomplete_records": total_records - complete_records,
            "completeness_ratio": completeness_ratio,
            "status": QualityControlService._get_completeness_status(completeness_ratio),
        }
    
    @staticmethod
    def _calculate_quality_score(
        missing_stats: Dict[str, Any],
        outlier_stats: Dict[str, Any],
        completeness_stats: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Calculate overall quality score (0-100)."""
        # Scoring weights
        completeness_weight = 0.4
        missing_weight = 0.4
        outlier_weight = 0.2
        
        # Completeness score
        completeness_score = completeness_stats["completeness_ratio"] * 100
        
        # Missing values score (inverse)
        missing_score = (1 - missing_stats["overall_missing_ratio"]) * 100
        
        # Outlier score (inverse average of outlier ratios)
        outlier_ratios = [
            stats["outlier_ratio"]
            for stats in outlier_stats.get("by_field", {}).values()
        ]
        avg_outlier_ratio = np.mean(outlier_ratios) if outlier_ratios else 0
        outlier_score = (1 - avg_outlier_ratio) * 100
        
        # Weighted overall score
        overall_score = (
            completeness_score * completeness_weight +
            missing_score * missing_weight +
            outlier_score * outlier_weight
        )
        
        return {
            "score": round(overall_score, 2),
            "completeness_score": round(completeness_score, 2),
            "missing_values_score": round(missing_score, 2),
            "outlier_score": round(outlier_score, 2),
            "rating": QualityControlService._get_quality_rating(overall_score),
        }
    
    @staticmethod
    def _get_missing_status(ratio: float) -> str:
        """Get status label for missing value ratio."""
        if ratio == 0:
            return "excellent"
        elif ratio < 0.05:
            return "good"
        elif ratio < 0.2:
            return "acceptable"
        elif ratio < 0.5:
            return "poor"
        else:
            return "critical"
    
    @staticmethod
    def _get_outlier_status(ratio: float) -> str:
        """Get status label for outlier ratio."""
        if ratio < 0.01:
            return "excellent"
        elif ratio < 0.05:
            return "good"
        elif ratio < 0.1:
            return "acceptable"
        else:
            return "high"
    
    @staticmethod
    def _get_completeness_status(ratio: float) -> str:
        """Get status label for completeness ratio."""
        if ratio >= 0.95:
            return "excellent"
        elif ratio >= 0.8:
            return "good"
        elif ratio >= 0.5:
            return "acceptable"
        else:
            return "poor"
    
    @staticmethod
    def _get_quality_rating(score: float) -> str:
        """Get overall quality rating."""
        if score >= 90:
            return "A"
        elif score >= 80:
            return "B"
        elif score >= 70:
            return "C"
        elif score >= 60:
            return "D"
        else:
            return "F"
