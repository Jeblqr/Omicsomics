"""
Advanced search service for unified data.

Provides full-text search and range queries on processed omics data.
"""

import logging
from typing import Dict, List, Optional, Any
import json

from app.schemas.unified_format import UnifiedData

logger = logging.getLogger(__name__)


class SearchQuery:
    """Search query builder."""

    def __init__(self):
        self.text_query: Optional[str] = None
        self.feature_id_query: Optional[str] = None
        self.range_queries: Dict[str, Dict[str, float]] = {}
        self.exact_matches: Dict[str, Any] = {}

    def with_text(self, query: str) -> "SearchQuery":
        """Add full-text search query."""
        self.text_query = query.lower()
        return self

    def with_feature_id(self, feature_id: str) -> "SearchQuery":
        """Search by feature ID (gene, protein, metabolite)."""
        self.feature_id_query = feature_id.lower()
        return self

    def with_range(
        self,
        field: str,
        min_val: Optional[float] = None,
        max_val: Optional[float] = None,
    ) -> "SearchQuery":
        """Add range query for numeric field."""
        if min_val is not None or max_val is not None:
            self.range_queries[field] = {}
            if min_val is not None:
                self.range_queries[field]["min"] = min_val
            if max_val is not None:
                self.range_queries[field]["max"] = max_val
        return self

    def with_exact_match(self, field: str, value: Any) -> "SearchQuery":
        """Add exact match filter."""
        self.exact_matches[field] = value
        return self


class DataSearchService:
    """Service for searching unified data."""

    @staticmethod
    def search(unified_data: UnifiedData, query: SearchQuery) -> List[Dict[str, Any]]:
        """
        Search unified data with query.

        Args:
            unified_data: Unified data to search
            query: Search query

        Returns:
            List of matching records
        """
        results = []

        for record in unified_data.records:
            if DataSearchService._matches_query(record, query):
                results.append(record.model_dump())

        return results

    @staticmethod
    def _matches_query(record: Any, query: SearchQuery) -> bool:
        """Check if record matches query."""
        record_dict = record.model_dump()

        # Get the values dictionary
        values = record_dict.get("values", {})
        search_dict = {**record_dict, **values}

        # Text search
        if query.text_query:
            if not DataSearchService._matches_text(search_dict, query.text_query):
                return False

        # Feature ID search
        if query.feature_id_query:
            record_id = record_dict.get("id", "")
            if not record_id or query.feature_id_query not in record_id.lower():
                return False

        # Range queries
        for field, range_spec in query.range_queries.items():
            if not DataSearchService._matches_range(search_dict, field, range_spec):
                return False

        # Exact matches
        for field, value in query.exact_matches.items():
            if search_dict.get(field) != value:
                return False

        return True

    @staticmethod
    def _matches_text(feature_dict: Dict[str, Any], query: str) -> bool:
        """Check if feature contains text query."""
        # Search in string fields
        searchable_fields = [
            "feature_id",
            "gene_name",
            "gene_id",
            "protein_id",
            "metabolite_id",
            "description",
            "annotation",
        ]

        for field in searchable_fields:
            value = feature_dict.get(field)
            if value and isinstance(value, str) and query in value.lower():
                return True

        # Search in nested structures
        for key, value in feature_dict.items():
            if isinstance(value, str) and query in value.lower():
                return True
            elif isinstance(value, dict):
                value_str = json.dumps(value).lower()
                if query in value_str:
                    return True

        return False

    @staticmethod
    def _matches_range(
        feature_dict: Dict[str, Any], field: str, range_spec: Dict[str, float]
    ) -> bool:
        """Check if feature value is within range."""
        value = feature_dict.get(field)

        if value is None:
            return False

        try:
            numeric_value = float(value)
        except (ValueError, TypeError):
            return False

        if "min" in range_spec and numeric_value < range_spec["min"]:
            return False

        if "max" in range_spec and numeric_value > range_spec["max"]:
            return False

        return True

    @staticmethod
    def get_searchable_fields(unified_data: UnifiedData) -> List[str]:
        """Get list of searchable fields from data."""
        if not unified_data.records:
            return []

        # Get headers from metadata
        fields = list(unified_data.headers)

        # Add 'id' field
        if "id" not in fields:
            fields.insert(0, "id")

        return fields

    @staticmethod
    def get_field_stats(unified_data: UnifiedData, field: str) -> Dict[str, Any]:
        """Get statistics for a numeric field."""
        values = []

        for record in unified_data.records:
            record_dict = record.model_dump()
            # Check both record level and values dict
            value = record_dict.get(field) or record_dict.get("values", {}).get(field)

            if value is not None:
                try:
                    values.append(float(value))
                except (ValueError, TypeError):
                    continue

        if not values:
            return {"field": field, "count": 0, "type": "non-numeric"}

        values.sort()
        n = len(values)

        return {
            "field": field,
            "count": n,
            "type": "numeric",
            "min": min(values),
            "max": max(values),
            "mean": sum(values) / n,
            "median": (
                values[n // 2]
                if n % 2 == 1
                else (values[n // 2 - 1] + values[n // 2]) / 2
            ),
        }
