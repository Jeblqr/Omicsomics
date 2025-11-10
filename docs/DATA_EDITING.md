# Data Editing Guide

## Overview

The Data Editing feature provides a comprehensive, interactive interface for performing in-place data transformations and operations on uploaded files. Users can preview changes, validate operations, and apply transformations with automatic backup and undo capabilities.

## Features

### Core Capabilities

- **15+ Operation Types**: Row operations, column operations, value transformations, and type conversions
- **Preview Before Apply**: See exact changes before committing
- **Validation**: Automatic validation of operations and parameters
- **Backup & Undo**: Automatic backup creation with revert capability
- **Session Management**: Save and resume editing sessions
- **Multi-Step Workflow**: Build complex transformation pipelines

### Operation Categories

1. **Row Operations**

   - Delete rows by index
   - Filter rows by conditions
   - Sort rows by columns
   - Remove duplicate rows

2. **Column Operations**

   - Add new columns with expressions
   - Delete columns
   - Rename columns
   - Reorder columns

3. **Value Operations**

   - Replace values (exact or pattern matching)
   - Fill missing values (multiple strategies)
   - Transform values with expressions

4. **Type Operations**
   - Convert data types (int, float, string, datetime, bool)
   - Parse dates with custom formats

## Architecture

### Database Models

#### EditSession

Tracks editing sessions with operations queue:

```python
- id: Primary key
- session_key: Unique identifier
- name: Session name
- description: Session description
- user_id: Owner user ID
- file_id: Source file ID
- status: Session status (draft/previewing/applying/applied/failed/reverted)
- operations: List of operations (JSON)
- preview_data: Preview results (JSON)
- validation_errors: Validation error messages
- output_file_id: Result file ID (after apply)
- backup_path: Backup file path
- rows_affected: Total rows affected
- created_at: Creation timestamp
- updated_at: Update timestamp
```

#### EditOperation

Individual operation records (for history):

```python
- id: Primary key
- session_id: Foreign key to EditSession
- operation_type: Type of operation
- parameters: Operation parameters (JSON)
- undo_data: Data for undo operation (JSON)
- rows_affected: Rows affected by this operation
- executed: Execution flag
- created_at: Creation timestamp
```

### Service Layer

**DataEditingService** provides:

#### Session Management

- `create_edit_session()`: Create new session with backup
- `get_edit_session()`: Get session details
- `list_edit_sessions()`: List sessions with filters
- `delete_edit_session()`: Delete session

#### Operation Management

- `add_operation()`: Add operation to queue
- `remove_operation()`: Remove operation from queue
- `preview_operations()`: Preview all operations (dry-run)
- `apply_operations()`: Apply all operations and create new file
- `revert_operations()`: Restore from backup

#### Data Operations

**Row Operations**:

- `delete_rows`: Delete rows by index range
- `filter_rows`: Filter by column conditions (eq, ne, gt, lt, contains)
- `sort_rows`: Sort by one or more columns
- `deduplicate`: Remove duplicate rows

**Column Operations**:

- `add_column`: Add new column with expression
- `delete_column`: Delete one or more columns
- `rename_column`: Rename column
- `reorder_columns`: Change column order

**Value Operations**:

- `replace_values`: Replace exact values or patterns
- `fill_missing`: Fill NaN values (value, mean, median, forward_fill, backward_fill)
- `transform_values`: Apply expression to values

**Type Operations**:

- `convert_type`: Convert column data type

### API Endpoints

All endpoints are prefixed with `/api/data-editing`.

#### Session Management

**POST /sessions**
Create a new editing session.

Request:

```json
{
  "name": "Clean Sales Data",
  "description": "Remove duplicates and handle missing values",
  "file_id": 123
}
```

Response: `201 Created` with session object

**GET /sessions/{session_id}**
Get session details including operations and preview.

Response:

```json
{
  "id": 1,
  "session_key": "edit_clean_sales_1699123456789",
  "name": "Clean Sales Data",
  "status": "draft",
  "operations": [
    {
      "operation_type": "deduplicate",
      "parameters": {}
    },
    {
      "operation_type": "fill_missing",
      "parameters": {
        "column": "price",
        "method": "mean"
      }
    }
  ],
  "preview_data": {...},
  "validation_errors": [],
  "created_at": "2024-11-10T10:00:00Z"
}
```

**GET /sessions**
List sessions with filters.

Query Parameters:

- `status`: Filter by status
- `skip`: Pagination offset
- `limit`: Page size

Response: Paginated session list

**DELETE /sessions/{session_id}**
Delete a session.

Response: `204 No Content`

#### Operations

**POST /sessions/{session_id}/operations**
Add an operation to the session.

Request:

```json
{
  "operation_type": "filter_rows",
  "parameters": {
    "column": "age",
    "operator": "gt",
    "value": 18
  }
}
```

Response: Updated session object

**DELETE /sessions/{session_id}/operations/{operation_index}**
Remove an operation from the queue.

Response: Updated session object

**POST /sessions/{session_id}/preview**
Preview all operations (dry-run).

Response:

```json
{
  "preview_data": {
    "sample": [...],  // First 10 rows
    "total_rows": 1000,
    "affected_rows": 250,
    "columns": ["col1", "col2"]
  },
  "validation_errors": [],
  "summary": {
    "operations_count": 3,
    "estimated_duration": 2.5
  }
}
```

**POST /sessions/{session_id}/apply**
Apply all operations and create new file.

Response:

```json
{
  "session_id": 1,
  "output_file_id": 456,
  "rows_affected": 250,
  "status": "applied",
  "duration": 2.3
}
```

**POST /sessions/{session_id}/revert**
Revert operations and restore from backup.

Response: Session object with status "reverted"

**GET /operation-types**
Get documentation for all operation types.

Response: Array of operation type documentation with examples

### Frontend

**DataEditingPage** features:

#### Layout

- **Left Sidebar**: Session list with status indicators
- **Main Area**: Stepper interface with 4 steps

#### Stepper Workflow

**Step 1: Select File**

- Create new session dialog
- Enter name, description, select file
- Automatic backup creation

**Step 2: Add Operations**

- Operations list with add/remove buttons
- Operation builder dialog with dynamic parameters
- Visual operation queue

**Step 3: Preview**

- Preview table (first 10 rows)
- Summary statistics
- Validation alerts
- Affected rows count

**Step 4: Apply**

- Success message with output file ID
- Download button for result
- Option to create new session

#### Operation Builder

Dynamic parameter inputs based on operation type:

- **filter_rows**: Column selector, operator dropdown, value input
- **delete_column**: Multi-column selector
- **rename_column**: Old name input, new name input
- **fill_missing**: Column selector, method dropdown, value input (if method=value)
- **convert_type**: Column selector, type dropdown

## Usage Examples

### Example 1: Clean and Filter Sales Data

```python
# Step 1: Create session
POST /api/data-editing/sessions
{
  "name": "Clean Sales Data",
  "description": "Remove duplicates and filter by date range",
  "file_id": 123
}

# Step 2: Add operations
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "deduplicate",
  "parameters": {}
}

POST /api/data-editing/sessions/1/operations
{
  "operation_type": "filter_rows",
  "parameters": {
    "column": "sale_date",
    "operator": "gt",
    "value": "2024-01-01"
  }
}

POST /api/data-editing/sessions/1/operations
{
  "operation_type": "sort_rows",
  "parameters": {
    "columns": ["sale_date"],
    "ascending": [false]
  }
}

# Step 3: Preview
POST /api/data-editing/sessions/1/preview

# Step 4: Apply
POST /api/data-editing/sessions/1/apply
```

### Example 2: Handle Missing Values

```python
# Create session
POST /api/data-editing/sessions
{
  "name": "Handle Missing Values",
  "file_id": 123
}

# Fill numeric column with mean
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "fill_missing",
  "parameters": {
    "column": "price",
    "method": "mean"
  }
}

# Fill categorical column with value
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "fill_missing",
  "parameters": {
    "column": "category",
    "method": "value",
    "fill_value": "Unknown"
  }
}

# Forward fill dates
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "fill_missing",
  "parameters": {
    "column": "date",
    "method": "forward_fill"
  }
}

# Preview and apply
POST /api/data-editing/sessions/1/preview
POST /api/data-editing/sessions/1/apply
```

### Example 3: Type Conversions

```python
# Create session
POST /api/data-editing/sessions
{
  "name": "Convert Data Types",
  "file_id": 123
}

# Convert to integer
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "convert_type",
  "parameters": {
    "column": "quantity",
    "target_type": "int"
  }
}

# Convert to datetime
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "convert_type",
  "parameters": {
    "column": "date",
    "target_type": "datetime"
  }
}

# Convert to float
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "convert_type",
  "parameters": {
    "column": "price",
    "target_type": "float"
  }
}

POST /api/data-editing/sessions/1/apply
```

### Example 4: Column Management

```python
# Create session
POST /api/data-editing/sessions
{
  "name": "Reorganize Columns",
  "file_id": 123
}

# Rename column
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "rename_column",
  "parameters": {
    "old_name": "qty",
    "new_name": "quantity"
  }
}

# Delete unnecessary columns
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "delete_column",
  "parameters": {
    "columns": ["temp_col1", "temp_col2"]
  }
}

# Add calculated column (in future version)
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "add_column",
  "parameters": {
    "name": "total",
    "expression": "quantity * price"
  }
}

POST /api/data-editing/sessions/1/apply
```

### Example 5: Complex Data Cleaning Pipeline

```python
# Create comprehensive cleaning pipeline
POST /api/data-editing/sessions
{
  "name": "Full Data Cleaning Pipeline",
  "description": "Complete data cleaning workflow",
  "file_id": 123
}

# 1. Remove duplicates
POST /api/data-editing/sessions/1/operations
{"operation_type": "deduplicate", "parameters": {}}

# 2. Filter valid rows
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "filter_rows",
  "parameters": {"column": "age", "operator": "gt", "value": 0}
}

# 3. Fill missing values
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "fill_missing",
  "parameters": {"column": "income", "method": "median"}
}

# 4. Replace invalid values
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "replace_values",
  "parameters": {
    "column": "status",
    "old_value": "N/A",
    "new_value": "Unknown"
  }
}

# 5. Convert types
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "convert_type",
  "parameters": {"column": "date", "target_type": "datetime"}
}

# 6. Sort results
POST /api/data-editing/sessions/1/operations
{
  "operation_type": "sort_rows",
  "parameters": {"columns": ["date"], "ascending": [true]}
}

# Preview all operations
POST /api/data-editing/sessions/1/preview

# Apply if preview looks good
POST /api/data-editing/sessions/1/apply
```

## Operation Reference

### filter_rows

Filter rows based on column conditions.

**Parameters**:

- `column` (string): Column name
- `operator` (string): Comparison operator
  - `eq`: Equal to
  - `ne`: Not equal to
  - `gt`: Greater than
  - `lt`: Less than
  - `contains`: Contains substring
- `value`: Comparison value (type depends on column)

**Example**:

```json
{
  "operation_type": "filter_rows",
  "parameters": {
    "column": "age",
    "operator": "gt",
    "value": 18
  }
}
```

### delete_rows

Delete rows by index range.

**Parameters**:

- `start_index` (int): Start index (inclusive)
- `end_index` (int): End index (exclusive)

**Example**:

```json
{
  "operation_type": "delete_rows",
  "parameters": {
    "start_index": 0,
    "end_index": 10
  }
}
```

### sort_rows

Sort rows by one or more columns.

**Parameters**:

- `columns` (list): List of column names
- `ascending` (list): List of boolean values (one per column)

**Example**:

```json
{
  "operation_type": "sort_rows",
  "parameters": {
    "columns": ["date", "amount"],
    "ascending": [false, true]
  }
}
```

### deduplicate

Remove duplicate rows.

**Parameters**:

- `subset` (list, optional): Columns to consider for duplicates
- `keep` (string, optional): Which duplicates to keep ("first", "last", "none")

**Example**:

```json
{
  "operation_type": "deduplicate",
  "parameters": {
    "subset": ["email"],
    "keep": "first"
  }
}
```

### delete_column

Delete one or more columns.

**Parameters**:

- `columns` (list): List of column names to delete

**Example**:

```json
{
  "operation_type": "delete_column",
  "parameters": {
    "columns": ["temp_col", "unused_col"]
  }
}
```

### rename_column

Rename a column.

**Parameters**:

- `old_name` (string): Current column name
- `new_name` (string): New column name

**Example**:

```json
{
  "operation_type": "rename_column",
  "parameters": {
    "old_name": "qty",
    "new_name": "quantity"
  }
}
```

### fill_missing

Fill missing (NaN) values.

**Parameters**:

- `column` (string): Column name
- `method` (string): Fill method
  - `value`: Fill with specific value
  - `mean`: Fill with column mean
  - `median`: Fill with column median
  - `forward_fill`: Fill with previous value
  - `backward_fill`: Fill with next value
- `fill_value` (optional): Value to use (required if method="value")

**Example**:

```json
{
  "operation_type": "fill_missing",
  "parameters": {
    "column": "price",
    "method": "mean"
  }
}
```

### replace_values

Replace values in a column.

**Parameters**:

- `column` (string): Column name
- `old_value`: Value to replace
- `new_value`: Replacement value

**Example**:

```json
{
  "operation_type": "replace_values",
  "parameters": {
    "column": "status",
    "old_value": "N/A",
    "new_value": "Unknown"
  }
}
```

### convert_type

Convert column data type.

**Parameters**:

- `column` (string): Column name
- `target_type` (string): Target type
  - `int`: Integer
  - `float`: Floating point
  - `string`: String/text
  - `datetime`: Datetime
  - `bool`: Boolean

**Example**:

```json
{
  "operation_type": "convert_type",
  "parameters": {
    "column": "date",
    "target_type": "datetime"
  }
}
```

## Best Practices

### Session Management

1. **Use Descriptive Names**: Give sessions clear, descriptive names
2. **Add Descriptions**: Document what transformations you're applying
3. **Preview Before Apply**: Always preview operations before applying
4. **Keep Sessions Small**: Break complex workflows into multiple sessions

### Operation Design

1. **Order Matters**: Operations are applied sequentially
2. **Filter Early**: Remove unnecessary rows early to improve performance
3. **Type Conversions Last**: Convert types after cleaning data
4. **Test Incrementally**: Add operations one at a time and preview

### Performance

1. **Large Files**: Preview only shows sample data (first 10 rows)
2. **Complex Operations**: May take longer on large datasets
3. **Memory Usage**: Be aware of memory limits for very large files

### Data Safety

1. **Automatic Backup**: Session creation automatically creates backup
2. **Revert Available**: Use revert to restore from backup
3. **Original Preserved**: Original file is never modified
4. **New File Created**: Apply creates a new file with results

## Troubleshooting

### Common Issues

**Preview shows no changes**:

- Check that operations are correctly configured
- Verify column names match exactly
- Check filter conditions are valid

**Operation fails during apply**:

- Review validation errors from preview
- Check data types are compatible
- Verify no circular dependencies

**Session stuck in "applying" status**:

- Large files may take time to process
- Check server logs for errors
- Consider breaking into smaller operations

**Cannot revert session**:

- Backup may have been deleted
- Only applied sessions can be reverted
- Create new session if backup is gone

### Debugging Tips

1. **Use Preview**: Always preview before applying
2. **Check Validation**: Review validation errors carefully
3. **Test with Sample**: Test operations on small dataset first
4. **Read Error Messages**: Error messages provide specific guidance

## API Reference Summary

| Endpoint                            | Method | Description        |
| ----------------------------------- | ------ | ------------------ |
| `/sessions`                         | POST   | Create session     |
| `/sessions/{id}`                    | GET    | Get session        |
| `/sessions`                         | GET    | List sessions      |
| `/sessions/{id}`                    | DELETE | Delete session     |
| `/sessions/{id}/operations`         | POST   | Add operation      |
| `/sessions/{id}/operations/{index}` | DELETE | Remove operation   |
| `/sessions/{id}/preview`            | POST   | Preview operations |
| `/sessions/{id}/apply`              | POST   | Apply operations   |
| `/sessions/{id}/revert`             | POST   | Revert operations  |
| `/operation-types`                  | GET    | Get operation docs |

## Future Enhancements

- **Custom Expressions**: User-defined transformation expressions
- **Conditional Operations**: Apply operations based on conditions
- **Batch Sessions**: Apply same operations to multiple files
- **Operation Templates**: Save and reuse operation sequences
- **Collaborative Editing**: Share sessions with team members
- **Change Tracking**: Detailed audit log of all changes
- **Undo/Redo**: Step-by-step undo and redo
- **Data Validation Rules**: Define and enforce data quality rules
