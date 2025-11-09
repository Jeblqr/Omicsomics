# API Response Format Standards

## 🎯 Unified API Response Formats

This document defines the **standardized response formats** for all API endpoints to ensure consistency across the platform.

---

## 📋 HTTP Status Codes

### Success Codes

- **200 OK** - Successful GET, PUT requests (returns data)
- **201 CREATED** - Successful POST requests (returns created resource)
- **204 NO CONTENT** - Successful DELETE requests (no response body)

### Error Codes

- **400 BAD REQUEST** - Invalid request data
- **401 UNAUTHORIZED** - Not authenticated
- **403 FORBIDDEN** - Not authorized to access resource
- **404 NOT FOUND** - Resource not found
- **422 UNPROCESSABLE ENTITY** - Validation error
- **500 INTERNAL SERVER ERROR** - Server error

---

## 🔧 CRUD Operations Standards

### CREATE (POST)

**Status Code:** `201 CREATED`

**Response Format:**

```json
{
  "id": 123,
  "name": "Resource Name",
  "description": "Description",
  "created_at": "2025-01-09T10:00:00Z",
  "updated_at": "2025-01-09T10:00:00Z",
  ...other_fields
}
```

**Examples:**

- `POST /projects/` → 201 + Project object
- `POST /runs/` → 201 + Run object
- `POST /data/upload` → 201 + DataFile object

---

### READ (GET)

#### Single Resource

**Status Code:** `200 OK`

**Response Format:**

```json
{
  "id": 123,
  "name": "Resource Name",
  ...fields
}
```

#### List Resources

**Status Code:** `200 OK`

**Response Format:**

```json
[
  {
    "id": 1,
    "name": "Resource 1",
    ...fields
  },
  {
    "id": 2,
    "name": "Resource 2",
    ...fields
  }
]
```

**Examples:**

- `GET /projects/` → 200 + Array of Projects
- `GET /projects/{id}` → 200 + Single Project
- `GET /runs/` → 200 + Array of Runs
- `GET /data/` → 200 + Array of DataFiles

---

### UPDATE (PUT)

**Status Code:** `200 OK`

**Response Format:**

```json
{
  "id": 123,
  "name": "Updated Name",
  "description": "Updated Description",
  "updated_at": "2025-01-09T11:00:00Z",
  ...other_fields
}
```

**Examples:**

- `PUT /projects/{id}` → 200 + Updated Project
- `PUT /runs/{id}` → 200 + Updated Run

---

### DELETE

**Status Code:** `204 NO CONTENT`

**Response:** No body (empty response)

**Examples:**

- `DELETE /projects/{id}` → 204 (no body)
- `DELETE /runs/{id}` → 204 (no body)
- `DELETE /data/{id}` → 204 (no body)
- `DELETE /pipelines/{id}` → 204 (no body)

**Frontend Handling:**

```typescript
// Correct way to handle 204 responses
await api.delete(`/projects/${id}`);
// No need to parse response - 204 means success
await fetchProjects(); // Refresh list
```

---

## 📊 Common Response Fields

### All Resources Should Include:

```typescript
{
  id: number; // Unique identifier
  created_at: string; // ISO 8601 timestamp
  updated_at: string; // ISO 8601 timestamp
}
```

### Projects

```typescript
{
  id: number;
  name: string;
  description: string;
  owner_id: number;
  created_at: string;
  updated_at: string;
}
```

### Runs

```typescript
{
  id: number;
  name: string;
  description: string;
  status: "pending" | "running" | "completed" | "failed" | "cancelled";
  pipeline_type?: string;
  pipeline_template_id?: string;
  custom_pipeline_id?: number;
  parameters?: Record<string, any>;
  input_files?: number[];           // Array of DataFile IDs
  output_files?: number[];
  progress?: number;                // 0-100
  logs?: string;
  error_message?: string;
  project_id: number;
  owner_id: number;
  started_at?: string;
  finished_at?: string;
  created_at: string;
  updated_at: string;
}
```

### DataFiles

```typescript
{
  id: number;
  filename: string;
  object_key: string;
  metadata_: Record<string, unknown>;
  size: number;
  checksum: string;
  project_id: number;
  run_id?: number;
  uploaded_by_id: number;
  created_at: string;
  updated_at: string;
}
```

### Custom Pipelines

```typescript
{
  id: number;
  name: string;
  description?: string;
  nodes: Array<{
    id: string;
    type: string;
    position: { x: number; y: number };
    data: Record<string, any>;
  }>;
  edges: Array<{
    id: string;
    source: string;
    target: string;
  }>;
  project_id: number;
  owner_id: number;
  created_at: string;
  updated_at: string;
}
```

---

## ⚠️ Error Response Format

### Standard Error Response

**Format:**

```json
{
  "detail": "Error message describing what went wrong"
}
```

**Examples:**

```json
// 404 Not Found
{
  "detail": "Project not found"
}

// 403 Forbidden
{
  "detail": "Not authorized"
}

// 422 Validation Error
{
  "detail": [
    {
      "loc": ["body", "name"],
      "msg": "field required",
      "type": "value_error.missing"
    }
  ]
}
```

---

## ✅ Implementation Checklist

### Backend (FastAPI)

- [x] All DELETE endpoints return `204 NO CONTENT`
- [x] All POST (create) endpoints return `201 CREATED`
- [x] All GET endpoints return `200 OK`
- [x] All PUT (update) endpoints return `200 OK`
- [x] Error responses include `detail` field
- [x] Consistent timestamp format (ISO 8601)
- [x] Response models match schemas

### Frontend (React/TypeScript)

- [x] Use TypeScript interfaces matching backend schemas
- [x] Handle 204 responses correctly (no body parsing)
- [x] Expect 201 for creates, 200 for reads/updates
- [x] Parse error responses from `error.response?.data?.detail`
- [x] Refresh lists after mutations (create/update/delete)

---

## 🔍 Issues Fixed

### Before (Inconsistent)

```typescript
// Projects delete returned 200 with body
DELETE /projects/123 → 200 {"detail": "Project deleted"}

// Runs delete returned 204 with no body
DELETE /runs/456 → 204 (no body)

// Create returned 200 instead of 201
POST /projects/ → 200 {...}
```

### After (Consistent)

```typescript
// All deletes return 204
DELETE /projects/123 → 204 (no body)
DELETE /runs/456 → 204 (no body)
DELETE /data/789 → 204 (no body)

// All creates return 201
POST /projects/ → 201 {...}
POST /runs/ → 201 {...}
POST /data/upload → 201 {...}
```

---

## 🧪 Testing

### Manual Test

```bash
# Test create (should return 201)
curl -X POST http://localhost:8000/projects/ \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"name":"Test","description":"Test"}' \
  -w "\nStatus: %{http_code}\n"

# Test delete (should return 204)
curl -X DELETE http://localhost:8000/projects/123 \
  -H "Authorization: Bearer $TOKEN" \
  -w "\nStatus: %{http_code}\n"
```

### Automated Test

```bash
python3 scripts/test_api_formats.py
```

---

## 📝 Notes

### Why 204 for DELETE?

- **RESTful Standard:** DELETE operations that succeed but don't return data should use 204
- **Bandwidth:** Saves bandwidth by not sending unnecessary response body
- **Clarity:** Client knows operation succeeded if status is 204, no need to parse body

### Why 201 for CREATE?

- **RESTful Standard:** 201 indicates a new resource was created
- **Semantics:** Distinguishes create from update (200)
- **Location Header:** Can include `Location` header with new resource URL

### Frontend Handling

```typescript
// DELETE - no response body
await api.delete(`/projects/${id}`);
// Success if no error thrown

// CREATE - parse response body
const response = await api.post("/projects/", data);
const newProject = response.data; // Available for 201

// UPDATE - parse response body
const response = await api.put(`/projects/${id}`, data);
const updated = response.data; // Available for 200
```

---

## 🎯 Summary

**Standardized Status Codes:**

- `201` - Create success (with body)
- `200` - Read/Update success (with body)
- `204` - Delete success (no body)
- `404` - Not found (with error detail)
- `403` - Forbidden (with error detail)

**All endpoints now follow REST conventions consistently!**
