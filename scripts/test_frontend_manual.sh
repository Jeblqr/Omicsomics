#!/bin/bash

# Manual Frontend Testing Script
# This script provides manual test instructions since Playwright requires additional setup

echo "==================================="
echo "Frontend Manual Testing Guide"
echo "==================================="
echo ""

# Check if services are running
echo "1. Checking if services are running..."
cd /home/jeblqr/data1/projects/Omicsomics/infrastructure
docker compose ps

echo ""
echo "2. Testing Backend Health..."
curl -s http://localhost:8001/healthz | jq '.'

echo ""
echo "==================================="
echo "Manual Testing Steps:"
echo "==================================="
echo ""
echo "📋 Step 1: Access Frontend"
echo "   URL: http://localhost:5173"
echo ""
echo "📋 Step 2: Login with Demo Account"
echo "   Email: demo@omicsomics.com"
echo "   Password: demo123456"
echo ""
echo "📋 Step 3: Test Each Module"
echo "   ✅ Dashboard - Check if text is visible"
echo "   ✅ Projects - Create/view projects"
echo "   ✅ Data - Upload files (warning message should be DARK, not blurred)"
echo "   ✅ Pipelines - Should show templates (not blank page)"
echo "   ✅ Custom Pipelines - Check functionality"
echo "   ✅ Runs - Check text visibility"
echo ""
echo "📋 Step 4: Verify Text Colors"
echo "   All text should have proper contrast:"
echo "   - Titles: Dark (#212529)"
echo "   - Descriptions: Medium gray (#6c757d)"
echo "   - Warnings: Golden brown (#856404)"
echo ""
echo "📋 Step 5: Test File Upload"
echo "   1. Select a project"
echo "   2. Go to Data page"
echo "   3. Click '+ Upload File'"
echo "   4. Select a test file from test_data/"
echo "   5. Verify upload succeeds"
echo ""
echo "==================================="
echo "Quick API Tests:"
echo "==================================="

# Get auth token
echo ""
echo "Getting auth token..."
TOKEN=$(curl -s -X POST http://localhost:8001/api/v1/auth/login/access-token \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=demo@omicsomics.com&password=demo123456" | jq -r '.access_token')

if [ -n "$TOKEN" ] && [ "$TOKEN" != "null" ]; then
  echo "✅ Login successful"
  
  echo ""
  echo "Testing /api/v1/pipelines/ endpoint..."
  PIPELINES=$(curl -s -H "Authorization: Bearer $TOKEN" http://localhost:8001/api/v1/pipelines/)
  PIPELINE_COUNT=$(echo "$PIPELINES" | jq '. | length')
  echo "✅ Found $PIPELINE_COUNT pipeline templates"
  
  echo ""
  echo "Testing /api/v1/projects/ endpoint..."
  PROJECTS=$(curl -s -H "Authorization: Bearer $TOKEN" http://localhost:8001/api/v1/projects/)
  PROJECT_COUNT=$(echo "$PROJECTS" | jq '. | length')
  echo "✅ Found $PROJECT_COUNT projects"
else
  echo "❌ Login failed"
fi

echo ""
echo "==================================="
echo "Summary of Recent Fixes:"
echo "==================================="
echo "✅ Data page: Fixed blurred warning text (#856404 color)"
echo "✅ Data page: Fixed title and description colors"
echo "✅ Pipelines page: Added friendly 401 error message"
echo "✅ Pipelines page: Fixed text colors"
echo "✅ Runs page: Fixed text visibility"
echo ""
echo "==================================="
echo "Notes:"
echo "==================================="
echo "- If Pipelines page is blank, check browser console (F12)"
echo "- Clear browser cache if old styles persist (Ctrl+Shift+R)"
echo "- Verify you're logged in before accessing modules"
echo ""
