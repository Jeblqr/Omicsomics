#!/usr/bin/env bash
# Test API endpoints

API_BASE="http://localhost:8001/api/v1"

echo "🔍 Testing Omicsomics API Endpoints"
echo "=================================="
echo ""

# Test health endpoint
echo "1. Testing /health (if exists)..."
curl -s "$API_BASE/../health" 2>/dev/null || echo "No /health endpoint"
echo ""

# Test auth endpoints
echo "2. Testing /auth endpoints..."
echo "   - POST /auth/register (unauthenticated, should work)"
echo "   - POST /auth/login (unauthenticated, should work)"
echo ""

# Create test user and get token
echo "3. Creating test user..."
REGISTER_RESPONSE=$(curl -s -X POST "$API_BASE/auth/register" \
  -H "Content-Type: application/json" \
  -d '{"username":"testuser_'$(date +%s)'","password":"testpass123"}')
echo "   Response: $REGISTER_RESPONSE"
echo ""

# Login to get token
echo "4. Logging in..."
LOGIN_RESPONSE=$(curl -s -X POST "$API_BASE/auth/login" \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=testuser_$(date +%s)&password=testpass123" 2>/dev/null || \
  curl -s -X POST "$API_BASE/auth/login" \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=admin&password=admin123")

TOKEN=$(echo "$LOGIN_RESPONSE" | grep -o '"access_token":"[^"]*"' | cut -d'"' -f4)

if [ -z "$TOKEN" ]; then
  echo "   ❌ Failed to get token. Response: $LOGIN_RESPONSE"
  echo "   Trying with default admin credentials..."
  LOGIN_RESPONSE=$(curl -s -X POST "$API_BASE/auth/login" \
    -H "Content-Type: application/x-www-form-urlencoded" \
    -d "username=admin&password=admin")
  TOKEN=$(echo "$LOGIN_RESPONSE" | grep -o '"access_token":"[^"]*"' | cut -d'"' -f4)
fi

if [ -z "$TOKEN" ]; then
  echo "   ❌ Still no token. Cannot test authenticated endpoints."
  exit 1
fi

echo "   ✅ Got token: ${TOKEN:0:20}..."
echo ""

# Test projects endpoint
echo "5. Testing /projects endpoints..."
curl -s "$API_BASE/projects/" -H "Authorization: Bearer $TOKEN" | head -20
echo ""

# Test pipelines endpoint (NEW)
echo "6. Testing /pipelines endpoints..."
curl -s "$API_BASE/pipelines/" -H "Authorization: Bearer $TOKEN" | head -50
echo ""

# Test runs endpoint
echo "7. Testing /runs endpoints..."
curl -s "$API_BASE/runs/" -H "Authorization: Bearer $TOKEN" | head -20
echo ""

# Test data endpoint
echo "8. Testing /data endpoints..."
curl -s "$API_BASE/data/" -H "Authorization: Bearer $TOKEN" | head -20
echo ""

echo "=================================="
echo "✅ API endpoint testing complete"
