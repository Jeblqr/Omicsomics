#!/bin/bash
# Run integration tests for Omicsomics

set -e

echo "🧪 Omicsomics Integration Tests"
echo "================================"

# Check if services are running
echo "Checking services..."
if ! curl -s http://localhost:8001/healthz > /dev/null 2>&1; then
    echo "❌ Backend is not running!"
    echo "Please start services with: ./docker-start.sh"
    exit 1
fi

echo "✓ Services are running"
echo ""

# Set API base URL
export API_BASE_URL="http://localhost:8001"

# Run pytest
echo "Running pytest integration tests..."
python3 -m pytest tests/integration/test_end_to_end.py -v --tb=short -s

exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo ""
    echo "✅ All integration tests passed!"
else
    echo ""
    echo "❌ Some tests failed. Check output above."
fi

exit $exit_code
