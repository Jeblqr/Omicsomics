#!/bin/bash

# Omicsomics 快速启动脚本
# 用于快速启动开发环境

set -e

echo "🚀 启动 Omicsomics 开发环境..."

# 颜色定义
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# 检查 Docker
if ! command -v docker &> /dev/null; then
    echo -e "${RED}❌ Docker 未安装${NC}"
    exit 1
fi

# 检查 micromamba
if ! command -v micromamba &> /dev/null; then
    echo -e "${RED}❌ Micromamba 未安装${NC}"
    exit 1
fi

echo -e "${GREEN}✓${NC} 环境检查通过"

# 启动 Docker 服务
echo "📦 启动 Docker 服务..."
cd infrastructure
docker compose up -d db minio

echo "⏳ 等待数据库就绪..."
sleep 5

# 创建测试数据库
echo "🗄️  创建测试数据库..."
docker compose exec -T db psql -U postgres -c "CREATE DATABASE omicsomics_test;" 2>/dev/null || echo "数据库可能已存在"

echo -e "${GREEN}✓${NC} Docker 服务已启动"

# 检查数据库连接
echo "🔍 检查数据库连接..."
if docker compose exec -T db pg_isready -U postgres &> /dev/null; then
    echo -e "${GREEN}✓${NC} 数据库连接正常"
else
    echo -e "${RED}❌ 数据库连接失败${NC}"
    exit 1
fi

cd ..

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}✨ Omicsomics 开发环境已就绪！${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 服务状态:"
echo "  • PostgreSQL: http://localhost:5432"
echo "  • MinIO: http://localhost:9000"
echo "  • MinIO Console: http://localhost:9001"
echo ""
echo "🧪 运行测试:"
echo "  cd backend"
echo '  TEST_DATABASE_URL="postgresql+asyncpg://postgres:postgres@localhost:5432/omicsomics_test" \'
echo "    micromamba run -n omicsomics-dev pytest tests/ -v"
echo ""
echo "🚀 启动后端 (可选):"
echo "  cd backend"
echo "  micromamba run -n omicsomics-dev uvicorn app.main:app --reload"
echo ""
echo "🎨 启动前端 (可选):"
echo "  cd frontend"
echo "  npm run dev"
echo ""
echo "🛑 停止服务:"
echo "  cd infrastructure && docker compose down"
echo ""
