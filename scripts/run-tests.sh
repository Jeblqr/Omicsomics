#!/bin/bash

# Omicsomics 测试运行脚本

set -e

# 颜色定义
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}🧪 Omicsomics 测试套件${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# 设置测试数据库URL
export TEST_DATABASE_URL="postgresql+asyncpg://postgres:postgres@localhost:5432/omicsomics_test"

cd "$(dirname "$0")/../backend"

# 检查环境
if ! micromamba env list | grep -q omicsomics-dev; then
    echo -e "${RED}❌ micromamba 环境 'omicsomics-dev' 不存在${NC}"
    exit 1
fi

# 显示选项
echo ""
echo "选择测试范围:"
echo "  1) 运行所有测试"
echo "  2) 仅运行通过的测试 (快速)"
echo "  3) 仅运行 GWAS 测试"
echo "  4) 仅运行认证测试"
echo "  5) 仅运行项目测试"
echo "  6) 自定义测试文件"
echo ""
read -p "请选择 (1-6): " choice

case $choice in
    1)
        echo -e "${YELLOW}运行所有测试...${NC}"
        micromamba run -n omicsomics-dev pytest tests/ -v --tb=short
        ;;
    2)
        echo -e "${YELLOW}运行快速测试...${NC}"
        micromamba run -n omicsomics-dev pytest tests/test_health.py tests/test_auth.py tests/test_projects.py -v
        ;;
    3)
        echo -e "${YELLOW}运行 GWAS 测试...${NC}"
        micromamba run -n omicsomics-dev pytest tests/test_gwas.py -v --tb=short
        ;;
    4)
        echo -e "${YELLOW}运行认证测试...${NC}"
        micromamba run -n omicsomics-dev pytest tests/test_auth.py -v
        ;;
    5)
        echo -e "${YELLOW}运行项目测试...${NC}"
        micromamba run -n omicsomics-dev pytest tests/test_projects.py -v
        ;;
    6)
        read -p "输入测试文件路径: " test_file
        echo -e "${YELLOW}运行 $test_file...${NC}"
        micromamba run -n omicsomics-dev pytest "$test_file" -v --tb=short
        ;;
    *)
        echo -e "${RED}无效选择${NC}"
        exit 1
        ;;
esac

echo ""
echo -e "${GREEN}✓ 测试完成${NC}"
