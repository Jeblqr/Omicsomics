#!/bin/bash
# Omicsomics 项目管理脚本

set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPOSE_FILE="$PROJECT_ROOT/infrastructure/docker-compose.yml"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 打印彩色消息
print_info() { echo -e "${BLUE}ℹ ${1}${NC}"; }
print_success() { echo -e "${GREEN}✓ ${1}${NC}"; }
print_warning() { echo -e "${YELLOW}⚠ ${1}${NC}"; }
print_error() { echo -e "${RED}✗ ${1}${NC}"; }

# 检查 Docker
check_docker() {
    if ! command -v docker &> /dev/null; then
        print_error "Docker 未安装"
        exit 1
    fi
    if ! docker info &> /dev/null; then
        print_error "Docker 未运行"
        exit 1
    fi
}

# 显示使用帮助
show_help() {
    cat << EOF
Omicsomics 项目管理工具

用法: ./manage.sh [命令]

命令:
  start         启动所有服务
  stop          停止所有服务
  restart       重启所有服务
  status        查看服务状态
  logs [服务]   查看日志 (可选指定服务名)
  ps            查看容器列表
  shell [服务]  进入容器 shell
  clean         清理所有容器和数据卷 (⚠️ 危险操作)
  init          初始化数据库和MinIO
  test          运行测试
  db:migrate    运行数据库迁移
  db:status     查看数据库迁移状态
  help          显示此帮助信息

示例:
  ./manage.sh start              # 启动所有服务
  ./manage.sh logs backend       # 查看后端日志
  ./manage.sh shell backend      # 进入后端容器
  ./manage.sh stop               # 停止所有服务

EOF
}

# 启动服务
start_services() {
    print_info "启动 Omicsomics 服务..."
    cd "$PROJECT_ROOT/infrastructure"
    
    docker compose up -d
    
    print_success "服务启动中..."
    print_info "等待服务健康检查..."
    sleep 5
    
    docker compose ps
    
    echo ""
    print_success "Omicsomics 已启动！"
    echo ""
    echo "访问地址:"
    echo "  前端:      http://localhost:5173"
    echo "  后端 API:  http://localhost:8001"
    echo "  API 文档:  http://localhost:8001/docs"
    echo "  MinIO:     http://localhost:9001"
    echo ""
    echo "查看日志: ./manage.sh logs"
}

# 停止服务
stop_services() {
    print_info "停止 Omicsomics 服务..."
    cd "$PROJECT_ROOT/infrastructure"
    docker compose down
    print_success "服务已停止"
}

# 重启服务
restart_services() {
    print_info "重启 Omicsomics 服务..."
    stop_services
    sleep 2
    start_services
}

# 查看状态
show_status() {
    print_info "服务状态:"
    cd "$PROJECT_ROOT/infrastructure"
    docker compose ps
}

# 查看日志
show_logs() {
    cd "$PROJECT_ROOT/infrastructure"
    if [ -z "$1" ]; then
        docker compose logs -f
    else
        docker compose logs -f "$1"
    fi
}

# 进入容器
enter_shell() {
    if [ -z "$1" ]; then
        print_error "请指定服务名: backend, frontend, db, redis, minio"
        exit 1
    fi
    
    cd "$PROJECT_ROOT/infrastructure"
    print_info "进入 $1 容器..."
    docker compose exec "$1" sh || docker compose exec "$1" bash
}

# 清理
clean_all() {
    print_warning "⚠️  这将删除所有容器、镜像和数据卷！"
    read -p "确认继续？(yes/no): " confirm
    
    if [ "$confirm" != "yes" ]; then
        print_info "已取消"
        exit 0
    fi
    
    print_info "清理所有资源..."
    cd "$PROJECT_ROOT/infrastructure"
    docker compose down -v --rmi all
    print_success "清理完成"
}

# 初始化
init_services() {
    print_info "初始化数据库和 MinIO..."
    cd "$PROJECT_ROOT/infrastructure"
    
    # 启动基础服务
    docker compose up -d db minio redis
    
    print_info "等待服务就绪..."
    sleep 10
    
    # 运行数据库迁移
    print_info "运行数据库迁移..."
    docker compose run --rm backend alembic upgrade head
    
    print_success "初始化完成"
}

# 运行测试
run_tests() {
    print_info "运行测试..."
    cd "$PROJECT_ROOT/infrastructure"
    docker compose run --rm backend pytest tests/ -v
}

# 运行数据库迁移
db_migrate() {
    print_info "运行数据库迁移..."
    cd "$PROJECT_ROOT/infrastructure"
    docker compose run --rm backend alembic upgrade head
}

# 查看数据库迁移状态
db_status() {
    print_info "查看数据库迁移状态..."
    cd "$PROJECT_ROOT/infrastructure"
    docker compose run --rm backend alembic current
}

# 主函数
main() {
    check_docker
    
    case "${1:-help}" in
        start)
            start_services
            ;;
        stop)
            stop_services
            ;;
        restart)
            restart_services
            ;;
        status)
            show_status
            ;;
        logs)
            show_logs "$2"
            ;;
        ps)
            show_status
            ;;
        shell)
            enter_shell "$2"
            ;;
        clean)
            clean_all
            ;;
        init)
            init_services
            ;;
        test)
            run_tests
            ;;
        db:migrate)
            db_migrate
            ;;
        db:status)
            db_status
            ;;
        help|--help|-h)
            show_help
            ;;
        *)
            print_error "未知命令: $1"
            echo ""
            show_help
            exit 1
            ;;
    esac
}

main "$@"
