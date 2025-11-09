# Security Policy

## 🔐 报告安全漏洞

如果您发现了安全漏洞，请通过私密渠道报告，**不要**在公开的 issue 中披露。

### 联系方式
- Email: security@yourdomain.com
- 或通过 GitHub Security Advisory

## ⚠️ 已知的历史安全问题

### Cloudflared Token 泄露 (2024-01)
- **状态**: 已修复
- **影响**: Cloudflared tunnel token 曾被错误提交到 git 历史
- **修复措施**:
  - 已将 `infrastructure/docker-compose.yml` 从 git 追踪中移除
  - 已添加到 `.gitignore`
  - 创建了 `.env.example` 配置模板
- **用户操作**: 
  - 如果你曾克隆过此仓库，请撤销旧的 cloudflared token
  - 使用新的 token 重新配置

## 🛡️ 安全最佳实践

### 1. 密钥管理
- ✅ **DO**: 使用 `.env` 文件存储密钥（已在 .gitignore 中）
- ✅ **DO**: 使用 `.env.example` 作为模板（不包含真实值）
- ❌ **DON'T**: 在代码中硬编码密钥
- ❌ **DON'T**: 提交 `.env` 文件到 git

### 2. Docker 配置
- ✅ **DO**: 使用 `docker-compose.example.yml` 作为模板
- ✅ **DO**: 复制到 `docker-compose.yml` 并填入真实值
- ❌ **DON'T**: 提交包含真实密钥的 `docker-compose.yml`

### 3. 环境变量
```bash
# 正确的方式
POSTGRES_PASSWORD=${POSTGRES_PASSWORD:-changeme}
CLOUDFLARED_TOKEN=${CLOUDFLARED_TOKEN}

# 错误的方式
POSTGRES_PASSWORD=my_real_password
CLOUDFLARED_TOKEN=eyJhIjoi...真实token...
```

### 4. 生产环境
- 使用强密码（至少 16 字符，包含大小写字母、数字、特殊字符）
- 定期轮换密钥
- 使用密钥管理服务（如 AWS Secrets Manager, Azure Key Vault）
- 启用访问日志和监控

## 📋 安全检查清单

### 部署前
- [ ] 所有密码已更改为强密码
- [ ] `.env` 文件不在 git 追踪中
- [ ] `docker-compose.yml` 不在 git 追踪中（或使用环境变量）
- [ ] 检查是否有硬编码的密钥
- [ ] 审查所有配置文件

### 部署后
- [ ] 更改所有默认密码
- [ ] 启用防火墙规则
- [ ] 配置 HTTPS/TLS
- [ ] 启用访问日志
- [ ] 设置监控和告警

### 定期检查
- [ ] 每月审查访问日志
- [ ] 每季度轮换密钥
- [ ] 检查依赖包的安全更新
- [ ] 审查用户权限

## 🔍 如何检查泄露

### 1. 检查 git 历史
```bash
# 搜索可能的密钥
git log -p | grep -i "password\|secret\|token\|key"

# 检查特定文件的历史
git log -p infrastructure/docker-compose.yml
```

### 2. 使用工具扫描
```bash
# 使用 gitleaks 扫描
docker run -v $(pwd):/repo zricethezav/gitleaks:latest detect --source /repo

# 使用 trufflehog
docker run --rm -v $(pwd):/repo trufflesecurity/trufflehog:latest filesystem /repo
```

### 3. GitHub 扫描
GitHub 会自动扫描推送的代码中的常见密钥模式，并发送 Secret Scanning Alert。

## 🚨 如果密钥泄露了怎么办

### 立即行动
1. **撤销泄露的密钥** - 在相应服务中撤销
2. **生成新密钥** - 创建新的密钥替换
3. **更新配置** - 在所有服务中更新为新密钥
4. **通知团队** - 告知所有相关人员

### 清理 Git 历史
```bash
# 方法 1: BFG Repo-Cleaner (推荐)
java -jar bfg.jar --delete-files sensitive-file.yml
cd repo.git
git reflog expire --expire=now --all
git gc --prune=now --aggressive
git push --force

# 方法 2: git-filter-repo
git filter-repo --path sensitive-file.yml --invert-paths
git push --force
```

⚠️ **注意**: 改写 git 历史后，所有协作者需要重新克隆仓库。

### 审计和监控
1. 检查访问日志，确认是否有未授权访问
2. 审查最近的所有操作
3. 加强监控和告警

## 📚 参考资源

- [OWASP Top 10](https://owasp.org/www-project-top-ten/)
- [GitHub Secret Scanning](https://docs.github.com/en/code-security/secret-scanning/about-secret-scanning)
- [BFG Repo-Cleaner](https://rtyley.github.io/bfg-repo-cleaner/)
- [git-filter-repo](https://github.com/newren/git-filter-repo)

---

**最后更新**: 2024-01  
**版本**: 1.0
