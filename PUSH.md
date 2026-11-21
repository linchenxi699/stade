# 发布到 GitHub / Gitee（一条命令）

仓库已在本地初始化并提交完毕（`main` 分支，含 42+ 个文件、约 2.8 MB、2 个 commit），
并且**不需要安装 git**：发布脚本用纯 Python 的 `dulwich` 直接与远端通信。

> 说明：本机网络对 `github.com` 的 443 端口不通（`codeload`、`api` 等子域可用），
> 因此**推荐优先发布到 Gitee**；GitHub 需要你的网络能访问 github.com，或你手动开启代理。

---

## 方案一：Gitee（推荐，当前网络可用）

1. 登录 Gitee，进入 **设置 → 私人令牌 → 生成新令牌**，勾选 `projects` 权限，复制令牌；
2. 在 PowerShell 执行：

```powershell
cd D:\dsh\大学课程\stade
python scripts\publish.py --platform gitee --owner 你的Gitee用户名
```

3. 按提示粘贴令牌（输入不回显）；
4. 脚本会：校验令牌 → 创建公开仓库 `stade` → 推送 `main` → 输出仓库链接。

---

## 方案二：GitHub（需先解决网络）

1. 在 <https://github.com/settings/developers> 新建一个 OAuth App：
   - Homepage URL 随便填（如 `https://example.com`）
   - Authorization callback URL 填 `http://localhost`
   - 记下 **Client ID**；
2. 执行：

```powershell
cd D:\dsh\大学课程\stade
python scripts\publish.py --platform github --client-id 你的ClientID
```

3. 按提示打开 <https://github.com/login/device>，输入屏幕上的 8 位码；
4. 脚本自动创建公开仓库并推送。

> 如果 `github.com` 仍不通：先把系统代理（如 Clash 的 7890 端口）打开，
> 然后 `set HTTPS_PROXY=http://127.0.0.1:7890` 再执行脚本。

---

## 方案三：网页导入（完全不用命令行）

1. 先把本目录打包：`Compress-Archive -Path D:\dsh\大学课程\stade\* -DestinationPath D:\stade.zip`；
2. 打开 <https://gitee.com/projects/import>，选择"导入 Git 仓库"或直接新建仓库后上传 zip；
3. 或在 GitHub 网页 **New repository → upload an existing file** 拖入解压后的文件。

---

## 发布后 30 秒检查

- [ ] 退出登录/无痕窗口打开链接，确认 README 正常渲染；
- [ ] 仓库可见性为 **Public**；
- [ ] 能看到 `stade/`、`tests/`、`docs/`、`data/demo/`、`data/raw/genesets/`；
- [ ] 点开 `data/demo/demo_counts.tsv.gz` 可下载；
- [ ] 把链接填入课程论文提交邮件正文。

## 常见问题

**Q：提示 `token does not authenticate as ...`？**
令牌与 `--owner` 不一致，或令牌权限不足（需要 `projects`）。

**Q：提示仓库已存在？**
换一个 `--repo` 名字，或先删掉远端同名仓库。

**Q：推送很慢？**
仓库只有约 2.8 MB，正常应在 1 分钟内完成；若卡住，检查网络或换 Gitee。

**Q：以后更新代码怎么推？**
改完文件后重新执行同一条命令即可（脚本会复用已有远端；若远端已存在同名仓库，
先在网页删除或改 `--repo` 名）。
