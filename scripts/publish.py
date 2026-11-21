"""Publish this repository to Gitee or GitHub from pure Python (no git binary needed).

The script uses ``dulwich`` for the Git wire protocol and the platform HTTP APIs
for repository creation, so it works on machines where ``git`` is not installed
and where github.com is only partially reachable.

Examples
--------
Gitee (recommended on networks where github.com is unreachable)::

    python scripts/publish.py --platform gitee --owner YOUR_GITEE_USER
    # prompts for a personal access token (Gitee -> 设置 -> 私人令牌, scope: projects)

GitHub (device flow, needs an OAuth App client id)::

    python scripts/publish.py --platform github --client-id YOUR_OAUTH_APP_ID
    # prints an 8-character code; open https://github.com/login/device

Credentials are never written to disk by this script; they stay in memory for
the duration of the push.
"""
from __future__ import annotations

import argparse
import getpass
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DESCRIPTION = (
    "StaDE: stability-aware differential expression and functional enrichment "
    "for public transcriptomes"
)
GITEE_API = "https://gitee.com/api/v5"


def ensure_dulwich() -> None:
    try:
        import dulwich  # noqa: F401
    except ImportError:
        print("[..] installing dulwich (pure-python git)")
        subprocess.run([sys.executable, "-m", "pip", "install", "--no-cache-dir", "dulwich"], check=True)


def push_main(url: str) -> None:
    """Push ``refs/heads/main`` to ``url`` using dulwich."""
    from dulwich import porcelain
    from dulwich.repo import Repo

    repo = Repo(str(ROOT))
    porcelain.push(repo, url, b"refs/heads/main:refs/heads/main")


# --------------------------------------------------------------------------- Gitee
def gitee_check(token: str, owner: str) -> dict:
    import requests

    r = requests.get(f"{GITEE_API}/user", params={"access_token": token}, timeout=60)
    if r.status_code != 200:
        raise SystemExit(f"[error] Gitee token rejected ({r.status_code}): {r.text[:300]}")
    login = r.json().get("login")
    if login != owner:
        raise SystemExit(f"[error] token belongs to {login!r}, not {owner!r}")
    return r.json()


def gitee_create(token: str, name: str, private: bool, description: str) -> dict:
    import requests

    r = requests.post(
        f"{GITEE_API}/user/repos",
        data={
            "access_token": token,
            "name": name,
            "description": description,
            "private": "true" if private else "false",
            "auto_init": "false",
        },
        timeout=60,
    )
    if r.status_code not in (200, 201):
        raise SystemExit(f"[error] Gitee repo creation failed ({r.status_code}): {r.text[:400]}")
    return r.json()


def gitee_create_or_get(token: str, owner: str, name: str, private: bool, description: str) -> dict:
    """Create the repository, or reuse it when it already exists."""
    try:
        return gitee_create(token, name, private, description)
    except SystemExit:
        import requests

        r = requests.get(
            f"{GITEE_API}/repos/{owner}/{name}", params={"access_token": token}, timeout=60
        )
        if r.status_code == 200:
            print(f"[i] repository {owner}/{name} already exists, reusing it")
            return r.json()
        raise


# -------------------------------------------------------------------------- GitHub
def github_device_login(client_id: str, scope: str = "repo") -> str:
    import requests

    r = requests.post(
        "https://github.com/login/device/code",
        headers={"Accept": "application/json"},
        data={"client_id": client_id, "scope": scope},
        timeout=60,
    )
    r.raise_for_status()
    d = r.json()
    print("\n=== GitHub device login ===")
    print(f"1. open  {d['verification_uri']}")
    print(f"2. enter code:  {d['user_code']}")
    print(f"   (expires in {int(d['expires_in']) // 60} min)\n")
    interval = int(d.get("interval", 5))
    deadline = time.time() + int(d["expires_in"])
    while time.time() < deadline:
        time.sleep(interval)
        rr = requests.post(
            "https://github.com/login/oauth/access_token",
            headers={"Accept": "application/json"},
            data={
                "client_id": client_id,
                "device_code": d["device_code"],
                "grant_type": "urn:ietf:params:oauth:grant-type:device_code",
            },
            timeout=60,
        ).json()
        if rr.get("access_token"):
            return rr["access_token"]
        err = rr.get("error")
        if err not in ("authorization_pending", "slow_down"):
            raise SystemExit(f"[error] device flow failed: {rr}")
        if err == "slow_down":
            interval += 5
    raise SystemExit("[error] device code expired; re-run the script")


def github_create(token: str, name: str, private: bool, description: str) -> dict:
    import requests

    h = {"Authorization": f"token {token}", "Accept": "application/vnd.github+json"}
    me = requests.get("https://api.github.com/user", headers=h, timeout=60).json()
    owner = me.get("login")
    if not owner:
        raise SystemExit(f"[error] cannot read GitHub user: {me}")
    r = requests.post(
        "https://api.github.com/user/repos",
        headers=h,
        json={"name": name, "description": description, "private": private, "auto_init": False},
        timeout=60,
    )
    if r.status_code not in (200, 201):
        existing = requests.get(f"https://api.github.com/repos/{owner}/{name}", headers=h, timeout=60)
        if existing.status_code != 200:
            raise SystemExit(f"[error] GitHub repo creation failed ({r.status_code}): {r.text[:400]}")
        print(f"[i] repository {owner}/{name} already exists, reusing it")
    else:
        print(f"[ok] created https://github.com/{owner}/{name}")
    return {"owner": owner, "name": name, "url": f"https://github.com/{owner}/{name}"}


def main() -> int:
    ap = argparse.ArgumentParser(description="Publish StaDE to Gitee or GitHub")
    ap.add_argument("--platform", choices=["gitee", "github"], required=True)
    ap.add_argument("--repo", default="stade")
    ap.add_argument("--owner", default=None, help="Gitee username (required for gitee)")
    ap.add_argument("--private", action="store_true")
    ap.add_argument("--description", default=DESCRIPTION)
    ap.add_argument("--client-id", default="", help="GitHub OAuth App client id (device flow)")
    ap.add_argument("--token", default="", help="personal access token (GitHub or Gitee); avoids the prompt")
    args = ap.parse_args()

    ensure_dulwich()
    if not (ROOT / ".git").exists():
        raise SystemExit("[error] .git not found; run this from the repository")

    if args.platform == "gitee":
        if not args.owner:
            raise SystemExit("[error] --owner is required for Gitee")
        token = args.token.strip()
        if not token:
            print(
                "Create a Gitee personal access token first:\n"
                "  Gitee -> 设置 -> 私人令牌 -> 生成新令牌 (scope: projects)\n"
                "It is used once, in memory, for repository creation and the push."
            )
            token = getpass.getpass("Gitee personal access token (input hidden): ").strip()
        if not token:
            raise SystemExit("[error] empty token")
        gitee_check(token, args.owner)
        info = gitee_create_or_get(token, args.owner, args.repo, args.private, args.description)
        print(f"[ok] created {info.get('html_url', args.repo)}")
        push_main(f"https://{args.owner}:{token}@gitee.com/{args.owner}/{args.repo}.git")
        url = f"https://gitee.com/{args.owner}/{args.repo}"
        print("\n" + "=" * 60)
        print("PUBLISHED:", url)
        print("=" * 60)
        return 0

    # ---- GitHub --------------------------------------------------------
    token = args.token.strip()
    if not token and args.client_id:
        token = github_device_login(args.client_id)
    if not token:
        print(
            "GitHub 需要先有一个 Personal Access Token（或 OAuth App 的 client id）。\n"
            "\n方式一（推荐，2 分钟）：\n"
            "  1) 打开 https://github.com/settings/tokens\n"
            "  2) Generate new token (classic) -> 勾选 repo -> 生成并复制\n"
            "  3) 运行:  python scripts/publish.py --platform github --token <你的token>\n"
            "\n方式二（设备码登录）：\n"
            "  1) 打开 https://github.com/settings/developers -> New OAuth App\n"
            "     Homepage 随便填，Callback 填 http://localhost\n"
            "  2) 运行:  python scripts/publish.py --platform github --client-id <ClientID>\n"
            "\n也可以发布到 Gitee: --platform gitee --owner <用户名>"
        )
        return 2

    out = github_create(token, args.repo, args.private, args.description)
    push_main(f"https://{out['owner']}:{token}@github.com/{out['owner']}/{out['name']}.git")
    print("\n" + "=" * 60)
    print("PUBLISHED:", out["url"])
    print("=" * 60)
    print("公开仓库链接（可直接发给老师）:", out["url"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
