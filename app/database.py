"""SQLite 持久层（标准库 sqlite3，无外部依赖）。"""
from __future__ import annotations

import os
import sqlite3
from pathlib import Path

SCHEMA = """
CREATE TABLE IF NOT EXISTS digests (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence    TEXT NOT NULL,
    topology    TEXT NOT NULL CHECK (topology IN ('linear', 'circular')),
    result_json TEXT NOT NULL,
    created_at  TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ', 'now'))
);

CREATE TABLE IF NOT EXISTS enzymes (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    digest_id  INTEGER NOT NULL REFERENCES digests (id) ON DELETE CASCADE,
    position   INTEGER NOT NULL,
    name       TEXT NOT NULL,
    pattern    TEXT NOT NULL,
    cut_top    INTEGER NOT NULL,
    cut_bottom INTEGER NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_enzymes_digest ON enzymes (digest_id, position);
"""


def db_path() -> str:
    """每次调用时读取环境变量，便于测试与容器配置覆盖。"""
    return os.environ.get("APP_DB_PATH", str(Path.cwd() / "data" / "digest.db"))


def get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(db_path())
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def init_db() -> None:
    path = Path(db_path())
    path.parent.mkdir(parents=True, exist_ok=True)
    with get_conn() as conn:
        conn.executescript(SCHEMA)
