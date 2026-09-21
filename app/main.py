"""FastAPI 应用入口。"""
from __future__ import annotations

import json
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException

from . import database
from .digest import compute_digest
from .models import DigestOut, DigestRequest, DigestSummary, EnzymeOut


@asynccontextmanager
async def lifespan(_: FastAPI):
    database.init_db()  # 容器启动即完成数据库初始化
    yield


app = FastAPI(
    title="DNA Restriction Digest Simulator",
    version="1.0.0",
    description="DNA 限制性酶切模拟：支持线性/环形序列、IUPAC 模糊碱基、"
                "双链切口偏移、跨原点识别、物理切口去重与片段生成。",
    lifespan=lifespan,
)


@app.get("/health", tags=["meta"])
def health() -> dict:
    return {"status": "ok"}


@app.get("/", tags=["meta"])
def root() -> dict:
    return {
        "service": "dna-restriction-digest",
        "docs": "/docs",
        "endpoints": ["POST /digests", "GET /digests", "GET /digests/{id}", "GET /health"],
    }


@app.post("/digests", response_model=DigestOut, status_code=201, tags=["digest"])
def create_digest(req: DigestRequest) -> DigestOut:
    """提交一次酶切模拟。

    所有输入先经 Pydantic 校验（序列字母表、拓扑、切口偏移范围），
    任一酶非法则整个请求以 422 拒绝，不会写入任何数据；
    计算与写库在同一事务内完成。
    """
    result = compute_digest(req.sequence, req.topology, req.enzymes)
    payload = json.dumps(
        {"cuts": result["cuts"], "fragments": result["fragments"]},
        ensure_ascii=False,
    )
    with database.get_conn() as conn:  # 事务：全部插入成功才提交
        cur = conn.execute(
            "INSERT INTO digests (sequence, topology, result_json) VALUES (?, ?, ?)",
            (req.sequence, req.topology, payload),
        )
        digest_id = cur.lastrowid
        conn.executemany(
            "INSERT INTO enzymes (digest_id, position, name, pattern, cut_top, cut_bottom)"
            " VALUES (?, ?, ?, ?, ?, ?)",
            [
                (digest_id, i, e.name, e.pattern, e.cut_top, e.cut_bottom)
                for i, e in enumerate(req.enzymes)
            ],
        )
        created_at = conn.execute(
            "SELECT created_at FROM digests WHERE id = ?", (digest_id,)
        ).fetchone()["created_at"]

    return DigestOut(
        id=digest_id,
        sequence=req.sequence,
        topology=req.topology,
        length=result["length"],
        enzymes=[EnzymeOut(**e.model_dump()) for e in req.enzymes],
        cuts=result["cuts"],
        fragments=result["fragments"],
        created_at=created_at,
    )


@app.get("/digests", response_model=list[DigestSummary], tags=["digest"])
def list_digests() -> list[DigestSummary]:
    """按提交顺序列出全部模拟记录（摘要）。"""
    with database.get_conn() as conn:
        rows = conn.execute(
            "SELECT d.id, d.topology, d.sequence, d.result_json, d.created_at,"
            "       COUNT(e.id) AS enzyme_count"
            " FROM digests d LEFT JOIN enzymes e ON e.digest_id = d.id"
            " GROUP BY d.id ORDER BY d.id"
        ).fetchall()
    return [
        DigestSummary(
            id=row["id"],
            topology=row["topology"],
            length=len(row["sequence"]),
            enzyme_count=row["enzyme_count"],
            fragment_count=len(json.loads(row["result_json"])["fragments"]),
            created_at=row["created_at"],
        )
        for row in rows
    ]


@app.get("/digests/{digest_id}", response_model=DigestOut, tags=["digest"])
def get_digest(digest_id: int) -> DigestOut:
    """按 ID 查询一次已保存的模拟结果。"""
    with database.get_conn() as conn:
        row = conn.execute(
            "SELECT * FROM digests WHERE id = ?", (digest_id,)
        ).fetchone()
        if row is None:
            raise HTTPException(status_code=404, detail=f"digest {digest_id} 不存在")
        enzyme_rows = conn.execute(
            "SELECT name, pattern, cut_top, cut_bottom FROM enzymes"
            " WHERE digest_id = ? ORDER BY position",
            (digest_id,),
        ).fetchall()

    payload = json.loads(row["result_json"])
    return DigestOut(
        id=row["id"],
        sequence=row["sequence"],
        topology=row["topology"],
        length=len(row["sequence"]),
        enzymes=[EnzymeOut(**dict(r)) for r in enzyme_rows],
        cuts=payload["cuts"],
        fragments=payload["fragments"],
        created_at=row["created_at"],
    )
