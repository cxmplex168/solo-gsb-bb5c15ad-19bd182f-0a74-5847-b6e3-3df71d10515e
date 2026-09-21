"""API 层测试：校验、原子性拒绝、持久化与结果查询。"""
import pytest
from fastapi.testclient import TestClient

ECORI = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("APP_DB_PATH", str(tmp_path / "test.db"))
    from app.main import app

    with TestClient(app) as c:
        yield c


def test_health(client):
    assert client.get("/health").json() == {"status": "ok"}


def test_post_and_get_roundtrip(client):
    body = {"sequence": "AAAGAATTCTTTT", "topology": "linear", "enzymes": [ECORI]}
    resp = client.post("/digests", json=body)
    assert resp.status_code == 201
    created = resp.json()
    assert created["id"] >= 1
    assert created["length"] == 13
    assert len(created["cuts"]) == 1
    assert [f["length"] for f in created["fragments"]] == [4, 9]

    got = client.get(f"/digests/{created['id']}")
    assert got.status_code == 200
    assert got.json() == created  # 查询结果与提交响应完全一致


def test_get_missing_returns_404(client):
    assert client.get("/digests/999").status_code == 404


def test_list_digests(client):
    client.post("/digests", json={"sequence": "AAAAAA", "topology": "linear", "enzymes": [ECORI]})
    client.post("/digests", json={"sequence": "AATTCG", "topology": "circular", "enzymes": [ECORI]})
    listing = client.get("/digests").json()
    assert len(listing) == 2
    assert listing[0]["topology"] == "linear"
    assert listing[1]["fragment_count"] == 1


def test_invalid_offset_rejects_entire_save(client):
    body = {
        "sequence": "AAAGAATTCTTTT",
        "topology": "linear",
        "enzymes": [
            ECORI,  # 合法
            {"name": "Bad", "pattern": "GAATTC", "cut_top": 99, "cut_bottom": 5},  # 非法偏移
        ],
    }
    assert client.post("/digests", json=body).status_code == 422
    assert client.get("/digests").json() == []  # 整次保存被拒绝，无任何残留


@pytest.mark.parametrize("cut", ["1", 1.5, True, None])
def test_non_integer_offset_rejected(client, cut):
    body = {
        "sequence": "AAAGAATTCTTTT",
        "topology": "linear",
        "enzymes": [{"name": "X", "pattern": "GAATTC", "cut_top": cut, "cut_bottom": 5}],
    }
    assert client.post("/digests", json=body).status_code == 422
    assert client.get("/digests").json() == []


def test_invalid_sequence_rejected(client):
    body = {"sequence": "ATGCX!*", "topology": "linear", "enzymes": [ECORI]}
    assert client.post("/digests", json=body).status_code == 422
    assert client.get("/digests").json() == []


def test_invalid_pattern_rejected(client):
    body = {
        "sequence": "AAAGAATTCTTTT",
        "topology": "linear",
        "enzymes": [{"name": "X", "pattern": "GAATZ", "cut_top": 1, "cut_bottom": 1}],
    }
    assert client.post("/digests", json=body).status_code == 422


def test_invalid_topology_rejected(client):
    body = {"sequence": "AAAA", "topology": "supercoiled", "enzymes": [ECORI]}
    assert client.post("/digests", json=body).status_code == 422


def test_sequence_whitespace_and_case_normalized(client):
    body = {"sequence": " aaagaattc\nTTTT ", "topology": "linear", "enzymes": [ECORI]}
    resp = client.post("/digests", json=body)
    assert resp.status_code == 201
    assert resp.json()["sequence"] == "AAAGAATTCTTTT"


def test_deterministic_ordering_across_requests(client):
    body = {
        "sequence": "CCGAATTCGGGATATCCCGGTACCAA",
        "topology": "circular",
        "enzymes": [
            {"name": "KpnI", "pattern": "GGTACC", "cut_top": 5, "cut_bottom": 1},
            ECORI,
            {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3},
        ],
    }
    r1 = client.post("/digests", json=body).json()
    r2 = client.post("/digests", json=body).json()
    assert r1["cuts"] == r2["cuts"]
    assert r1["fragments"] == r2["fragments"]
    # 片段长度守恒
    assert sum(f["length"] for f in r1["fragments"]) == len(body["sequence"])
