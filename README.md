# DNA 限制性酶切模拟 API

基于 **Python + FastAPI + SQLite** 的限制性内切酶酶切模拟服务。提交一段线性或环形
DNA 序列（仅含 IUPAC 碱基）以及若干酶的识别序列与正/负链切口偏移，服务会在
**正反两个方向** 查找识别位点（允许位点重叠、含模糊碱基、可跨越环形序列首尾），
对同一物理切口去重后按分子拓扑生成片段，返回每个片段的坐标、长度以及两端的
平端/黏端类型与突出端序列。

## 快速启动

只需一条命令（依赖安装在镜像构建时完成，数据库在容器启动时自动初始化）：

```bash
docker compose up --build
```

启动后：

| 地址 | 说明 |
| --- | --- |
| `http://localhost:8000` | API 根路径（服务信息） |
| `http://localhost:8000/docs` | Swagger 交互式文档 |
| `http://localhost:8000/health` | 健康检查 |

停止：`docker compose down`（加 `-v` 可同时删除数据卷）。

## 配置

通过环境变量配置（`docker-compose.yml` 中已给出默认值）：

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `APP_PORT` | `8000` | 容器内监听端口；宿主机端口在 compose 的 `ports` 中映射（默认 `8000:8000`） |
| `APP_HOST` | `0.0.0.0` | 监听地址 |
| `APP_DB_PATH` | `/data/digest.db`（容器）/ `./data/digest.db`（本地） | SQLite 数据库文件路径 |

SQLite 数据通过 compose 命名卷 `digest-data` 持久化，重启容器数据不丢失。

## API 一览

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `POST` | `/digests` | 提交一次酶切模拟，计算并保存，返回完整结果（201） |
| `GET` | `/digests` | 列出全部模拟记录（摘要，按提交顺序） |
| `GET` | `/digests/{id}` | 按 ID 查询已保存的完整结果（不存在返回 404） |
| `GET` | `/health` | 健康检查 |

### 提交示例

```bash
curl -X POST http://localhost:8000/digests \
  -H 'Content-Type: application/json' \
  -d '{
    "sequence": "CCGAATTCGGGATATCCCGGTACCAA",
    "topology": "circular",
    "enzymes": [
      {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5},
      {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3},
      {"name": "KpnI",  "pattern": "GGTACC", "cut_top": 5, "cut_bottom": 1}
    ]
  }'
```

响应（节选）：

```json
{
  "id": 1,
  "topology": "circular",
  "length": 26,
  "cuts": [
    {"top_cut": 3,  "bottom_cut": 7,  "end_type": "5prime", "overhang": "AATT", "sources": [...]},
    {"top_cut": 13, "bottom_cut": 13, "end_type": "blunt",  "overhang": "",     "sources": [...]},
    {"top_cut": 23, "bottom_cut": 19, "end_type": "3prime", "overhang": "GTAC", "sources": [...]}
  ],
  "fragments": [
    {"index": 0, "start": 3,  "end": 13, "length": 10, "spans_origin": false,
     "left_end": {"type": "5prime", "overhang": "AATT"}, "right_end": {"type": "blunt", "overhang": ""}},
    {"index": 1, "start": 13, "end": 23, "length": 10, "spans_origin": false, "...": "..."},
    {"index": 2, "start": 23, "end": 3,  "length": 6,  "spans_origin": true,  "...": "..."}
  ]
}
```

### 查询示例

```bash
curl http://localhost:8000/digests/1     # 完整结果
curl http://localhost:8000/digests       # 摘要列表
```

## 坐标与建模约定

- **序列坐标**：0-based 半开区间 `[start, end)`；切口坐标 `c` 表示切在第 `c-1` 位与
  第 `c` 位碱基之间（即切口左侧的碱基数）。
- **切口偏移**：`cut_top` / `cut_bottom` 均相对于识别窗口起点 `s`（正链坐标系）：
  正链切口 = `s + cut_top`，负链切口 = `s + cut_bottom`。`0` 表示切在识别序列第一个
  碱基之前；允许超出识别区（如 BsaI 等 Type IIS 酶可写作 `cut_top=7, cut_bottom=11`）。
  反向匹配（识别序列的反向互补落在正链上）时自动换算为正链坐标。
- **末端类型**：`top == bottom` 为 `blunt`；`top < bottom` 为 `5prime`（5' 突出）；
  `top > bottom` 为 `3prime`（3' 突出）。`overhang` 为突出端序列（按正链 5'→3' 报告）。
  无切口时保留整条分子：线性分子末端记为 `blunt`，环形分子记为 `circular`。
- **环形序列**：识别窗口可跨越首尾（取模读取），切口坐标对序列长度取模；
  片段 `spans_origin=true` 表示该片段跨越坐标原点。
- **线性序列**：切口必须完整落在分子内部 `(0, n)`，落在分子之外（如末端附近的
  Type IIS 切口）不形成双链断裂，予以忽略。
- **模糊碱基**：识别序列与底物序列都可含 IUPAC 简并碱基
  （`R Y S W K M B D H V N`），按"集合相交"判定匹配。
- **去重**：正、负链切口坐标都相同的切口视为同一物理切口，只保留一条，
  所有命中它的酶/链/位点记录在 `sources` 中（回文位点的双向匹配、
  不同酶产生相同切口都会被合并）。
- **确定性**：切口按 `(top_cut, bottom_cut)` 排序，片段按坐标顺序生成，
  `sources` 按酶提交顺序追加——相同输入的输出顺序必然一致。

## 输入校验与原子性

以下情况整个请求以 **422** 拒绝，且不会写入任何数据（整次保存被回绝）：

- 序列含非 IUPAC 字符、为空，或 `topology` 不是 `linear` / `circular`；
- 酶的识别序列含非 IUPAC 字符；
- 切口偏移不是整数（字符串、浮点、布尔一律拒绝）；
- 切口偏移绝对值超过识别序列长度的 2 倍（`|offset| ≤ 2 × len(pattern)`）。

合法请求的计算与写库在同一事务中完成，不会出现半截保存。

## 本地开发

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload          # 默认使用 ./data/digest.db
pytest                                 # 运行测试（28 个用例）
```

## 项目结构

```
app/
  main.py       # FastAPI 路由：POST/GET /digests、/health
  digest.py     # 核心逻辑：IUPAC 匹配、双向扫描、切口去重、片段生成
  models.py     # Pydantic 请求/响应模型与校验
  database.py   # SQLite 连接与建表（启动时初始化）
tests/
  test_digest.py  # 核心逻辑单元测试
  test_api.py     # API 行为测试（校验、原子拒绝、查询、确定性）
Dockerfile
docker-compose.yml
requirements.txt
```
