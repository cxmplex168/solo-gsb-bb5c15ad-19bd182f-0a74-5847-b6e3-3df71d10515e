"""核心酶切模拟逻辑。

约定
----
* 序列坐标为 0-based 半开区间；切口坐标 c 表示在第 c-1 位与第 c 位碱基之间切开
  （即坐标等于切口左侧的碱基数）。
* 酶的 cut_top / cut_bottom 为相对于识别窗口起点 s 的偏移：
  正链切口坐标 = s + cut_top，负链切口坐标 = s + cut_bottom（均为正链坐标系）。
  偏移为 0 表示切在识别序列第一个碱基之前，为 len(pattern) 表示切在最后一个碱基之后，
  允许超出识别区（如 Type IIS 酶），取值范围由 API 层校验。
* 反向匹配（识别序列的反向互补落在正链上）时，切口坐标换算为：
  top = s + L - cut_bottom，bottom = s + L - cut_top。
* 模糊碱基按 IUPAC 集合相交判定匹配（如 R 与 N 相交即视为匹配）。
* 同一物理切口（正、负链切口坐标都相同）只保留一次，来源酶列表合并。
* 线性分子：切口必须完整落在分子内部 (0, n) 才生效；环形分子：切口坐标对 n 取模。
"""
from __future__ import annotations

IUPAC_SETS = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}

BLUNT = "blunt"
OVERHANG_5 = "5prime"
OVERHANG_3 = "3prime"
CIRCULAR_END = "circular"


def reverse_complement(pattern: str) -> str:
    return "".join(_COMPLEMENT[b] for b in reversed(pattern))


def _window_matches(pattern: str, window: str) -> bool:
    """IUPAC 集合相交即匹配（允许识别序列与底物双方都含模糊碱基）。"""
    return all(IUPAC_SETS[p] & IUPAC_SETS[w] for p, w in zip(pattern, window))


def _find_raw_cuts(sequence: str, topology: str, pattern: str,
                   cut_top: int, cut_bottom: int) -> list[tuple[int, int, str, int]]:
    """枚举一条酶的全部候选切口，返回 (top_raw, bottom_raw, strand, site_start)。

    坐标尚未取模 / 过滤。环形序列允许识别窗口跨越首尾（取模读碱基），
    每个起点只扫描一次，天然允许位点相互重叠。
    """
    n = len(sequence)
    size = len(pattern)
    rc = reverse_complement(pattern)
    starts = range(n) if topology == "circular" else range(0, n - size + 1)
    raw: list[tuple[int, int, str, int]] = []
    for s in starts:
        if topology == "circular":
            window = "".join(sequence[(s + i) % n] for i in range(size))
        else:
            window = sequence[s:s + size]
        if _window_matches(pattern, window):
            raw.append((s + cut_top, s + cut_bottom, "+", s))
        if _window_matches(rc, window):
            raw.append((s + size - cut_bottom, s + size - cut_top, "-", s))
    return raw


def _end_info(sequence: str, topology: str, t_raw: int, b_raw: int) -> dict:
    """根据正/负链切口坐标推断末端类型与突出端序列（按正链 5'->3' 报告）。"""
    n = len(sequence)
    if t_raw == b_raw:
        return {"type": BLUNT, "overhang": ""}
    lo, hi = (t_raw, b_raw) if t_raw < b_raw else (b_raw, t_raw)
    if topology == "circular":
        overhang = "".join(sequence[i % n] for i in range(lo, hi))
    else:
        overhang = sequence[lo:hi]
    return {"type": OVERHANG_5 if t_raw < b_raw else OVERHANG_3, "overhang": overhang}


def _build_fragments(n: int, topology: str, cuts: list[dict]) -> list[dict]:
    """按拓扑把有序切口装配成片段。片段边界以正链切口坐标为准。"""
    if not cuts:
        # 无切口：保留整条分子
        end_type = CIRCULAR_END if topology == "circular" else BLUNT
        end = {"type": end_type, "overhang": ""}
        stop = 0 if topology == "circular" else n
        return [{
            "index": 0, "start": 0, "end": stop, "length": n,
            "spans_origin": False,
            "left_end": dict(end), "right_end": dict(end),
        }]

    fragments: list[dict] = []
    if topology == "linear":
        boundaries = [0] + [c["top_cut"] for c in cuts] + [n]
        ends = ([{"type": BLUNT, "overhang": ""}]
                + [c["end"] for c in cuts]
                + [{"type": BLUNT, "overhang": ""}])
        for i in range(len(cuts) + 1):
            start, stop = boundaries[i], boundaries[i + 1]
            fragments.append({
                "index": i, "start": start, "end": stop,
                "length": stop - start, "spans_origin": False,
                "left_end": ends[i], "right_end": ends[i + 1],
            })
    else:
        k = len(cuts)
        for i, cut in enumerate(cuts):
            nxt = cuts[(i + 1) % k]
            start, stop = cut["top_cut"], nxt["top_cut"]
            length = (stop - start) % n
            if k == 1:
                length = n  # 单切口开环，片段即整条分子
            fragments.append({
                "index": i, "start": start, "end": stop,
                "length": length,
                "spans_origin": start + length > n,
                "left_end": cut["end"], "right_end": nxt["end"],
            })
    return fragments


def compute_digest(sequence: str, topology: str, enzymes) -> dict:
    """对一条序列执行多酶联合酶切，返回 {"length", "cuts", "fragments"}。

    输出顺序完全由输入决定：切口按 (top_cut, bottom_cut) 排序，
    来源按酶提交顺序与扫描顺序追加，片段按坐标顺序生成。
    """
    n = len(sequence)
    merged: dict[tuple[int, int], dict] = {}
    for idx, ez in enumerate(enzymes):
        for t_raw, b_raw, strand, site_start in _find_raw_cuts(
            sequence, topology, ez.pattern, ez.cut_top, ez.cut_bottom
        ):
            if topology == "circular":
                key = (t_raw % n, b_raw % n)
            else:
                if not (0 < t_raw < n and 0 < b_raw < n):
                    continue  # 线性分子：切口须完整落在分子内部
                key = (t_raw, b_raw)
            cut = merged.get(key)
            if cut is None:
                cut = merged[key] = {
                    "top_cut": key[0],
                    "bottom_cut": key[1],
                    "end": _end_info(sequence, topology, t_raw, b_raw),
                    "sources": [],
                }
            cut["sources"].append({
                "enzyme": ez.name,
                "enzyme_index": idx,
                "strand": strand,
                "site_start": site_start,
            })

    cuts = sorted(merged.values(), key=lambda c: (c["top_cut"], c["bottom_cut"]))
    return {
        "length": n,
        "cuts": [
            {
                "top_cut": c["top_cut"],
                "bottom_cut": c["bottom_cut"],
                "end_type": c["end"]["type"],
                "overhang": c["end"]["overhang"],
                "sources": c["sources"],
            }
            for c in cuts
        ],
        "fragments": _build_fragments(n, topology, cuts),
    }
