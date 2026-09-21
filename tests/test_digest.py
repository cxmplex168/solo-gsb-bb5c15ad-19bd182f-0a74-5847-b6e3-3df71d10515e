"""核心酶切逻辑单元测试。"""
from types import SimpleNamespace

from app.digest import compute_digest


def ez(name, pattern, cut_top, cut_bottom):
    return SimpleNamespace(name=name, pattern=pattern, cut_top=cut_top, cut_bottom=cut_bottom)


ECORI = ez("EcoRI", "GAATTC", 1, 5)      # 5' 突出 AATT
ECORV = ez("EcoRV", "GATATC", 3, 3)      # 平端
KPNI = ez("KpnI", "GGTACC", 5, 1)        # 3' 突出 GTAC
BSAI = ez("BsaI", "GGTCTC", 7, 11)       # Type IIS，切在识别区外


def test_ecori_linear_sticky_5prime():
    #          0123456789012
    seq = "AAAGAATTCTTTT"  # 位点起点 3，切口 top=4 bottom=8
    out = compute_digest(seq, "linear", [ECORI])
    assert len(out["cuts"]) == 1  # 回文位点正反链匹配去重为同一物理切口
    cut = out["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (4, 8)
    assert cut["end_type"] == "5prime"
    assert cut["overhang"] == "AATT"
    assert {s["strand"] for s in cut["sources"]} == {"+", "-"}

    frags = out["fragments"]
    assert [(f["start"], f["end"], f["length"]) for f in frags] == [(0, 4, 4), (4, 13, 9)]
    assert frags[0]["left_end"] == {"type": "blunt", "overhang": ""}      # 分子天然末端
    assert frags[0]["right_end"] == {"type": "5prime", "overhang": "AATT"}
    assert frags[1]["left_end"] == {"type": "5prime", "overhang": "AATT"}
    assert frags[1]["right_end"] == {"type": "blunt", "overhang": ""}


def test_ecorv_blunt():
    seq = "TTGATATCTT"  # 位点起点 2，切口 (5, 5)
    out = compute_digest(seq, "linear", [ECORV])
    assert out["cuts"][0]["end_type"] == "blunt"
    assert out["cuts"][0]["overhang"] == ""
    assert [f["length"] for f in out["fragments"]] == [5, 5]


def test_kpnI_3prime_overhang():
    #          012345678901
    seq = "AAGGTACCAA"  # 位点起点 2，切口 top=7 bottom=3
    out = compute_digest(seq, "linear", [KPNI])
    cut = out["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (7, 3)
    assert cut["end_type"] == "3prime"
    assert cut["overhang"] == "GTAC"


def test_circular_site_spanning_origin():
    seq = "AATTCG"  # 环形：从位置 5 起窗口为 G+AATTC，跨越首尾
    out = compute_digest(seq, "circular", [ECORI])
    assert len(out["cuts"]) == 1
    cut = out["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (0, 4)  # 坐标对 6 取模
    assert cut["overhang"] == "AATT"
    assert len(out["fragments"]) == 1  # 单切口开环 -> 一条全长片段
    frag = out["fragments"][0]
    assert frag["length"] == 6
    assert frag["left_end"]["type"] == "5prime"
    assert frag["right_end"]["type"] == "5prime"


def test_ambiguous_pattern_and_sequence():
    # R=A/G, Y=C/T；底物中的 N 与 Y 相交也算匹配
    # 回文识别序列 + 居中平端偏移 (3,3)：正反两个方向匹配到同一物理切口
    enz = ez("Deg", "RGATCY", 3, 3)
    out1 = compute_digest("CAGATCTG", "linear", [enz])
    out2 = compute_digest("CAGATNTG", "linear", [enz])
    for out in (out1, out2):
        assert len(out["cuts"]) == 1
        assert out["cuts"][0]["top_cut"] == 4
        assert [f["length"] for f in out["fragments"]] == [4, 4]


def test_overlapping_sites():
    out = compute_digest("AAAA", "linear", [ez("X", "AA", 1, 1)])
    assert [c["top_cut"] for c in out["cuts"]] == [1, 2, 3]  # 重叠位点全部保留
    assert [f["length"] for f in out["fragments"]] == [1, 1, 1, 1]


def test_reverse_strand_non_palindromic():
    # BsaI 识别 GGTCTC(7/11)；底物只含其反向互补 GAGACC
    #          012345678901234567
    seq = "AAAAAAGAGACCAAAAAA"  # GAGACC 起点 6
    out = compute_digest(seq, "linear", [BSAI])
    assert len(out["cuts"]) == 1
    cut = out["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (1, 5)  # s+L-cut_bottom, s+L-cut_top
    assert cut["end_type"] == "5prime"
    assert cut["sources"][0]["strand"] == "-"


def test_type_iis_cut_outside_site_on_circle_wraps():
    out = compute_digest("GGTCTC", "circular", [BSAI])
    assert len(out["cuts"]) == 1
    cut = out["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (1, 5)  # 7、11 对 6 取模
    assert cut["overhang"] == "GTCT"


def test_linear_cut_outside_molecule_ignored():
    # 线性分子末端附近的 Type IIS 切口落在分子之外 -> 不产生双链断裂
    out = compute_digest("GGTCTC", "linear", [BSAI])
    assert out["cuts"] == []
    assert len(out["fragments"]) == 1
    assert out["fragments"][0]["length"] == 6


def test_no_cut_keeps_whole_molecule():
    lin = compute_digest("AAAAAA", "linear", [ECORI])
    assert lin["cuts"] == []
    assert lin["fragments"][0]["left_end"]["type"] == "blunt"
    assert lin["fragments"][0]["length"] == 6

    cir = compute_digest("AAAAAA", "circular", [ECORI])
    assert cir["cuts"] == []
    frag = cir["fragments"][0]
    assert frag["length"] == 6
    assert frag["left_end"]["type"] == "circular"
    assert frag["right_end"]["type"] == "circular"


def test_same_physical_cut_deduped_across_enzymes():
    other = ez("EcoRI-copy", "GAATTC", 1, 5)
    out = compute_digest("AAAGAATTCTTTT", "linear", [ECORI, other])
    assert len(out["cuts"]) == 1
    assert {s["enzyme"] for s in out["cuts"][0]["sources"]} == {"EcoRI", "EcoRI-copy"}


def test_multi_enzyme_fragment_order_and_lengths():
    #          01234567890123456789
    seq = "CCGAATTCGGGATATCCC"  # EcoRI 位点起点 2，EcoRV 位点起点 10
    out = compute_digest(seq, "linear", [ECORI, ECORV])
    assert [(c["top_cut"], c["bottom_cut"]) for c in out["cuts"]] == [(3, 7), (13, 13)]
    assert [f["length"] for f in out["fragments"]] == [3, 10, 5]
    assert out["fragments"][0]["right_end"]["type"] == "5prime"
    assert out["fragments"][1]["right_end"]["type"] == "blunt"


def test_circular_two_cuts_two_fragments():
    #          01234567890123
    seq = "CGAATTCCGAATTC"  # 两个 EcoRI 位点，起点 1 和 8 -> 正链切口 2 和 9
    out = compute_digest(seq, "circular", [ECORI])
    assert len(out["cuts"]) == 2
    assert [f["length"] for f in out["fragments"]] == [7, 7]
    assert all(f["left_end"]["type"] == "5prime" for f in out["fragments"])
    assert sum(f["length"] for f in out["fragments"]) == len(seq)


def test_deterministic_output():
    seq = "CCGAATTCGGGATATCCCGGTACCAA"
    enzymes = [KPNI, ECORI, ECORV]
    first = compute_digest(seq, "linear", enzymes)
    second = compute_digest(seq, "linear", enzymes)
    assert first == second
