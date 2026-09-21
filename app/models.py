"""请求 / 响应模型与输入校验。

所有校验在写入数据库之前完成：任一酶或序列非法都会使整个请求被 422 拒绝，
从而保证"非法输入不会留下半截保存结果"。
"""
from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, StrictInt, field_validator, model_validator

_IUPAC_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
_IUPAC_MSG = "只允许 IUPAC 碱基代码 (A C G T R Y S W K M B D H V N)"


class EnzymeIn(BaseModel):
    """一种限制性内切酶：识别序列 + 正/负链切口偏移。"""

    model_config = ConfigDict(str_strip_whitespace=True)

    name: str = Field(min_length=1, max_length=100, description="酶名称，如 EcoRI")
    pattern: str = Field(min_length=1, max_length=200, description="IUPAC 识别序列")
    cut_top: StrictInt = Field(description="正链切口相对识别区起点的偏移")
    cut_bottom: StrictInt = Field(description="负链切口相对识别区起点的偏移（正链坐标系）")

    @field_validator("pattern")
    @classmethod
    def _pattern_iupac(cls, v: str) -> str:
        v = v.upper()
        if not _IUPAC_RE.fullmatch(v):
            raise ValueError(f"识别序列非法：{_IUPAC_MSG}")
        return v

    @model_validator(mode="after")
    def _offsets_in_range(self) -> "EnzymeIn":
        limit = 2 * len(self.pattern)
        for label, value in (("cut_top", self.cut_top), ("cut_bottom", self.cut_bottom)):
            if abs(value) > limit:
                raise ValueError(
                    f"{label}={value} 非法：偏移绝对值不得超过识别序列长度的 2 倍 "
                    f"(允许范围 [-{limit}, {limit}])"
                )
        return self


class DigestRequest(BaseModel):
    sequence: str = Field(min_length=1, description="DNA 序列，仅含 IUPAC 碱基，允许空白字符分隔")
    topology: Literal["linear", "circular"]
    enzymes: list[EnzymeIn] = Field(min_length=1, max_length=64)

    @field_validator("sequence")
    @classmethod
    def _sequence_iupac(cls, v: str) -> str:
        seq = re.sub(r"\s+", "", v).upper()
        if not seq:
            raise ValueError("序列不能为空")
        if not _IUPAC_RE.fullmatch(seq):
            raise ValueError(f"序列非法：{_IUPAC_MSG}")
        return seq


# ---------- 响应模型 ----------

class CutSource(BaseModel):
    enzyme: str
    enzyme_index: int
    strand: str
    site_start: int


class EndInfo(BaseModel):
    type: str = Field(description="blunt / 5prime / 3prime / circular")
    overhang: str = Field(description="突出端序列（正链 5'->3'），平端为空串")


class CutOut(BaseModel):
    top_cut: int
    bottom_cut: int
    end_type: str
    overhang: str
    sources: list[CutSource]


class FragmentOut(BaseModel):
    index: int
    start: int
    end: int
    length: int
    spans_origin: bool
    left_end: EndInfo
    right_end: EndInfo


class EnzymeOut(BaseModel):
    name: str
    pattern: str
    cut_top: int
    cut_bottom: int


class DigestOut(BaseModel):
    id: int
    sequence: str
    topology: str
    length: int
    enzymes: list[EnzymeOut]
    cuts: list[CutOut]
    fragments: list[FragmentOut]
    created_at: str


class DigestSummary(BaseModel):
    id: int
    topology: str
    length: int
    enzyme_count: int
    fragment_count: int
    created_at: str
