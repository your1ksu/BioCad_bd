"""Контракты данных BIOCODE core.

Внутри пакета данные ходят типизированными dataclass'ами, а не «сырыми»
dict/DataFrame — это делает интерфейсы явными и тестируемыми.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any


@dataclass(frozen=True)
class SequenceRecord:
    """Одна последовательность антитела с метаданными.

    ``seq`` — рабочая нуклеотидная последовательность БЕЗ гэпов (какая именно —
    задаётся ``RunConfig.working_seq``). ``regions`` — границы FR/CDR в координатах
    ``seq`` (если аннотация доступна). ``germline`` — зародышевая V(D)J этой записи
    (для построения outgroup/предка).
    """
    id: str
    seq: str
    v_gene: str = ""
    j_gene: str = ""
    d_gene: str | None = None
    v_call: str = ""
    j_call: str = ""
    d_call: str | None = None
    locus: str = ""
    productive: bool = True
    isotype: str | None = None
    regions: dict[str, tuple[int, int]] | None = None
    germline: str | None = None
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.id:
            raise ValueError("SequenceRecord.id пустой")

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass
class Group:
    """Клональная группа (общий V+J) — вход для одного дерева."""
    key: str
    v_gene: str
    j_gene: str
    records: list[SequenceRecord] = field(default_factory=list)
    outgroup: SequenceRecord | None = None  # germline-предок = корень дерева

    @property
    def size(self) -> int:
        return len(self.records)

    def all_records(self) -> list[SequenceRecord]:
        """Записи группы + outgroup (если есть) — то, что уходит в MSA."""
        recs = list(self.records)
        if self.outgroup is not None:
            recs.append(self.outgroup)
        return recs


@dataclass
class TreeResult:
    """Результат построения дерева одним методом (ML или Bayes)."""
    method: str                      # "iqtree" | "mrbayes"
    newick: str                      # дерево с опорами на ветвях
    model: str | None = None         # модель замен (напр. GTR+F+G4)
    supports: dict[str, dict] = field(default_factory=dict)   # clade_key → {ufboot,alrt|posterior}
    ancestral: dict[str, str] | None = None                   # node → предковая последовательность
    tool_version: str = ""
    log_path: str = ""
    runtime_s: float = 0.0
    outgroup: str | None = None


@dataclass(frozen=True)
class Mutation:
    """Замена на ветви parent→child (эволюционное событие)."""
    group: str
    branch: str                      # "<parent>→<child>"
    position: int                    # позиция в MSA (0-based)
    ref: str                         # исходный символ/кодон/аминокислота
    alt: str                         # производный
    kind: str = "unknown"            # "synonymous" | "replacement" | "unknown"
    region: str = "other"            # FR1..FR4 | CDR1..CDR3 | other
    support: dict[str, float] = field(default_factory=dict)
    branch_internal: bool = False    # мутация на внутренней ветви (наследуется кладой) → persistence

    def as_row(self) -> dict:
        r = asdict(self)
        r["support"] = ";".join(f"{k}={v}" for k, v in self.support.items())
        return r
