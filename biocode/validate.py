"""Валидация входа и выравниваний.

Никаких «тихих» пропусков: каждая проблема — это :class:`Issue` с уровнем
(``error``/``warn``/``info``), кодом и сообщением. Закрывает пункты backlog Miro:
уведомление об одной последовательности / странностях в файле, доп. проверка
длины сиквенсов в файле выравнивания.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from statistics import median

from .model import SequenceRecord

_ACGTN = set("ACGTN")
_STOPS = {"TAA", "TAG", "TGA"}


@dataclass(frozen=True)
class Issue:
    level: str          # "error" | "warn" | "info"
    code: str
    message: str
    seq_id: str | None = None

    def __str__(self) -> str:
        who = f" [{self.seq_id}]" if self.seq_id else ""
        return f"{self.level.upper()}: {self.code}{who} — {self.message}"


def _has_internal_stop(seq: str) -> bool:
    n = len(seq) - len(seq) % 3
    for i in range(0, n - 3, 3):          # без последнего (терминального) кодона
        if seq[i:i + 3] in _STOPS:
            return True
    return False


def validate_records(records: list[SequenceRecord], *, min_group_size: int = 4,
                     length_outlier_frac: float = 0.5) -> list[Issue]:
    """Проверки набора записей ДО группировки/выравнивания."""
    issues: list[Issue] = []
    if not records:
        return [Issue("error", "empty_input", "нет ни одной последовательности после фильтров")]
    if len(records) == 1:
        issues.append(Issue("warn", "single_sequence",
                            "в наборе одна последовательность — дерево построить нельзя",
                            records[0].id))

    # дубликаты id
    id_counts = Counter(r.id for r in records)
    for sid, c in id_counts.items():
        if c > 1:
            issues.append(Issue("error", "duplicate_id", f"id встречается {c} раз", sid))

    # смешанные локусы
    loci = {r.locus for r in records if r.locus}
    if len(loci) > 1:
        issues.append(Issue("warn", "mixed_loci", f"в наборе смешаны локусы: {sorted(loci)}"))

    # символы, стоп-кодоны, длины
    lengths = [r.length for r in records]
    med = median(lengths) if lengths else 0
    for r in records:
        bad = set(r.seq) - _ACGTN
        if bad:
            issues.append(Issue("warn", "non_acgt",
                                f"нестандартные символы: {sorted(bad)}", r.id))
        if _has_internal_stop(r.seq):
            issues.append(Issue("info", "internal_stop",
                                "внутренний стоп-кодон (в рамке 0)", r.id))
        if med and r.length < med * length_outlier_frac:
            issues.append(Issue("warn", "short_sequence",
                                f"длина {r.length} << медианы {med} (возможна обрезка)", r.id))
    return issues


def validate_alignment(msa: dict[str, str]) -> list[Issue]:
    """Проверка файла выравнивания: все строки одной длины (backlog Miro)."""
    issues: list[Issue] = []
    if not msa:
        return [Issue("error", "empty_alignment", "пустое выравнивание")]
    lengths = {sid: len(s) for sid, s in msa.items()}
    uniq = set(lengths.values())
    if len(uniq) > 1:
        w = Counter(lengths.values()).most_common(1)[0][0]
        for sid, L in lengths.items():
            if L != w:
                issues.append(Issue("error", "align_length_mismatch",
                                    f"длина {L} ≠ ширины выравнивания {w}", sid))
    return issues


def has_errors(issues: list[Issue]) -> bool:
    return any(i.level == "error" for i in issues)
