"""Реконструкция предков: парсинг IQ-TREE ``<prefix>.state``.

IQ-TREE ``--asr`` пишет посайтовые предковые состояния внутренних узлов в
TSV-файл ``.state`` с шапкой ``Node  Site  State  p_A  p_C  p_G  p_T`` (после
строк-комментариев, начинающихся с ``#``). Здесь собираем по узлам полные
предковые последовательности (в координатах MSA).
"""
from __future__ import annotations

from pathlib import Path


def parse_state_file(path: str | Path) -> dict[str, str]:
    """``<prefix>.state`` → {node_name: предковая_последовательность} (в колонках MSA)."""
    path = Path(path)
    per_node: dict[str, list[tuple[int, str]]] = {}
    header: list[str] | None = None
    node_i = site_i = state_i = -1
    for line in path.read_text(errors="ignore").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = [p.strip() for p in parts]
            try:
                node_i = header.index("Node")
                site_i = header.index("Site")
                state_i = header.index("State")
            except ValueError:
                # неожиданная шапка — не гадаем, отдаём пусто
                return {}
            continue
        if len(parts) <= max(node_i, site_i, state_i):
            continue
        node = parts[node_i]
        try:
            site = int(parts[site_i])
        except ValueError:
            continue
        per_node.setdefault(node, []).append((site, parts[state_i].upper()))
    return {node: "".join(st for _, st in sorted(sites))
            for node, sites in per_node.items()}
