#!/usr/bin/env python3
"""Где мутации (FR/CDR) и соотношение dN/dS — по выходу шага Дениса.

Вход:  mutations.tsv (выход scripts/analyze_mutations/analyze_mutations.py).
       Либо один файл, либо папка — тогда рекурсивно берутся все mutations.tsv
       (структура mutation_tables/<клада>/mutations.tsv), и каждая клада
       считается отдельно + сводка по всем.
Выход: таблица в консоль; --tsv <path> — та же таблица файлом.

Считает две вещи:

1. РАСПРЕДЕЛЕНИЕ ПО РЕГИОНАМ — сколько мутаций попало в каждый FR/CDR.
   Работает всегда.

2. dN/dS (Nei–Gojobori) — отношение частоты несинонимных замен к синонимным
   с поправкой на то, что несинонимных сайтов в коде объективно больше (~3:1).
   dN/dS > 1 — положительный отбор (ожидаем в CDR — антиген-связывающие петли),
   dN/dS < 1 — очищающий отбор (ожидаем в FR — каркас держит структуру).

   ВАЖНО — сейчас dN/dS посчитать НЕЛЬЗЯ из-за бага в analyze_mutations.py:
   там мутация записывается только если изменилась аминокислота
   (`if q_aa != s_aa`), поэтому синонимные замены в mutations.tsv не попадают
   вообще, а колонка is_silent всегда 'no'. Без синонимных нет знаменателя.
   Скрипт это детектирует и говорит прямо, а не печатает бессмысленный dN/dS.
   Однострочный фикс описан в readme.md рядом; после него dN/dS считается
   автоматически, менять этот скрипт не нужно.

Ограничение оценки сайтов: N_sites/S_sites считаются по ref_codon НАБЛЮДАЕМЫХ
мутаций — полной последовательности гермлайна в mutations.tsv нет. Это оценка
композиции по мутировавшим кодонам, а не по всему V-гену; при неравномерном
SHM (hotspot-мотивы WRCY/RGYW) она смещена. Для строгого dN/dS нужны все
кодоны гермлайна — см. readme.md.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from math import log
from pathlib import Path

CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

# порядок регионов как в IMGT-нумерации, от начала V до конца J
REGION_ORDER = ['FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT',
                'FR3-IMGT', 'CDR3-IMGT (germline)', 'CDR3-IMGT', 'FR4-IMGT']


def translate(codon: str) -> str:
    return CODON_TABLE.get(codon.upper().replace('U', 'T'), 'X')


def is_cdr(region: str) -> bool:
    return region.upper().startswith('CDR')


def codon_sites(codon: str) -> tuple[float, float] | None:
    """Nei–Gojobori: (несинонимных сайтов, синонимных сайтов) кодона, n+s=3.

    Для каждой из 3 позиций смотрим 3 возможные замены: доля тех, что не меняют
    аминокислоту, — вклад в синонимные сайты.
    """
    codon = codon.upper().replace('U', 'T')
    aa = translate(codon)
    if aa == 'X' or aa == '*':
        return None
    s = 0.0
    for i in range(3):
        for base in 'ACGT':
            if base == codon[i]:
                continue
            mutated = codon[:i] + base + codon[i + 1:]
            if translate(mutated) == aa:
                s += 1.0 / 3.0
    return 3.0 - s, s


def jukes_cantor(p: float) -> float | None:
    """Коррекция p-дистанции на множественные замены. None если не определена."""
    if p <= 0:
        return 0.0
    x = 1.0 - (4.0 / 3.0) * p
    if x <= 0:
        return None          # насыщение: дистанция не оценивается
    return -0.75 * log(x)


class RegionStats:
    def __init__(self) -> None:
        self.n_count = 0          # несинонимных замен (Nd)
        self.s_count = 0          # синонимных замен (Sd)
        self.n_sites = 0.0        # несинонимных сайтов
        self.s_sites = 0.0        # синонимных сайтов

    def add(self, ref_codon: str, silent: bool) -> None:
        if silent:
            self.s_count += 1
        else:
            self.n_count += 1
        sites = codon_sites(ref_codon)
        if sites is not None:
            self.n_sites += sites[0]
            self.s_sites += sites[1]

    @property
    def total(self) -> int:
        return self.n_count + self.s_count

    def dnds(self) -> tuple[float | None, str]:
        """(значение, пояснение). None если посчитать нельзя."""
        if self.s_count == 0:
            return None, "нет синонимных"
        if self.n_sites <= 0 or self.s_sites <= 0:
            return None, "нет сайтов"
        pn = self.n_count / self.n_sites
        ps = self.s_count / self.s_sites
        dn, ds = jukes_cantor(pn), jukes_cantor(ps)
        if dn is None or ds is None:
            return None, "насыщение"
        if ds == 0:
            return None, "dS=0"
        return dn / ds, ""


def read_mutations(path: Path) -> list[dict]:
    """mutations.tsv → строки. Пустые строки-заглушки (последовательность без
    мутаций) analyze_mutations.py пишет с пустым region — отбрасываем их."""
    rows = []
    with open(path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if not (row.get('region') or '').strip():
                continue
            rows.append(row)
    return rows


def collect(rows: list[dict]) -> tuple[dict[str, RegionStats], int, int]:
    """Строки → (статистика по регионам, число последовательностей, строк без кодона)."""
    stats: dict[str, RegionStats] = defaultdict(RegionStats)
    seqs = set()
    no_codon = 0
    for r in rows:
        seqs.add(r.get('sequence_id', ''))
        region = r['region'].strip()
        ref_codon = (r.get('ref_codon') or '').strip()
        silent = (r.get('is_silent') or '').strip().lower() == 'yes'
        if len(ref_codon) != 3:
            no_codon += 1
            continue
        stats[region].add(ref_codon, silent)
    return stats, len(seqs), no_codon


def sort_regions(regions) -> list[str]:
    known = [r for r in REGION_ORDER if r in regions]
    rest = sorted(r for r in regions if r not in REGION_ORDER)
    return known + rest


def format_report(stats: dict[str, RegionStats], n_seqs: int, title: str) -> list[str]:
    lines = []
    total_muts = sum(s.total for s in stats.values())
    lines.append("=" * 78)
    lines.append(f"{title}")
    lines.append("=" * 78)
    lines.append(f"последовательностей: {n_seqs}    мутаций: {total_muts}")
    lines.append("")
    lines.append(f"{'регион':<22} {'мутаций':>8} {'доля':>7} {'N':>5} {'S':>5} "
                 f"{'N/S':>6} {'dN/dS':>8}")
    lines.append("-" * 78)

    for region in sort_regions(stats):
        st = stats[region]
        share = st.total / total_muts * 100 if total_muts else 0.0
        omega, note = st.dnds()
        if omega is not None:
            omega_str = f"{omega:.3f}"
        else:
            omega_str = "—"
        ns = f"{st.n_count / st.s_count:.2f}" if st.s_count else "—"
        mark = "  ← CDR" if is_cdr(region) else ""
        lines.append(f"{region:<22} {st.total:>8} {share:>6.1f}% {st.n_count:>5} "
                     f"{st.s_count:>5} {ns:>6} {omega_str:>8}{mark}")

    lines.append("-" * 78)

    cdr = RegionStats()
    fr = RegionStats()
    for region, st in stats.items():
        target = cdr if is_cdr(region) else fr
        target.n_count += st.n_count
        target.s_count += st.s_count
        target.n_sites += st.n_sites
        target.s_sites += st.s_sites

    for label, st in (("ВСЕ CDR", cdr), ("ВСЕ FR", fr)):
        omega, note = st.dnds()
        omega_str = f"{omega:.3f}" if omega is not None else "—"
        ns = f"{st.n_count / st.s_count:.2f}" if st.s_count else "—"
        share = st.total / total_muts * 100 if total_muts else 0.0
        lines.append(f"{label:<22} {st.total:>8} {share:>6.1f}% {st.n_count:>5} "
                     f"{st.s_count:>5} {ns:>6} {omega_str:>8}")

    total_s = sum(s.s_count for s in stats.values())
    lines.append("")
    if total_s == 0 and total_muts > 0:
        lines.append("dN/dS НЕ ПОСЧИТАН: во входном файле нет ни одной синонимной мутации")
        lines.append("(во всех строках is_silent='no'). Это не свойство данных, а баг")
        lines.append("analyze_mutations.py: там `if q_aa != s_aa` отсекает синонимные")
        lines.append("замены ещё до записи, поэтому знаменатель dS всегда пуст.")
        lines.append("Фикс — одна строка, см. readme.md рядом с этим скриптом.")
    else:
        lines.append("dN/dS > 1 — положительный отбор (ожидается в CDR);")
        lines.append("dN/dS < 1 — очищающий отбор (ожидается в FR).")
        lines.append("N_sites/S_sites оценены по ref_codon наблюдаемых мутаций — см. docstring.")
    return lines


def write_tsv(stats: dict[str, RegionStats], path: Path) -> None:
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['region', 'mutations', 'nonsynonymous', 'synonymous',
                    'n_sites', 's_sites', 'n_over_s', 'dnds'])
        for region in sort_regions(stats):
            st = stats[region]
            omega, _ = st.dnds()
            w.writerow([
                region, st.total, st.n_count, st.s_count,
                round(st.n_sites, 3), round(st.s_sites, 3),
                round(st.n_count / st.s_count, 4) if st.s_count else '',
                round(omega, 4) if omega is not None else '',
            ])


def find_inputs(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    if path.is_dir():
        return sorted(path.rglob('mutations.tsv'))
    return []


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('input', type=Path,
                    help='mutations.tsv или папка mutation_tables/ (выход analyze_mutations.sh)')
    ap.add_argument('--tsv', type=Path, default=None, help='записать сводную таблицу в TSV')
    ap.add_argument('--per-clade', action='store_true',
                    help='печатать отчёт по каждой кладе отдельно (иначе только сводный)')
    args = ap.parse_args(argv)

    inputs = find_inputs(args.input)
    if not inputs:
        print(f"Не найдено mutations.tsv: {args.input}", file=sys.stderr)
        return 1

    all_rows: list[dict] = []
    for tsv in inputs:
        rows = read_mutations(tsv)
        all_rows.extend(rows)
        if args.per_clade and len(inputs) > 1:
            stats, n_seqs, _ = collect(rows)
            if stats:
                print("\n".join(format_report(stats, n_seqs, f"КЛАДА: {tsv.parent.name}")))
                print()

    stats, n_seqs, no_codon = collect(all_rows)
    if not stats:
        print("Мутаций не найдено (все строки пустые).", file=sys.stderr)
        return 1

    title = ("СВОДКА ПО ВСЕМ КЛАДАМ" if len(inputs) > 1
             else f"МУТАЦИИ: {inputs[0].parent.name or inputs[0].name}")
    print("\n".join(format_report(stats, n_seqs, title)))
    if no_codon:
        print(f"\n(пропущено строк без корректного ref_codon: {no_codon})")

    if args.tsv:
        write_tsv(stats, args.tsv)
        print(f"\nТаблица записана: {args.tsv}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
