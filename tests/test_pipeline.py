#!/usr/bin/env python3
"""E2E-прогон mrbayes + clade_search на локальных фикстурах — ТОЛЬКО консоль.

Берёт tests/fixtures/*_aligned.fasta (детерминированные выравнивания с заранее
известными кладами), прогоняет через актуальные скрипты пайплайна и печатает
пути к результатам (.nex.con.tre, report.json). Сравнение с ground truth — в
tests/test_fixtures.py; здесь только «прошло/не прошло по возврату».

Требует conda-окружение с MrBayes (`mb`) и biopython. Визуализация (HTML)
намеренно не импортируется — см. tests/visualize_tree.py.

Раньше этот файл тянул данные из внешнего репозитория BIOCAD.bigchallenges и
вызывал mrbayes/run_mrbayes.py + clades/confident_clades_report.py — их в
текущей раскладке нет; переписан на self-contained фикстуры и нумерованные
скрипты.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
FIXTURES_DIR = Path(__file__).parent / "fixtures"
MRBAYES_SCRIPT = PROJECT_ROOT / "scripts" / "04b_build_trees_mrbayes" / "build_trees_mrbayes.py"
CLADES_SCRIPT = PROJECT_ROOT / "scripts" / "05_clade_search" / "clade_search.py"


def collect_groups(work_dir: Path, n_groups: int) -> dict[str, Path]:
    """Скопировать несколько *_aligned.fasta из фикстур в рабочую папку входа."""
    in_dir = work_dir / "aligned_sequences"
    in_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Path] = {}
    for src in sorted(FIXTURES_DIR.glob("*_aligned.fasta"))[:n_groups]:
        group_key = src.stem.replace("_aligned", "")
        dst = in_dir / src.name
        dst.write_bytes(src.read_bytes())
        result[group_key] = dst
        n_seqs = src.read_text().count(">")
        print(f"  ✓ {group_key}: {n_seqs} таксонов")
    return result


def run_test_cycle(n_groups: int = 3) -> int:
    print("=" * 70)
    print("E2E: MrBayes + уверенные клады (локальные фикстуры)")
    print("=" * 70)

    for script in (MRBAYES_SCRIPT, CLADES_SCRIPT):
        if not script.is_file():
            print(f"Ошибка: не найден скрипт {script}", file=sys.stderr)
            return 1
    if not shutil.which("mb"):
        print("Ошибка: не найден MrBayes (`mb`). Активируйте окружение пайплайна.",
              file=sys.stderr)
        return 1

    work_dir = Path(__file__).parent / "_pipeline_run"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    aligned_in = work_dir / "aligned_sequences"
    mrbayes_out = work_dir / "mrbayes"
    report_path = work_dir / "groups" / "report.json"

    print(f"\n1. Загрузка {n_groups} групп из фикстур...")
    groups = collect_groups(work_dir, n_groups)
    if not groups:
        print("Ошибка: фикстуры не найдены", file=sys.stderr)
        return 1

    print("\n2. build_trees_mrbayes.py ...")
    r = subprocess.run(
        [sys.executable, str(MRBAYES_SCRIPT), str(aligned_in),
         "--out", str(mrbayes_out), "--mb-ngen", "800000"],
        capture_output=True, text=True)
    print(r.stdout)
    if r.returncode != 0:
        print(f"Ошибка mrbayes: {r.stderr}", file=sys.stderr)
        return 1

    print("3. clade_search.py ...")
    r = subprocess.run(
        [sys.executable, str(CLADES_SCRIPT), "--mrbayes-dir", str(mrbayes_out),
         "--posterior-min", "0.95", "--out", str(report_path)],
        capture_output=True, text=True)
    print(r.stdout)
    if r.returncode != 0:
        print(f"Ошибка clade_search: {r.stderr}", file=sys.stderr)
        return 1

    print("=" * 70)
    print("РЕЗУЛЬТАТЫ:")
    print("=" * 70)
    for group_key in groups:
        con_tre = mrbayes_out / f"{group_key}.nex.con.tre"
        if con_tre.is_file():
            print(f"  {group_key}: {con_tre}")
    if report_path.is_file():
        print(f"\nОтчёт: {report_path}")
    return 0


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    raise SystemExit(run_test_cycle(n))
