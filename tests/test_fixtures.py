#!/usr/bin/env python3
"""Тесты на статических fixture-файлах с заранее известными мутациями.

В отличие от test_synthetic.py (который раньше генерировал MSA случайно в
рантайме), здесь используются СТАТИЧЕСКИЕ fasta-файлы из tests/fixtures/ —
их можно открыть глазами и проверить: последовательность построена из
повторяющегося паттерна ACGT, поэтому любая мутация видна невооружённым
взглядом (замена буквы ломает паттерн — например "...ACGTACTTACGT..." вместо
"...ACGTACGTACGT...").

Формат fixture: <name>_aligned.fasta — точно формат, который производит
Alina/MSA_final.py (mafft --auto), т.е. прямой вход mrbayes/run_mrbayes.py.
Рядом лежит <name>.expected.json — задокументированный ground truth: корневая
последовательность, мутации каждого листа (позиция, ref, alt) и ожидаемые
уверенные клады с обоснованием.

Никакого рантайм-рандома: fixture-файлы построены один раз детерминированным
скриптом (см. git-историю tests/fixtures/), сам тест их только читает.

Только консоль — без HTML/визуализации (см. tests/visualize_tree.py, который
не импортируется отсюда).
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

FIXTURES_DIR = Path(__file__).parent / "fixtures"
# Корень проекта = родитель tests/ (там лежат mrbayes/, clades/, biocode/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent


def find_pipeline_dir(name: str) -> Path:
    """Найти папку 'mrbayes' или 'clades': <root>/<name> либо <root>/scripts/<name>
    (локальная раскладка и раскладка ветки Nikita на GitHub отличаются)."""
    for candidate in (PROJECT_ROOT / name, PROJECT_ROOT / "scripts" / name):
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        f"Не найдена папка '{name}' ни в {PROJECT_ROOT / name}, "
        f"ни в {PROJECT_ROOT / 'scripts' / name}")


def load_fixture(name: str) -> tuple[Path, dict]:
    fasta = FIXTURES_DIR / f"{name}_aligned.fasta"
    expected = FIXTURES_DIR / f"{name}.expected.json"
    if not fasta.is_file():
        raise FileNotFoundError(f"Fixture не найден: {fasta}")
    if not expected.is_file():
        raise FileNotFoundError(f"Expected-файл не найден: {expected}")
    return fasta, json.loads(expected.read_text(encoding="utf-8"))


def run_case(name: str, python_bin: str, work_dir: Path) -> bool:
    """Прогнать один fixture через реальный pipeline (run_mrbayes.py →
    confident_clades_report.py) и сверить найденные клады с expected.json."""
    print(f"\n{'=' * 70}\nFIXTURE: {name}\n{'=' * 70}")

    fasta_path, expected = load_fixture(name)
    n_taxa = len(expected["leaves"])
    print(f"Корень: {expected['root_len']} нт, таксонов: {n_taxa}")
    for leaf, info in expected["leaves"].items():
        shared = info["shared_positions"] or "—"
        private = info["private_positions"] or "—"
        print(f"  {leaf}: группа={info['shared_group']}, общие позиции={shared}, приватные={private}")

    in_dir = work_dir / f"in_{name}"
    in_dir.mkdir(parents=True, exist_ok=True)
    (in_dir / fasta_path.name).write_bytes(fasta_path.read_bytes())

    mrbayes_dir = work_dir / f"mrbayes_{name}"
    clades_dir = work_dir / f"clades_{name}"
    clades_dir.mkdir(parents=True, exist_ok=True)

    mrbayes_scripts_dir = find_pipeline_dir("mrbayes")
    clades_scripts_dir = find_pipeline_dir("clades")

    print("\n→ mrbayes/run_mrbayes.py ...")
    r = subprocess.run(
        [python_bin, str(mrbayes_scripts_dir / "run_mrbayes.py"),
         str(in_dir), "--out", str(mrbayes_dir), "--mb-ngen", "800000"],
        capture_output=True, text=True, cwd=str(PROJECT_ROOT))
    print("  " + r.stdout.strip().splitlines()[-1] if r.stdout.strip() else "  (нет вывода)")
    if r.returncode != 0:
        print(f"  ОШИБКА: {r.stderr[-500:]}", file=sys.stderr)
        return False

    print("→ clades/confident_clades_report.py ...")
    r = subprocess.run(
        [python_bin, str(clades_scripts_dir / "confident_clades_report.py"),
         "--mrbayes-dir", str(mrbayes_dir), "--posterior-min", "0.95",
         "--out", str(clades_dir / "report.json")],
        capture_output=True, text=True, cwd=str(PROJECT_ROOT))
    print("  " + (r.stdout.strip().splitlines()[0] if r.stdout.strip() else "(нет вывода)"))
    if r.returncode != 0:
        print(f"  ОШИБКА: {r.stderr[-500:]}", file=sys.stderr)
        return False

    report = json.loads((clades_dir / "report.json").read_text(encoding="utf-8"))
    group_key = next(iter(report), None)
    found = report.get(group_key, {}).get("mrbayes", {}).get("clades", []) if group_key else []
    found_sets = [frozenset(c["leaves"]) for c in found]

    print("\nСравнение с ожидаемым:")
    ok = True
    for exp in expected["expected_confident_clades"]:
        exp_set = frozenset(exp["leaves"])
        if exp_set in found_sets:
            idx = found_sets.index(exp_set)
            print(f"  ✓ {sorted(exp['leaves'])} найдена (posterior={found[idx]['posterior']}) — {exp['reason']}")
        else:
            print(f"  ✗ {sorted(exp['leaves'])} НЕ найдена (ожидалась) — {exp['reason']}")
            ok = False

    expected_sets = [frozenset(e["leaves"]) for e in expected["expected_confident_clades"]]
    for i, fs in enumerate(found_sets):
        if fs not in expected_sets:
            print(f"  ⚠ неожиданная клада {sorted(fs)} (posterior={found[i]['posterior']})")
            ok = False

    print(f"\n{'PASS' if ok else 'FAIL'}: {name}")
    return ok


def main(python_bin: str | None = None) -> int:
    python_bin = python_bin or sys.executable

    check = subprocess.run([python_bin, "-c", "import Bio"], capture_output=True)
    if check.returncode != 0:
        print("Ошибка: нужен biopython (pip install biopython)", file=sys.stderr)
        return 1
    try:
        find_pipeline_dir("mrbayes")
        find_pipeline_dir("clades")
    except FileNotFoundError as e:
        print(f"Ошибка: {e}", file=sys.stderr)
        return 1

    fixture_names = sorted(p.stem.replace("_aligned", "")
                          for p in FIXTURES_DIR.glob("*_aligned.fasta"))
    if not fixture_names:
        print(f"Fixture-файлы не найдены в {FIXTURES_DIR}", file=sys.stderr)
        return 1

    print("=" * 70)
    print(f"ТЕСТЫ НА СТАТИЧЕСКИХ FIXTURE (заранее известные мутации)")
    print(f"Fixtures: {len(fixture_names)} — {', '.join(fixture_names)}")
    print("=" * 70)

    work_dir = Path(__file__).parent / "_fixture_runs"
    work_dir.mkdir(exist_ok=True)

    results = {name: run_case(name, python_bin, work_dir) for name in fixture_names}

    print("\n" + "=" * 70)
    print("ИТОГ")
    print("=" * 70)
    n_ok = sum(results.values())
    for name, ok in results.items():
        print(f"  {'✓' if ok else '✗'} {name}")
    print(f"\nПройдено: {n_ok}/{len(results)}")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    py = sys.argv[1] if len(sys.argv) > 1 else None
    raise SystemExit(main(py))
