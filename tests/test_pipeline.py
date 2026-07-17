#!/usr/bin/env python3
"""E2E-тест pipeline на реальных данных — ТОЛЬКО консоль, без HTML/GUI.

1. Загружает несколько тестовых групп из BIOCAD.bigchallenges (реальные данные
   Ксюши/Алины: anotherpipeline/part_*/*_aligned.fasta)
2. Запускает 04b_build_trees_mrbayes/build_trees_mrbayes.py
3. Запускает 05_clade_search/clade_search.py
4. Печатает статистику и пути к результатам (.nex.con.tre, report.json)

Визуализация (HTML) намеренно НЕ импортируется и не запускается здесь — см.
tests/visualize_tree.py, отдельный скрипт для рендеринга дерева из уже
посчитанных результатов. Это разделение сделано специально: данный файл
должен работать в чистой консоли даже без biopython/HTML-кода.

Расположение mrbayes/groups в репозитории отличается между локальной копией
(<root>/04b_build_trees_mrbayes/, <root>/05_clade_search/) и веткой Nikita на GitHub
(<root>/scripts/04b_build_trees_mrbayes/, <root>/scripts/05_clade_search/) — find_pipeline_dir() ищет
оба варианта, чтобы тест работал в любой раскладке.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

# Этот файл лежит в tests/, корень проекта — на уровень выше
PROJECT_ROOT = Path(__file__).resolve().parent.parent


def find_pipeline_dir(name: str) -> Path:
    """Найти папку '04b_build_trees_mrbayes' или '05_clade_search': <root>/<name> либо <root>/scripts/<name>."""
    for candidate in (PROJECT_ROOT / name, PROJECT_ROOT / "scripts" / name):
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        f"Не найдена папка '{name}' ни в {PROJECT_ROOT / name}, "
        f"ни в {PROJECT_ROOT / 'scripts' / name}")


def download_test_groups(n_groups: int = 3) -> dict[str, Path]:
    """Загрузить несколько тестовых групп из BIOCAD.bigchallenges.

    Возвращает {group_key: fasta_path} для групп разного размера.
    """
    repo_path = PROJECT_ROOT / "BIOCAD.bigchallenges"
    if not repo_path.is_dir():
        print(f"Ошибка: репозиторий не найден: {repo_path}", file=sys.stderr)
        return {}

    anotherp = repo_path / "anotherpipeline"
    available_files = []
    for part in ["part_1", "part_2", "part_3", "part_4"]:
        part_dir = anotherp / part
        if part_dir.is_dir():
            for fasta in sorted(part_dir.glob("*_aligned.fasta"))[:20]:
                group_key = fasta.stem.replace("_aligned", "")
                available_files.append((group_key, fasta))

    in_dir = PROJECT_ROOT / "aligned_sequences"
    in_dir.mkdir(exist_ok=True)
    result = {}

    for group_key, src in available_files[:n_groups]:
        dst = in_dir / f"{group_key}_aligned.fasta"
        dst.write_bytes(src.read_bytes())
        result[group_key] = dst
        n_seqs = src.read_text().count(">")
        print(f"  ✓ {group_key}: {n_seqs} таксонов")

    return result


def run_test_cycle(n_groups: int = 3) -> int:
    """Запустить полный цикл: загрузка → MrBayes → уверенные клады. Возвращает exit code."""
    print("=" * 70)
    print("E2E ТЕСТ: nexus + MrBayes + уверенные клады (реальные данные)")
    print("=" * 70)

    try:
        mrbayes_scripts_dir = find_pipeline_dir("04b_build_trees_mrbayes")
        groups_scripts_dir = find_pipeline_dir("05_clade_search")
    except FileNotFoundError as e:
        print(f"Ошибка: {e}", file=sys.stderr)
        return 1

    print(f"\n1. Загрузка {n_groups} тестовых групп из BIOCAD.bigchallenges...")
    groups = download_test_groups(n_groups)
    if not groups:
        print("Ошибка: не удалось загрузить группы", file=sys.stderr)
        return 1

    print("\n2. Запуск 04b_build_trees_mrbayes/build_trees_mrbayes.py (это займёт ~30-50s на группу)...")
    mrbayes_out = mrbayes_scripts_dir  # результаты рядом со скриптом, как в readme
    result = subprocess.run(
        [sys.executable, str(mrbayes_scripts_dir / "build_trees_mrbayes.py"),
         str(PROJECT_ROOT / "aligned_sequences"), "--out", str(mrbayes_out),
         "--mb-ngen", "2000000"],
        capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"Ошибка mrbayes: {result.stderr}", file=sys.stderr)
        return 1

    print("3. Запуск 05_clade_search/clade_search.py...")
    report_path = groups_scripts_dir / "report.json"
    result = subprocess.run(
        [sys.executable, str(groups_scripts_dir / "clade_search.py"),
         "--mrbayes-dir", str(mrbayes_out), "--posterior-min", "0.95",
         "--out", str(report_path)],
        capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"Ошибка clades: {result.stderr}", file=sys.stderr)
        return 1

    print("=" * 70)
    print("РЕЗУЛЬТАТЫ:")
    print("=" * 70)
    for group_key in groups:
        con_tre = mrbayes_out / f"{group_key}.nex.con.tre"
        print(f"\n{group_key}:")
        if con_tre.is_file():
            print(f"  Дерево: {con_tre}")
    if report_path.is_file():
        print(f"\nОтчёт всех групп: {report_path}")

    print("\nДля визуализации дерева (опционально, отдельно):")
    print(f"  python visualize_tree.py <group_key> {mrbayes_out} --report {report_path}")
    return 0


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    raise SystemExit(run_test_cycle(n))
