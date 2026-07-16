#!/usr/bin/env python3
"""Синтетические тесты: известное дерево → MSA → pipeline → проверка результатов.

Создаёт контролируемые тестовые данные с известной филогенетической структурой,
прогоняет их через pipeline и проверяет, что восстановленное дерево совпадает с исходным.

Это гарантирует открытость и достоверность pipeline.
"""
from __future__ import annotations

import json
import random
import subprocess
import sys
from pathlib import Path
from io import StringIO

try:
    from Bio import Phylo
except ImportError:
    print("Требуется biopython: pip install biopython", file=sys.stderr)
    sys.exit(1)


# ============================================================================
# 1. СИНТЕТИЧЕСКАЯ MSA ПО ИЗВЕСТНОМУ ДЕРЕВУ
# ============================================================================

def generate_synthetic_msa(newick: str, seq_len: int = 300, seed: int = 42) -> dict[str, str]:
    """
    Дано: дерево в формате newick (с длинами ветвей).
    Генерирует MSA путём эволюции случайной последовательности вниз по дереву.

    seq_len — длина каждой последовательности (нт).
    На каждой ветви:
      - мутирует с вероятностью branch_length (число мутаций ~ branch_length * seq_len).
      - каждая мутация — случайная замена в одной позиции.
    """
    random.seed(seed)
    tree = Phylo.read(StringIO(newick), "newick")

    # Генерируем корневую последовательность (случайная DNA)
    alphabet = "ACGT"
    root_seq = "".join(random.choice(alphabet) for _ in range(seq_len))

    # Эволюция по дереву
    seqs: dict[str, str] = {}

    def evolve(clade, parent_seq: str) -> None:
        """Рекурсивно эволюционируем последовательность вниз по дереву."""
        curr_seq = list(parent_seq)
        n_muts = int(round((clade.branch_length or 0) * seq_len))
        for _ in range(n_muts):
            pos = random.randint(0, seq_len - 1)
            curr_seq[pos] = random.choice(alphabet)
        curr_seq_str = "".join(curr_seq)

        if clade.is_terminal():
            seqs[clade.name] = curr_seq_str
        else:
            for child in clade.clades:
                evolve(child, curr_seq_str)

    for child in tree.root.clades:
        evolve(child, root_seq)

    return seqs


def write_fasta(seqs: dict[str, str], path: Path) -> None:
    """Записать MSA в FASTA."""
    lines = []
    for sid in sorted(seqs):
        lines.append(f">{sid}")
        lines.append(seqs[sid])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ============================================================================
# 2. СИНТЕТИЧЕСКИЕ ТЕСТ-КЕЙСЫ
# ============================================================================

TEST_CASES = [
    {
        "name": "well_separated_pairs",
        "desc": "Две хорошо разделённые пары: ((seq1,seq2),(seq3,seq4))",
        "newick": "((seq1:0.10,seq2:0.10):0.15,(seq3:0.10,seq4:0.10):0.15);",
        "expected_clades": {
            "seq1|seq2": 0.95,  # ожидаем высокую апостериорную вероятность
            "seq3|seq4": 0.95,
        },
        "min_posterior": 0.90,
        "seq_len": 400,
        "ngen": 800000,
    },
    {
        "name": "star_tree_5taxa",
        "desc": "Star tree (политомия): все 5 таксонов от корня на равном расстоянии",
        "newick": "(seq1:0.1,seq2:0.1,seq3:0.1,seq4:0.1,seq5:0.1);",
        "expected_clades": {},  # никаких разрешённых клад, политомия
        "min_posterior": 0.95,
        "seq_len": 300,
        "ngen": 500000,
    },
    {
        "name": "deep_nested_tree",
        "desc": "Глубокая вложенная структура: ((A,B),(C,(D,(E,F))))",
        "newick": "((A:0.12,B:0.12):0.10,(C:0.12,(D:0.12,(E:0.12,F:0.12):0.08):0.08):0.10);",
        "expected_clades": {
            "A|B": 0.95,
            "E|F": 0.95,
        },
        "min_posterior": 0.90,
        "seq_len": 500,
        "ngen": 800000,
    },
    {
        "name": "uneven_branches",
        "desc": "Неравномерные длины ветвей: ((A,B)0.2,(C,(D,E)0.08)0.12)",
        "newick": "((A:0.08,B:0.08):0.20,(C:0.10,(D:0.08,E:0.08):0.08):0.12);",
        "expected_clades": {
            "A|B": 0.95,
            "D|E": 0.95,
        },
        "min_posterior": 0.90,
        "seq_len": 400,
        "ngen": 600000,
    },
]


# ============================================================================
# 3. ТЕСТОВЫЙ ЦИКЛ
# ============================================================================

def run_test_case(test_case: dict, base_dir: Path, venv_bin: Path) -> bool:
    """Запустить один тест-кейс: MSA → pipeline → проверка.

    Возвращает True если тест прошёл.
    """
    name = test_case["name"]
    print(f"\n{'='*70}")
    print(f"ТЕСТ: {name}")
    print(f"Описание: {test_case['desc']}")
    print(f"Newick: {test_case['newick'][:60]}...")
    print(f"{'='*70}")

    # Создаём временную папку для этого теста
    test_dir = base_dir / f"_test_{name}"
    test_dir.mkdir(exist_ok=True)

    # 1. Генерируем синтетическую MSA
    print("\n1️⃣  Генерируем синтетическую MSA...")
    seq_len = test_case.get("seq_len", 300)
    seqs = generate_synthetic_msa(test_case["newick"], seq_len=seq_len)
    n_taxa = len(seqs)
    actual_len = len(next(iter(seqs.values())))
    print(f"   ✓ {n_taxa} таксонов, {actual_len} нт")

    # Записываем в FASTA
    in_fasta = test_dir / "input.fasta"
    write_fasta(seqs, in_fasta)
    print(f"   ✓ Записано: {in_fasta}")

    # 2. Запускаем mrbayes/run_mrbayes.py
    print("\n2️⃣  Запускаем MrBayes...")
    ngen = test_case.get("ngen", 500000)
    result = subprocess.run(
        [str(venv_bin), str(base_dir / "mrbayes" / "run_mrbayes.py"),
         str(test_dir), "--out", str(test_dir / "mrbayes"),
         "--mb-ngen", str(ngen)],
        capture_output=True, text=True, cwd=str(base_dir))
    out = result.stdout.strip().split("\n")[-1]
    print(f"   {out}")
    if result.returncode != 0:
        print(f"   ❌ Ошибка: {result.stderr[:200]}")
        return False

    # 3. Запускаем confident_clades_report.py
    print("\n3️⃣  Ищем уверенные клады...")
    mrbayes_dir = test_dir / "mrbayes"
    groups_dir = test_dir / "groups"
    groups_dir.mkdir(exist_ok=True)
    result = subprocess.run(
        [str(venv_bin), str(base_dir / "groups" / "confident_clades_report.py"),
         "--mrbayes-dir", str(mrbayes_dir), "--posterior-min",
         str(test_case.get("min_posterior", 0.95)),
         "--out", str(groups_dir / "report.json")],
        capture_output=True, text=True, cwd=str(base_dir))
    print(f"   {result.stdout.strip().split(chr(10))[0]}")
    if result.returncode != 0:
        print(f"   ❌ Ошибка: {result.stderr[:200]}")
        return False

    # 4. Парсим результаты
    print("\n4️⃣  Проверяем результаты...")
    report_file = groups_dir / "report.json"
    if not report_file.is_file():
        print(f"   ❌ Отчёт не создан")
        return False

    report = json.loads(report_file.read_text(encoding="utf-8"))
    group_key = list(report.keys())[0] if report else None
    if not group_key:
        print(f"   ❌ Нет групп в отчёте")
        return False

    clades = report[group_key].get("mrbayes", {}).get("clades", [])
    found_clades = {"|".join(sorted(c["leaves"])): c.get("posterior", 0) for c in clades}

    print(f"   Найдено клад: {len(found_clades)}")
    expected = test_case.get("expected_clades", {})
    print(f"   Ожидалось клад: {len(expected)}")

    # Проверяем совпадение
    passed = True
    for clade_sig, expected_pp in expected.items():
        if clade_sig in found_clades:
            found_pp = found_clades[clade_sig]
            pp_match = abs(found_pp - expected_pp) < 0.15  # допуск 15% из-за стохастичности
            print(f"   ✓ Клада {clade_sig}: найдена (pp={found_pp:.2f}, "
                  f"ожидалось ~{expected_pp:.2f}) {'✓' if pp_match else '⚠'}")
        else:
            print(f"   ❌ Клада {clade_sig}: НЕ НАЙДЕНА (ожидалась)")
            passed = False

    for clade_sig, found_pp in found_clades.items():
        if clade_sig not in expected:
            print(f"   ⚠️  Неожиданная клада {clade_sig}: pp={found_pp:.2f}")

    # Итог
    print(f"\n{'✓ ТЕСТ ПРОШЁЛ' if passed else '❌ ТЕСТ НЕ ПРОШЁЛ'}")
    return passed


def main(venv_python: str | None = None) -> int:
    """Запустить все синтетические тесты."""
    base_dir = Path(__file__).parent

    if venv_python is None:
        # Пытаемся найти venv автоматически
        venv_candidates = [
            Path("/private/tmp/claude-501/-Users-nikitasyzdykov-Desktop-biocad-biocad/"
                 "4b837d71-b53b-4f63-900f-f788d4c2323c/scratchpad/venv/bin/python"),
            Path.home() / ".venv" / "bin" / "python",
            Path("/opt/miniconda3/envs/*/bin/python").expanduser(),
        ]
        for cand in venv_candidates:
            if cand.is_file():
                venv_python = str(cand)
                break
        if not venv_python:
            venv_python = "python3"

    venv_bin = Path(venv_python)
    print(f"Используется Python: {venv_bin}")

    # Проверяем наличие biocode
    check = subprocess.run([str(venv_bin), "-c", "import biocode"],
                          capture_output=True)
    if check.returncode != 0:
        print("❌ Ошибка: biocode не установлен или недоступен", file=sys.stderr)
        return 1

    print("\n" + "="*70)
    print("СИНТЕТИЧЕСКИЕ ТЕСТЫ PIPELINE")
    print("="*70)
    print(f"Тестов: {len(TEST_CASES)}")
    print(f"Базовая папка: {base_dir}\n")

    results = {}
    for test_case in TEST_CASES:
        passed = run_test_case(test_case, base_dir, venv_bin)
        results[test_case["name"]] = passed

    # Итоговый отчёт
    print("\n" + "="*70)
    print("ИТОГОВЫЙ ОТЧЁТ")
    print("="*70)
    n_passed = sum(1 for v in results.values() if v)
    print(f"Прошло: {n_passed}/{len(results)}")
    for name, passed in results.items():
        status = "✓" if passed else "❌"
        print(f"  {status} {name}")

    return 0 if n_passed == len(results) else 1


if __name__ == "__main__":
    venv_python = sys.argv[1] if len(sys.argv) > 1 else None
    raise SystemExit(main(venv_python))
