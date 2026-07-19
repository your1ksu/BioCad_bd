#!/usr/bin/env python3
"""Build IQ-TREE trees for aligned FASTA files (ускоренная версия).

Та же команда IQ-TREE (-m MFP -B 1000), те же деревья — меняется только
планировщик, чтобы 16 ядер были загружены всё время:

  * бюджет потоков = число ядер; каждая группа берёт из бюджета столько потоков,
    сколько ей полезно (крупные — больше, мелкие — 1), через семафор. Это не
    даёт переподписки CPU и не оставляет ядра простаивать на «хвосте».
  * группы запускаются ОТ КРУПНЫХ К МЕЛКИМ (longest-processing-time-first):
    самая тяжёлая группа стартует первой и считается параллельно с длинным
    хвостом мелких, а не в самом конце.
  * группы < 4 таксонов пропускаются (IQ-TREE не строит UFBoot на <4).
"""

import argparse
import os
import shutil
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

FASTA_EXTS = {".fa", ".fasta", ".fas", ".aln"}


def discover_fasta_files(input_dir: Path) -> list[Path]:
    return sorted(p for p in input_dir.iterdir()
                  if p.is_file() and p.suffix.lower() in FASTA_EXTS)


def count_sequences(fasta_path: Path) -> int:
    with open(fasta_path) as f:
        return sum(1 for line in f if line.startswith(">"))


def detect_nproc() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        n = os.cpu_count()
        return n if n else 1


def threads_for(nseqs: int, nproc: int) -> int:
    """Сколько потоков полезно группе: мелкой хватает 1, крупной больше.
    IQ-TREE плохо масштабируется на коротких выравниваниях, поэтому потоки
    даём только по мере роста числа таксонов."""
    if nseqs >= 80:
        return min(8, nproc)
    if nseqs >= 40:
        return min(4, nproc)
    if nseqs >= 20:
        return 2
    return 1


class ThreadBudget:
    """Семафор на N потоков: задача берёт нужное число, ядра не переподписываются."""
    def __init__(self, total: int):
        self._cv = threading.Condition()
        self._free = total
        self._total = total

    def acquire(self, k: int):
        k = min(k, self._total)
        with self._cv:
            while self._free < k:
                self._cv.wait()
            self._free -= k
        return k

    def release(self, k: int):
        with self._cv:
            self._free += k
            self._cv.notify_all()


def run_iqtree(fasta: Path, name: str, subdir: Path, want_threads: int,
               budget: ThreadBudget) -> None:
    t = budget.acquire(want_threads)
    try:
        subdir.mkdir(parents=True, exist_ok=True)
        cmd = ["iqtree", "-s", str(fasta), "-m", "MFP", "-B", "1000",
               "-T", str(t), "--prefix", str(subdir / name), "-redo"]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    finally:
        budget.release(t)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build IQ-TREE trees for aligned FASTA files.")
    p.add_argument("-i", "--input", default="aligned_sequences")
    p.add_argument("-o", "--output", default="trees")
    p.add_argument("--budget", type=int, default=0,
                   help="суммарный бюджет потоков (по умолчанию: число ядер)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    input_dir, output_dir = Path(args.input), Path(args.output)
    nproc = detect_nproc()
    budget_n = args.budget if args.budget > 0 else nproc

    if not input_dir.is_dir():
        print(f"Error: directory {input_dir} not found", file=sys.stderr)
        sys.exit(1)
    if not shutil.which("iqtree"):
        print("Error: iqtree not found. Activate the pipeline environment.", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_files = discover_fasta_files(input_dir)
    if not fasta_files:
        print("No FASTA files found in", input_dir)
        sys.exit(0)

    tasks = []
    for fasta in fasta_files:
        name = fasta.stem
        nseqs = count_sequences(fasta)
        if nseqs < 4:
            print(f"[{name}] Skipping: only {nseqs} sequences (need >=4)")
            continue
        tasks.append((nseqs, fasta, name, output_dir / name))
    tasks.sort(reverse=True)                 # крупные первыми

    budget = ThreadBudget(budget_n)
    total, completed = len(tasks), 0
    # воркеров много (=число задач или ядер×несколько), реальную загрузку держит семафор
    pool_workers = min(len(tasks), nproc * 4) or 1
    print(f"IQ-TREE: {total} групп, бюджет потоков={budget_n}, крупнейшая={tasks[0][0]} таксонов")

    with ThreadPoolExecutor(max_workers=pool_workers) as ex:
        futs = {ex.submit(run_iqtree, fp, nm, sd, threads_for(n, nproc), budget): nm
                for n, fp, nm, sd in tasks}
        for fut in as_completed(futs):
            nm = futs[fut]
            try:
                fut.result()
            except Exception as e:
                print(f"[FAILED] {nm}: {e}", file=sys.stderr)
                continue
            completed += 1
            if completed % 25 == 0 or completed == total:
                print(f"[{completed}/{total}] built")

    print(f"All trees built. Results are in {output_dir}")


if __name__ == "__main__":
    main()
