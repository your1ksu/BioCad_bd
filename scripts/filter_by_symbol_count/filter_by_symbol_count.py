#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from shared.utils import count_sequences, FASTA_EXTS


def filter_and_link_files(target_dir: Path, output_dir: Path,
                           min_count: int, max_count: int) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(
        p for p in target_dir.iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTS
    )

    linked = 0
    for fpath in fasta_files:
        n = count_sequences(fpath)
        if min_count < n < max_count:
            dest = output_dir / fpath.name
            if not dest.exists():
                dest.symlink_to(fpath.resolve())
            print(f"  {fpath.name}: {n} seqs → OK")
            linked += 1
        else:
            print(f"  {fpath.name}: {n} seqs → пропущен")

    print(f"\nОтфильтровано: {linked} из {len(fasta_files)} файлов")
    return linked


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Фильтрация fasta-файлов по числу последовательностей."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Папка с fasta-файлами")
    parser.add_argument("-o", "--output", required=True,
                        help="Куда сохранить отфильтрованные fasta (symlinks)")
    parser.add_argument("--min", type=int, default=5,
                        help="Минимальное число seqs (строго больше, по умолч. 5)")
    parser.add_argument("--max", type=int, default=200,
                        help="Максимальное число seqs (строго меньше, по умолч. 200)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    target_dir = Path(args.input)
    output_dir = Path(args.output)

    if not target_dir.is_dir():
        print(f"Ошибка: папка не найдена: {target_dir}", file=sys.stderr)
        sys.exit(1)

    filter_and_link_files(target_dir, output_dir, args.min, args.max)


if __name__ == "__main__":
    main()
