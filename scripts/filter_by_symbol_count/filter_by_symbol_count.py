import argparse
import os
import shutil
import sys

# ============ ПУТИ ============
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# paths.py лежит на уровень выше (в scripts/), а не рядом с этим файлом
SCRIPTS_DIR = os.path.dirname(BASE_DIR)
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from paths import get_paths  # noqa: E402

SYMBOL = ">"
MIN_COUNT = 5   # строго больше
MAX_COUNT = 100  # строго меньше

DEFAULT_SUBFOLDERS = ["grouped_by_germlines", "v"]


def count_symbol(file_path, symbol=SYMBOL):
    """Считает, сколько раз символ встречается в файле."""
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    return content.count(symbol)


def filter_and_copy_files(target_dir, output_dir, min_count=MIN_COUNT, max_count=MAX_COUNT):
    """Копирует из target_dir в output_dir только те файлы, где количество
    символов '>' строго больше min_count и строго меньше max_count.
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"Папка для отфильтрованных файлов: {output_dir}")

    files = [
        f for f in os.listdir(target_dir)
        if os.path.isfile(os.path.join(target_dir, f))
    ]

    copied_count = 0
    for file_name in files:
        file_path = os.path.join(target_dir, file_name)
        try:
            count_symbols = count_symbol(file_path)
        except OSError as e:
            print(f"Не удалось прочитать файл {file_name}: {e}")
            continue

        if min_count < count_symbols < max_count:
            dest_path = os.path.join(output_dir, file_name)
            shutil.copy2(file_path, dest_path)
            print(f"Скопирован: {file_name} (найдено '{SYMBOL}': {count_symbols})")
            copied_count += 1

    print(f"\nГотово! Всего скопировано файлов: {copied_count}")
    return copied_count


def parse_args():
    parser = argparse.ArgumentParser(
        description="Фильтрация fasta-файлов по количеству символов '>' (числу последовательностей)."
    )
    parser.add_argument(
        "-k", "--key",
        default="BCR",
        help="Название подпапки внутри results/ (по умолчанию: BCR)",
    )
    parser.add_argument(
        "--subfolders",
        nargs="+",
        default=DEFAULT_SUBFOLDERS,
        help=(
            "Цепочка подпапок внутри results/<key>/ до папки с файлами "
            f"(по умолчанию: {' '.join(DEFAULT_SUBFOLDERS)})"
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    paths = get_paths(args.key)

    # вход берём из results/<key>/..., как и остальные скрипты пайплайна
    target_dir = os.path.join(paths["output_dir"], *args.subfolders)

    if not os.path.isdir(target_dir):
        print(f"Ошибка: путь {target_dir} не найден. Проверь --key и --subfolders.")
        return

    # новая папка создаётся рядом с целевой, с суффиксом _filtered
    last_folder_name = args.subfolders[-1] + "_filtered"
    output_dir = os.path.join(os.path.dirname(target_dir), last_folder_name)

    filter_and_copy_files(target_dir, output_dir)


if __name__ == "__main__":
    main()
