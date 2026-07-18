import argparse
import os
import shutil


SYMBOL = ">"
MIN_COUNT = 5
MAX_COUNT = 100



def count_symbol(file_path, symbol=SYMBOL):
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    return content.count(symbol)


def filter_and_copy_files(target_dir, output_dir, min_count=MIN_COUNT, max_count=MAX_COUNT):
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
    parser.add_argument("-i", "--input", required=True,
                        help="Папка с исходными FASTA-файлами")
    parser.add_argument("-o", "--output", required=True,
                        help="Папка для отфильтрованных FASTA-файлов")
    return parser.parse_args()


def main():
    args = parse_args()
    target_dir = args.input
    output_dir = args.output

    if not os.path.isdir(target_dir):
        print(f"Ошибка: входная папка {target_dir} не найдена.")
        return

<<<<<<< HEAD
    # новая папка создаётся на уровне results/<key>/ — то есть рядом с самой
    # папкой grouped_by_germlines, а не внутри неё рядом с vj
    last_folder_name = args.subfolders[-1] + "_filtered"
    output_dir = os.path.join(paths["output_dir"], last_folder_name)

=======
>>>>>>> main
    filter_and_copy_files(target_dir, output_dir)


if __name__ == "__main__":
    main()
