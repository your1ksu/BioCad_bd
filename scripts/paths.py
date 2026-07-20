"""
Определение путей ввода/вывода по ключу.

Идея:
    - входные данные всегда лежат в  data/<key>/
    - результаты всегда сохраняются в results/<key>/
    - IMGT germline-справочники (IGHV.fasta, IGHJ.fasta и т.д.) всегда лежат
      в data/imgt/ — этот путь НЕ зависит от ключа и не меняется.

<key> — это просто название подпапки, которое передаёт пользователь
(например: имя батча, дата эксперимента, ник человека и т.п.)

Использование из командной строки:
    python paths.py -k batch1
    python paths.py --key ivanova_donor3

Использование из другого скрипта:
    from paths import get_paths
    paths = get_paths("batch1")
    input_dir = paths["input_dir"]      # .../data/batch1
    output_dir = paths["output_dir"]    # .../results/batch1
    imgt_dir = paths["imgt_dir"]        # .../data/imgt   (всегда)
"""
import argparse
import os

# Скрипт лежит в scripts/, поэтому корень репозитория — на уровень выше
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)

DATA_ROOT = os.path.join(REPO_ROOT, "data")
RESULTS_ROOT = os.path.join(REPO_ROOT, "results")

# IMGT-справочники — всегда по умолчанию тут (структура как в самом IMGT:
# вид -> локус), ключ на это не влияет и никогда не передаётся руками
IMGT_DIR = os.path.join(DATA_ROOT, "IMGT", "Homo_sapiens", "IG")


def get_paths(key, create_output=True):
    """Возвращает словарь с путями input_dir / output_dir / imgt_dir для данного ключа.

    key: str — название подпапки внутри data/ и results/
    create_output: bool — создавать ли results/<key>, если её ещё нет
    """
    if not key or not str(key).strip():
        raise ValueError("Ключ (key) не может быть пустым")

    input_dir = os.path.join(DATA_ROOT, key)
    output_dir = os.path.join(RESULTS_ROOT, key)

    if create_output:
        os.makedirs(output_dir, exist_ok=True)

    return {
        "input_dir": input_dir,
        "output_dir": output_dir,
        "imgt_dir": IMGT_DIR,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Определяет пути data/<key> и results/<key> по ключу."
    )
    parser.add_argument(
        "-k", "--key",
        required=True,
        help="Название подпапки внутри data/ и results/ (например: batch1)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    paths = get_paths(args.key)

    print(f"Входные данные:    {paths['input_dir']}")
    print(f"Результаты:        {paths['output_dir']}")
    print(f"IMGT-справочники:  {paths['imgt_dir']}  (всегда по умолчанию)")

    if not os.path.isdir(paths["input_dir"]):
        print(f"\n⚠ Внимание: папка {paths['input_dir']} ещё не существует.")
    if not os.path.isdir(paths["imgt_dir"]):
        print(f"⚠ Внимание: папка {paths['imgt_dir']} ещё не существует "
              f"(положи туда IGHV.fasta, IGHJ.fasta и т.д.)")


if __name__ == "__main__":
    main()
