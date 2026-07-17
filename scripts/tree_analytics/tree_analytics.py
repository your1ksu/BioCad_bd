import argparse
import json
import os
import sys

import matplotlib.pyplot as plt
import seaborn as sns

# ============ ПУТИ ============
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# paths.py лежит на уровень выше (в scripts/), а не рядом с этим файлом
SCRIPTS_DIR = os.path.dirname(BASE_DIR)
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from paths import get_paths  # noqa: E402

INPUT_FILENAME = "clades_report.json"
OUTPUT_FILENAME = "tree_analytics.png"


def load_json_data(file_path):
    """Загружает данные из JSON-файла."""
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


def extract_statistics(data):
    """Извлекает по каждому гену: число клад, размер (size), глубину (depth)
    и среднюю генетическую дистанцию предок-лист (ancestor_to_leaves.mean).

    depth и ancestor_to_leaves — готовые поля из отчёта, вычислять их самим
    не нужно (в отличие от старой версии, которая пересчитывала depth сама).
    """
    genes = []
    clades_counts = []
    sizes = []
    depths = []
    ancestor_distances = []

    for gene, content in data.items():
        clades_list = content.get("mrbayes", {}).get("clades", [])

        genes.append(gene)
        clades_counts.append(len(clades_list))

        for clade in clades_list:
            sizes.append(clade.get("size", 0))
            depths.append(clade.get("depth", 0))

            ancestor_to_leaves = clade.get("ancestor_to_leaves") or {}
            mean_distance = ancestor_to_leaves.get("mean")
            if mean_distance is not None:
                ancestor_distances.append(mean_distance)

    return genes, clades_counts, sizes, depths, ancestor_distances


def plot_results(genes, clades_counts, sizes, depths, ancestor_distances, output_path):
    """Строит 4 диаграммы на основе собранных параметров и сохраняет в файл."""
    sns.set_theme(style="whitegrid")

    fig, axes = plt.subplots(1, 4, figsize=(22, 5))

    sns.barplot(x=genes, y=clades_counts, ax=axes[0], palette="viridis")
    axes[0].set_title("Количество клад по генам")
    axes[0].set_xlabel("Ген")
    axes[0].set_ylabel("Количество клад")
    axes[0].tick_params(axis="x", rotation=45)

    sns.histplot(sizes, ax=axes[1], kde=True, color="crimson", bins=10)
    axes[1].set_title("Распределение размера клады (size)")
    axes[1].set_xlabel("Размер клады (число листьев)")
    axes[1].set_ylabel("Частота")

    sns.countplot(x=depths, ax=axes[2], palette="magma")
    axes[2].set_title("Распределение глубины клады (depth)")
    axes[2].set_xlabel("Глубина")
    axes[2].set_ylabel("Количество клад")

    sns.histplot(ancestor_distances, ax=axes[3], kde=True, color="teal", bins=10)
    axes[3].set_title("Дистанция предок→листья (mean)")
    axes[3].set_xlabel("Замен на сайт (substitutions/site)")
    axes[3].set_ylabel("Частота")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"График сохранён в {output_path}")
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser(description="Анализ и визуализация статистики по кладам.")
    parser.add_argument(
        "-k", "--key",
        default="BCR",
        help="Название подпапки внутри results/ (по умолчанию: BCR)",
    )
    parser.add_argument(
        "--input-file",
        default=INPUT_FILENAME,
        help=f"Имя входного JSON-файла внутри results/<key>/ (по умолчанию: {INPUT_FILENAME})",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    paths = get_paths(args.key)

    # clades_report.json — это результат предыдущих шагов пайплайна,
    # поэтому берём его из results/<key>/, а не из data/<key>/
    input_path = os.path.join(paths["output_dir"], args.input_file)
    output_path = os.path.join(paths["output_dir"], OUTPUT_FILENAME)

    try:
        raw_data = load_json_data(input_path)
    except FileNotFoundError:
        print(f"Ошибка: не найден файл {input_path}")
        return
    except json.JSONDecodeError:
        print("Ошибка: некорректный JSON-формат.")
        return

    genes, clades_counts, sizes, depths, ancestor_distances = extract_statistics(raw_data)

    print(f"Обработано генов: {len(genes)}")
    print(f"Всего клад найдено: {len(sizes)}")
    if depths:
        print(f"Максимальная глубина клады: {max(depths)}")
    if ancestor_distances:
        print(f"Средняя дистанция предок-лист по всем кладам: "
              f"{sum(ancestor_distances) / len(ancestor_distances):.5f}")

    plot_results(genes, clades_counts, sizes, depths, ancestor_distances, output_path)


if __name__ == "__main__":
    main()
