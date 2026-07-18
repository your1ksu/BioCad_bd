import argparse
import json
import os

import matplotlib.pyplot as plt
import seaborn as sns


OUTPUT_FILENAME = "tree_analytics.png"


def load_json_data(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


def extract_statistics(data):
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
    parser.add_argument("-i", "--input", required=True,
                        help="Путь к входному JSON-файлу (clades_report.json)")
    parser.add_argument("-o", "--output", default=OUTPUT_FILENAME,
                        help=f"Путь для сохранения графика (по умолчанию: {OUTPUT_FILENAME})")
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = args.input
    output_path = args.output

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
