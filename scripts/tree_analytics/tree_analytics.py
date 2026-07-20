import argparse
import json

import matplotlib.pyplot as plt
import seaborn as sns

OUTPUT_FILENAME = "tree_analytics.png"


def load_json_data(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


def extract_statistics(data, source="both"):
    genes = []
    clades_counts = []
    sizes = []
    depths = []
    ancestor_distances = []

    sources_to_read = []
    if source in ("both", "mrbayes"):
        sources_to_read.append("mrbayes")
    if source in ("both", "iqtree"):
        sources_to_read.append("iqtree")

    for gene, content in data.items():
        total_clades = 0
        for src in sources_to_read:
            clades_list = content.get(src, {}).get("clades", [])
            total_clades += len(clades_list)
            for clade in clades_list:
                sizes.append(clade.get("size", 0))
                if src == "mrbayes":
                    depths.append(clade.get("depth", 0))
                    if clade.get("ancestor_to_leaves"):
                        ancestor_distances.append(clade["ancestor_to_leaves"].get("mean", 0))
                elif src == "iqtree":
                    if clade.get("ancestor_to_leaves"):
                        ancestor_distances.append(clade["ancestor_to_leaves"].get("mean", 0))
        genes.append(gene)
        clades_counts.append(total_clades)

    return genes, clades_counts, sizes, depths, ancestor_distances


def plot_results(genes, clades_counts, sizes, depths, ancestor_distances, output_path):
    sns.set_theme(style="whitegrid")

    fig, axes = plt.subplots(1, 4, figsize=(22, 5))

    # 589 генов не помещаются на одну ось (столбики сливаются, подписи наслаиваются),
    # поэтому показываем топ по числу клад горизонтальными столбиками — так длинные
    # имена V–J читаются без поворота.
    top_n = 20
    ranked = sorted(
        ((g, c) for g, c in zip(genes, clades_counts) if c > 0),
        key=lambda gc: gc[1], reverse=True,
    )[:top_n]
    top_genes = [g for g, _ in ranked]
    top_counts = [c for _, c in ranked]
    sns.barplot(x=top_counts, y=top_genes, hue=top_genes, ax=axes[0],
                palette="viridis", legend=False, orient="h")
    axes[0].set_title(f"Топ-{top_n} генов по числу клад (из {len(genes)})")
    axes[0].set_xlabel("Количество клад")
    axes[0].set_ylabel("Ген (V–J)")

    if sizes:
        sns.histplot(sizes, ax=axes[1], kde=True, color="crimson", bins=10)
        axes[1].set_title("Распределение размера клады (size)")
        axes[1].set_xlabel("Размер клады (число листьев)")
        axes[1].set_ylabel("Частота")

    if depths:
        # discrete-гистограмма по целочисленной оси: countplot схлопывал пропущенные
        # значения глубины (10, 14, 16…28), из-за чего одиночная глубина 29 рисовалась
        # вплотную к 15. Здесь ось числовая — разрыв виден, выброс стоит на своём месте.
        sns.histplot(x=depths, discrete=True, ax=axes[2], color="purple")
        axes[2].set_title("Распределение глубины клады (depth)")
        axes[2].set_xlabel("Глубина (внутренних узлов в кладе)")
        axes[2].set_ylabel("Количество клад")
    else:
        axes[2].text(0.5, 0.5, "Нет данных (только для MrBayes)",
                     ha="center", va="center", transform=axes[2].transAxes)
        axes[2].set_title("Глубина клады (depth)")

    if ancestor_distances:
        # переводим ветвевую дистанцию (замен/сайт, доли ~0.001–0.36) в % дивергенции —
        # стандартную и читаемую единицу: медиана ~1.6 %, хвост до ~36 %.
        ancestor_divergence_pct = [d * 100 for d in ancestor_distances]
        sns.histplot(ancestor_divergence_pct, ax=axes[3], kde=True, color="teal", bins=10)
        axes[3].set_title("Дивергенция предок→листья (mean)")
        axes[3].set_xlabel("Дивергенция от предка, % (замен на 100 сайтов)")
        axes[3].set_ylabel("Частота")
    else:
        axes[3].text(0.5, 0.5, "Нет данных (только для IQ-TREE)",
                     ha="center", va="center", transform=axes[3].transAxes)
        axes[3].set_title("Дивергенция предок→листья")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"График сохранён в {output_path}")


def parse_args():
    parser = argparse.ArgumentParser(description="Анализ и визуализация статистики по кладам.")
    parser.add_argument("-i", "--input", required=True,
                        help="Путь к входному JSON-файлу (report.json)")
    parser.add_argument("-o", "--output", default=OUTPUT_FILENAME,
                        help=f"Путь для сохранения графика (по умолчанию: {OUTPUT_FILENAME})")
    parser.add_argument("--source", choices=["both", "mrbayes", "iqtree"], default="both",
                        help="Источник данных: both (оба), mrbayes, iqtree (default: both)")
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

    genes, clades_counts, sizes, depths, ancestor_distances = extract_statistics(
        raw_data, source=args.source)

    print(f"Обработано генов: {len(genes)}")
    print(f"Всего клад найдено: {len(sizes)}")
    if ancestor_distances:
        print(f"Средняя дистанция предок-лист по всем кладам: "
              f"{sum(ancestor_distances) / len(ancestor_distances):.5f}")

    plot_results(genes, clades_counts, sizes, depths, ancestor_distances, output_path)


if __name__ == "__main__":
    main()
