import argparse
import json
import os

import matplotlib.pyplot as plt
import seaborn as sns

OUTPUT_FILENAME = "tree_analytics.png"
TOP_N_GENES = 20


def load_json_data(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


def extract_statistics(data):
    genes = []
    clades_counts = []
    sizes = []
    depths = []
    ancestor_distances = []
    posteriors = []
    ufboots = []

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

            posterior = clade.get("posterior")
            if posterior is not None:
                posteriors.append(posterior)

            ufboot = clade.get("ufboot")
            if ufboot is not None:
                ufboots.append(ufboot)

        # UFBoot берём из IQ-TREE-клад (у mrbayes-клад его нет) — для 2-й линии
        # на 5-м графике: сравнение поддержки ML (IQ-TREE) и Bayes (MrBayes)
        for clade in content.get("iqtree", {}).get("clades", []):
            ub = clade.get("ufboot")
            if ub is not None:
                ufboots.append(ub)

    return genes, clades_counts, sizes, depths, ancestor_distances, posteriors, ufboots


def plot_results(genes, clades_counts, sizes, depths, ancestor_distances, posteriors, ufboots, output_path):
    sns.set_theme(style="whitegrid")

    fig, axes = plt.subplots(1, 5, figsize=(26, 6))

    # 1. Топ-20 генов по числу клад — горизонтальные столбики, чтобы длинные
    # имена генов помещались и читались (все 500+ генов на одну ось не влезают).
    gene_pairs = sorted(zip(genes, clades_counts), key=lambda p: p[1], reverse=True)
    top = gene_pairs[:TOP_N_GENES]
    top_genes = [p[0] for p in top]
    top_counts = [p[1] for p in top]

    n = len(top_genes)
    colors = [plt.cm.viridis(i / max(1, n - 1)) for i in range(n)]  # верхний (больше клад) — тёмный
    axes[0].barh(range(n), top_counts, color=colors)
    axes[0].set_yticks(range(n))
    axes[0].set_yticklabels(top_genes, fontsize=8)
    axes[0].invert_yaxis()                        # самый частый ген — сверху
    axes[0].set_title(f"Топ-{n} генов по числу клад (из {len(genes)})")
    axes[0].set_xlabel("Количество клад")
    axes[0].set_ylabel("Ген (V–J)")

    # 2. Размер клады — почти всё до ~10-15, длинный хвост редких крупных клад
    # растягивал график, поэтому обрезаем ось.
    size_cap = 15
    axes[1].hist([s for s in sizes if s <= size_cap], bins=range(0, size_cap + 2),
                 color="crimson", edgecolor="white")
    axes[1].set_title("Сколько антител в кладе")
    axes[1].set_xlabel("Размер клады\n(антител в кладе, до 15)")
    axes[1].set_ylabel("Сколько клад")
    axes[1].set_xticks(range(0, size_cap + 1))
    axes[1].set_xlim(0, size_cap)

    # 3. Глубина клады — та же история, обрезаем по тому же принципу.
    depth_cap = 10
    depths_capped = [d for d in depths if d <= depth_cap]
    sns.countplot(x=depths_capped, ax=axes[2], palette="magma", order=range(0, depth_cap + 1))
    axes[2].set_title("Сложность родословной")
    axes[2].set_xlabel("Глубина клады\n(уровней ветвления, до 10)")
    axes[2].set_ylabel("Число клад")

    # 4. Генетическая дистанция — это доля нуклеотидов, отличных от предка.
    sns.histplot(ancestor_distances, ax=axes[3], kde=True, color="teal", bins=10)
    axes[3].set_title("Отличие потомков от предка")
    axes[3].set_xlabel("Отличие от предка\n(доля изменённых нуклеотидов)")
    axes[3].set_ylabel("Сколько клад")

    # 5. Уверенность в кладе. Если есть ufboot (IQ-TREE) — накладываем вторую
    # линию для сравнения двух методов; если пуст — показываем только posterior.
    sns.histplot(posteriors, ax=axes[4], kde=True, color="darkorange",
                 bins=10, stat="density", label="MrBayes (posterior)")
    if ufboots:
        ufboots_scaled = [u / 100 for u in ufboots]    # ufboot 0-100 → шкала 0-1
        sns.histplot(ufboots_scaled, ax=axes[4], kde=True, color="steelblue",
                     bins=10, stat="density", label="IQ-TREE (UFBoot/100)")
        axes[4].legend()
    else:
        axes[4].text(0.02, 0.95,
                     "UFBoot (IQ-TREE) в отчёте не заполнен —\nвторая линия появится, когда он будет посчитан",
                     transform=axes[4].transAxes, fontsize=8, va="top", color="gray")
    axes[4].set_title("Уверенность в кладе")
    axes[4].set_xlabel("Поддержка клады\n(0–1, ближе к 1 — надёжнее)")
    axes[4].set_ylabel("Плотность")

    fig.tight_layout(w_pad=2.0)
    plt.savefig(output_path, dpi=150)
    print(f"График сохранён в {output_path}")


def parse_args():
    parser = argparse.ArgumentParser(description="Анализ и визуализация статистики по кладам.")
    parser.add_argument("-i", "--input", required=True,
                         help="Путь к входному JSON-файлу (report.json)")
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

    genes, clades_counts, sizes, depths, ancestor_distances, posteriors, ufboots = extract_statistics(raw_data)

    print(f"Обработано генов: {len(genes)}")
    print(f"Всего клад найдено: {len(sizes)}")
    if depths:
        print(f"Максимальная глубина клады: {max(depths)}")
    if ancestor_distances:
        print(f"Средняя дистанция предок-лист по всем кладам: "
              f"{sum(ancestor_distances) / len(ancestor_distances):.5f}")
    if posteriors:
        print(f"Средний posterior по всем кладам: {sum(posteriors) / len(posteriors):.4f}")
    if not ufboots:
        print("UFBoot (IQ-TREE) в этом отчёте не заполнен — на 5-м графике будет только posterior.")

    plot_results(genes, clades_counts, sizes, depths, ancestor_distances, posteriors, ufboots, output_path)


if __name__ == "__main__":
    main()
