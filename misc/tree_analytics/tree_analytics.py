import json
import matplotlib.pyplot as plt
import seaborn as sns


# 1. ПОДТЯГИВАНИЕ ФАЙЛА
def load_json_data(file_path):
    """Загружает данные из JSON-файла."""
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


# 2.3 ФУНКЦИЯ ВЫЧИСЛЕНИЯ ГЛУБИНЫ И ОТНОШЕНИЙ ПРЕДОК-ПОТОМОК
def calculate_clade_depths(clades_list):
    """Вычисляет глубину вложенности каждой клады на основе пересечения листьев (leaves).

    Глубина определяется как количество предков (родительских клад) над текущей
    кладой.
    """
    if not clades_list:
        return []

    # Превращаем списки листьев в set (множества) для быстрой проверки вложений
    clades_sets = []
    for i, c in enumerate(clades_list):
        clades_sets.append({"index": i, "leaves": set(c.get("leaves", []))})

    depths = [0] * len(clades_list)

    # Ищем, сколько "родительских" клад находится строго выше текущей клады
    for i, current in enumerate(clades_sets):
        ancestors_count = 0
        for j, other in enumerate(clades_sets):
            if i != j:
                # Если листья текущей клады являются строгим подмножеством другой клады,
                # значит 'other' — это предок для 'current'
                if current["leaves"].issubset(other["leaves"]):
                    ancestors_count += 1
        depths[i] = ancestors_count

    return depths


# 2. САМИ ФУНКЦИИ (АНАЛИЗ ДАННЫХ)
def extract_statistics(data):
    """Извлекает количество клад, размеры (size) и глубину узлов."""
    genes = []
    clades_counts = []
    all_sizes = []
    all_depths = []

    for gene, content in data.items():
        mrbayes_data = content.get("mrbayes", {})
        clades_list = mrbayes_data.get("clades", [])

        # 2.2 Длина списка по ключу "clades"
        clades_count = len(clades_list)
        genes.append(gene)
        clades_counts.append(clades_count)

        # 2.1 Статистика по параметру "size"
        for clade in clades_list:
            all_sizes.append(clade.get("size", 0))

        # 2.3 Вычисление глубины внутренних узлов
        gene_depths = calculate_clade_depths(clades_list)
        all_depths.extend(gene_depths)

    return genes, clades_counts, all_sizes, all_depths


# 3. ВЫВОД И СТРОЕНИЕ ДИАГРАММ
def plot_results(genes, clades_counts, all_sizes, all_depths):
    """Строит 3 диаграммы на основе собранных параметров."""
    sns.set_theme(style="whitegrid")

    # Создаем холст для трех графиков (1 строка, 3 колонки)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # График 1: Количество клад (2.2)
    sns.barplot(x=genes, y=clades_counts, ax=axes[0], palette="viridis")
    axes[0].set_title("Количество клад (clades) по генам", fontsize=12)
    axes[0].set_xlabel("Идентификатор гена")
    axes[0].set_ylabel("Количество словарей-клад")
    axes[0].tick_params(axis="x", rotation=45)

    # График 2: Распределение размеров (2.1)
    sns.histplot(all_sizes, ax=axes[1], kde=True, color="crimson", bins=10)
    axes[1].set_title("Распределение параметра 'size'", fontsize=12)
    axes[1].set_xlabel("Размер клады (size)")
    axes[1].set_ylabel("Частота")

    # График 3: Глубина внутренних узлов / Удаленность потомков (2.3)
    sns.countplot(x=all_depths, ax=axes[2], palette="magma")
    axes[2].set_title("Глубина внутренних узлов (Уровень предка)", fontsize=12)
    axes[2].set_xlabel("Глубина (0 = корень/главный предок)")
    axes[2].set_ylabel("Количество узлов на этом уровне")

    plt.tight_layout()
    plt.show()


# Главный блок выполнения
if __name__ == "__main__":
    INPUT_FILE_PATH = "report2.json"

    try:
        # Шаг 1: Загрузка
        raw_data = load_json_data(INPUT_FILE_PATH)

        # Шаг 2: Обработка
        genes_list, counts, sizes, depths = extract_statistics(raw_data)

        # Вывод текстовой статистики
        print(f"Обработано генов: {len(genes_list)}")
        print(f"Всего клад найдено: {len(sizes)}")
        if depths:
            print(
                f"Максимальная глубина потомка в дереве: {max(depths)} уровней от предка"
            )

        # Шаг 3: Визуализация
        plot_results(genes_list, counts, sizes, depths)

    except FileNotFoundError:
        print(f"Ошибка: Не найден файл {INPUT_FILE_PATH}")
    except json.JSONDecodeError:
        print("Ошибка: Некорректный JSON-формат.")
