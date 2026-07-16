# visualize_trees.sh

Визуализирует филогенетические деревья (IQ-TREE и MrBayes) в интерактивный HTML.

## Требования

- [Miniconda](https://docs.anaconda.com/miniconda/)

## Использование

```bash
./visualize_trees.sh <input_dir> <output_dir>
```

Примеры:
```bash
./visualize_trees.sh trees trees_visualization
./visualize_trees.sh /path/to/iqtree_output /path/to/html_output
./visualize_trees.sh mrbayes_output mrbayes_viz
```

## Что делает

- Находит все файлы деревьев рекурсивно во входной директории
- **Автоматически определяет тип** по расширению:
  - **IQ-TREE**: `.treefile` — парсит `UFBoot/SH-aLRT`
  - **MrBayes**: `.con.tre`, `.t`, `.nex`, `.tre`, `.tree` — парсит posterior probability
- Рисует дерево через toytree с цветовой кодировкой поддержки
- Сохраняет интерактивный HTML (SVG + zoom/pan)
- Сохраняет структуру подпапок: `output_dir/family_name/family_name.html`

## Поддерживаемые форматы

| Источник | Расширения | Метрика поддержки | Формат в узле |
|---|---|---|---|
| IQ-TREE | `.treefile` | UFBoot / SH-aLRT | `95/90` или `95/90/85` |
| MrBayes | `.con.tre`, `.t`, `.nex`, `.tre`, `.tree` | Posterior probability | `[posterior=0.95]` или `0.95` |

## Цвета узлов

### IQ-TREE (среднее UFBoot + SH-aLRT)

| Уровень | Цвет | Условие |
|---|---|---|
| Высокая | 🟢 зелёный | ≥ 80 |
| Средняя | 🟡 жёлтый | ≥ 50 |
| Низкая | 🔴 красный | < 50 |

### MrBayes (Posterior Probability)

| Уровень | Цвет | Условие |
|---|---|---|
| Высокая | 🟢 зелёный | ≥ 0.95 |
| Средняя | 🟡 жёлтый | ≥ 0.75 |
| Низкая | 🔴 красный | < 0.75 |

## Выходные файлы

```
output_dir/
├── family1/
│   └── family1.html
├── family2/
│   └── family2.html
└── ...
```

## Зависимости

Устанавливаются автоматически в conda-окружении `trees_building_env`:
- `toytree`
- `toyplot`

## Запуск тестов

```bash
pytest tests/test_visualize_trees/test_visualize_trees.py -v
```