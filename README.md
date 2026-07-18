# BioCad_bd — филогенетический конвейер по антителам (BCR)

Конвейер обработки B-клеточных рецепторов (BCR): от сырого AIRR-подобного TSV
(либо готового FASTA) до эволюционных деревьев, уверенных клад и таблиц мутаций.
Каждый этап — самостоятельный скрипт, соединённые в единую цепочку через файлы
на диске (tsv → fasta → fasta → nexus/newick → json/html/tsv).

## Конвейер

| # | Этап | Папка в `scripts/` | Вход | Выход |
|---|---|---|---|---|
| 1 | Фильтрация сиквенсов | `01_filter_sequences/` | `data/<batch>/BCR_data.tsv` (или FASTA) | `<out>/BCR_data_filtered.tsv` |
| 2 | Группировка по гермлайнам | `02_group_by_germlines/` | `<out>/BCR_data_filtered.tsv` | `<out>/grouped_by_germlines/{v,d,j,vj}/*.fasta` |
| 3 | Множественное выравнивание (MAFFT) | `03_multiple_alignment/` | `<out>/grouped_by_germlines/vj/*.fasta` | `aligned_sequences/*_aligned.fasta`, `manifest.tsv` |
| 4 | ML-деревья (IQ-TREE) | `04a_build_trees_iqtree/` | `aligned_sequences/` | `trees/<group>/<group>.treefile` |
| 5 | Байесовские деревья (MrBayes) | `04b_build_trees_mrbayes/` | `aligned_sequences/*.fasta` | `mrbayes/<group>.nex.con.tre` |
| 6 | Уверенные клады | `05_clade_search/` | `mrbayes/*.nex.con.tre` + `trees/*/*.treefile` | `groups/report.json` |
| 7 | Анализ мутаций | `06_analyze_mutations/` | fasta по кладам | `<out>/<clade>/mutations.tsv`, `mutations_summary.tsv` |

Вспомогательные скрипты (не являются шагами конвейера):

| Этап | Папка в `scripts/` | Вход | Выход |
|---|---|---|---|
| Визуализация деревьев | `visualize_trees/` | `trees/` или `mrbayes/` | `<out>/<group>/<group>.html` |
| Аналитика клад | `tree_analytics/` | `groups/report.json` | `tree_analytics.png` |

Все скрипты принимают пути явными аргументами (`-i`/`-o`/`-r`).

## Требования

- [Miniconda](https://docs.anaconda.com/miniconda/) или Anaconda
- Python 3.11 (ставится через conda)

Все зависимости описаны в `environment.yml` и устанавливаются одной командой:

```bash
conda env create -f environment.yml
conda activate biocad_bcr_pipeline_environment
```

В окружение входят: `pandas`, `biopython`, `mafft`, `iqtree`, `mrbayes`, `igblast`,
`blast`, `toytree`, `toyplot`, `matplotlib`, `seaborn`.

## Быстрый старт

Поместите сырые данные в `data/<batch>/BCR_data.tsv` (AIRR-формат).
Гермлайновые справочники IMGT — в `data/references/`.

### Весь конвейер одной командой

```bash
python3 scripts/run_pipeline.py -k <batch>
```

Пример:

```bash
python3 scripts/run_pipeline.py -k BCR
```

Можно указать явные пути и кастомный выход:

```bash
python3 scripts/run_pipeline.py -i data/BCR/BCR_data.tsv -o results/report_$(date +%d%m%Y)
```

Частичный запуск:

```bash
# только шаги 3-6 (пропустить 1-2)
python3 scripts/run_pipeline.py -k BCR --skip filter group

# только фильтрация + группировка
python3 scripts/run_pipeline.py -k BCR --only filter group

# IQ-TREE + MrBayes параллельно (экспериментально)
python3 scripts/run_pipeline.py -k BCR --parallel-trees
```

### Поэтапно

```bash
export BATCH=BCR
export OUT=results/report_$(date +%d%m%Y_%H%M%S)

# 1: фильтрация
python3 scripts/01_filter_sequences/filter_sequences.py \
  -i data/$BATCH/BCR_data.tsv \
  -o $OUT/BCR_data_filtered.tsv \
  -r data/references

# 2: группировка по гермлайнам
python3 scripts/02_group_by_germlines/group_by_germlines.py \
  -i $OUT/BCR_data_filtered.tsv \
  -o $OUT/grouped_by_germlines \
  -r data/references

# 3: выравнивание (MAFFT)
python3 scripts/03_multiple_alignment/multiple_alignment.py \
  -i $OUT/grouped_by_germlines/vj \
  -o aligned_sequences

# 4: ML-деревья (IQ-TREE)
python3 scripts/04a_build_trees_iqtree/build_trees_iqtree.py \
  -i aligned_sequences -o trees

# 5: байесовские деревья (MrBayes)
python3 scripts/04b_build_trees_mrbayes/build_trees_mrbayes.py \
  aligned_sequences --out mrbayes

# 6: уверенные клады
python3 scripts/05_clade_search/clade_search.py \
  --mrbayes-dir mrbayes --iqtree-dir trees --out groups/report.json

# 7: анализ мутаций
python3 scripts/06_analyze_mutations/run_mutations.py \
  -i fasta_from_clades -o mutation_tables -r data/references

# визуализация деревьев (отдельно)
python3 scripts/visualize_trees/visualize_trees.py -i trees -o trees_visualization
python3 scripts/visualize_trees/visualize_trees.py -i mrbayes -o mrbayes_visualization

# аналитика клад (отдельно)
python3 scripts/tree_analytics/tree_analytics.py -i groups/report.json
```

## Структура репозитория

```
data/
  <batch>/BCR_data.tsv        — сырой вход (AIRR-подобный TSV)
  references/                 — гермлайновые справочники IMGT (IGHV.fasta, …)
results/
  example_report/             — пример результатов полного прогона
scripts/
  <этап>/                     — код этапа + README
  run_pipeline.py             — master-раннер всего конвейера
tests/
  test_<этап>/                — pytest для большинства этапов
  test_fixtures.py            — консольные тесты для шагов 5-6
  test_pipeline.py            — E2E на реальных данных
  fixtures/                   — статические тестовые данные
  visualize_tree.py           — опциональная визуализация результатов тестов
misc/                         — зарезервировано
```

## Тестирование

```bash
conda activate biocad_bcr_pipeline_environment

# pytest (шаги 1-4, 7)
pytest tests/ -v

# или выборочно
pytest tests/test_build_trees_iqtree/ -v
pytest tests/test_analyze_mutations/ -v
pytest tests/test_visualize_trees/ -v

# консольные тесты (шаги 5-6, без pytest)
python tests/test_fixtures.py
python tests/test_pipeline.py 1
```

Тестовые отчёты и разбор найденных багов — в [tests/TEST_REPORT.md](tests/TEST_REPORT.md).

## Известные пробелы

- **`biocode/` не опубликован.** Шаги 5 и 6 импортируют пакет `biocode`
  (вендорен из `BIOCAD.bigchallenges@main`) — без него падают на импорте.
- **Разрыв между шагами 6 и 7.** `groups/report.json` (шаг 6) хранит клады как
  списки id в JSON, а шаг 7 ожидает отдельные fasta-файлы на каждую кладу.
  Скрипт-прослойка (`report.json` → `fasta_from_clades/*.fasta`) пока не написан.
- `.gitignore` и `gitignore.txt` дублируются; активный — `.gitignore`.

## Участники

| Участник | Этапы |
|---|---|
| Ксюша | фильтрация, группировка по гермлайнам |
| Алина | множественное выравнивание (MSA) |
| Денис | ML-деревья (IQ-TREE), визуализация, анализ мутаций |
| Никита | NEXUS/MrBayes, уверенные клады |
