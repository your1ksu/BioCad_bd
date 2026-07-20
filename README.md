# BioCad_bd — филогенетический конвейер по антителам (BCR repertoire)

Конвейер обработки B-клеточных рецепторов (BCR): от репертуара антител (AIRR TSV
или FASTA) до филогенетических деревьев клональных групп, набора надёжных клад и
таблиц мутаций с разметкой по регионам FR/CDR (IMGT). Все шаги оркеструются одним
скриптом `scripts/run_pipeline.py`.

## Быстрый старт

```bash
# по ключу (data/<key>/BCR_data.tsv -> results/report_<timestamp>/)
python scripts/run_pipeline.py -k BCR

# по явным путям входа/выхода
python scripts/run_pipeline.py -i data/BCR/BCR_data.tsv -o results/output
```

При первом запуске оркестратор проверяет/создаёт conda-окружение
`biocad_bcr_pipeline_environment` из `environment.yml` и выполняет все шаги внутри него.

## Шаги пайплайна

`ALL_STEPS = [filter, group, filter_groups, msa, iqtree, mrbayes, viz, clades, mutations]`

| Шаг | Скрипт | Вход → выход |
|---|---|---|
| `filter` | `scripts/01_filter_sequences/filter_sequences.py` | AIRR TSV/FASTA → `BCR_data_filtered.tsv` |
| `group` | `scripts/02_group_by_germlines/group_by_germlines.py` | filtered TSV → `grouped_by_germlines/vj/<V_J>.fasta` |
| `filter_groups` | `scripts/filter_by_symbol_count/filter_by_symbol_count.py` | `vj/` → `vj_filtered/` (симлинки групп нужного размера) |
| `msa` | `scripts/03_multiple_alignment/multiple_alignment.py` | группы → `aligned_sequences/<V_J>_aligned.fasta` (MAFFT) |
| `iqtree` | `scripts/04a_build_trees_iqtree/build_trees_iqtree.py` | выравнивания → `trees/<V_J>/<V_J>.treefile` (ML, UFBoot+SH-aLRT) |
| `mrbayes` | `scripts/04b_build_trees_mrbayes/build_trees_mrbayes.py` | выравнивания → `mrbayes/<V_J>.nex(.con.tre)` (Bayes) |
| `viz` | `scripts/visualize_trees/visualize_trees.py` | деревья → `trees_visualization/…/*.html` (toytree) |
| `clades` | `scripts/05_clade_search/clade_search.py` | деревья + выравнивания → `groups/report.json` + `groups/clades/*.fa` |
| `mutations` | `scripts/06_analyze_mutations/run_mutations.py` | FASTA клад → `mutation_tables/<клада>/mutations.tsv(+summary)` (IgBLAST) |

Управление составом:

```bash
python scripts/run_pipeline.py -k BCR --only filter group filter_groups msa
python scripts/run_pipeline.py -k BCR --skip mrbayes --parallel-trees
```

Зависимости шагов: `mutations` подтягивает `clades`; `filter_groups` — `group`.
Переопределения из CLI: `--grouping-strategy`, `--min-group-size`, `--max-group-size`,
`--iqtree-model`, `--parallel-trees`, `--gpu-mb-bin`.

### CPU / GPU для MrBayes

По умолчанию MrBayes считается на CPU (параллельно по группам, каждая в своей
изолированной подпапке, со stoprule по сходимости). Для ускорения крупных групп
на GPU (A100, BEAGLE-CUDA) укажите путь к GPU-бинарю `mb`:

```bash
python scripts/run_pipeline.py -k BCR --gpu-mb-bin /path/to/mrbayes-gpu/bin/mb
```

Тогда группы с числом таксонов ≥ `gpu_min_taxa` (по умолчанию 60) идут на GPU
последовательно, мелкие — параллельно на CPU. Без `--gpu-mb-bin` — всё на CPU.
GPU-бинарь MrBayes с BEAGLE-CUDA нужно собрать отдельно (в conda его нет).

## Структура репозитория

```
scripts/
  run_pipeline.py              оркестратор всех шагов
  shared/config.py             PipelineConfig (все параметры, сериализуется в config.json)
  shared/utils.py              общие утилиты (FASTA, кодоны, ядра)
  01_filter_sequences/ … 06_analyze_mutations/   шаги пайплайна (см. таблицу)
  filter_by_symbol_count/      отсев групп по размеру (шаг 2b)
  visualize_trees/             HTML-визуализация деревьев
  mutation_stats/              dN/dS и распределение мутаций по регионам (отдельно)
  SHM/                         поиск SHM hotspot/coldspot мотивов (отдельно)
  tree_analytics/              сводная аналитика по report.json (отдельно)
  verify_by_amino/             трансляция в аминокислоты + сверка с IMGT-справочником
  paths.py                     единые пути data/<key>/ ↔ results/<key>/ (для verify_by_amino)
  synthetic_data/              генератор синтетики + валидация разметки FR/CDR
  MACSE/                       кодон-осознанное выравнивание (альтернатива MAFFT)
data/
  BCR/BCR_data.tsv             входной репертуар (AIRR)
  references/                  germline IMGT (IGHV/…), базы IgBLAST (all_V/D/J), human_gl.aux
  example.fasta, example.tsv   маленькие входы для быстрых прогонов
tests/                         pytest + консольные E2E (см. tests/README.md)
results/                       результаты прогонов (в .gitignore)
```

## Выходная директория одного прогона

```
results/report_<timestamp>/
├── config.json                       снимок конфигурации
├── logs/                             логи по шагам
├── BCR_data_filtered.tsv
├── grouped_by_germlines/vj[/_filtered]
├── aligned_sequences/*_aligned.fasta + manifest.tsv
├── trees/<V_J>/<V_J>.treefile        (IQ-TREE)
├── mrbayes/<V_J>.nex.con.tre         (MrBayes)
├── trees_visualization/…/*.html
├── groups/report.json + groups/clades/*.fa
└── mutation_tables/<клада>/mutations.tsv (+ mutations_summary.tsv)
```

## Требования

- Python 3.11, conda/miniconda.
- Окружение `biocad_bcr_pipeline_environment` (`environment.yml`): `pandas`, `biopython`,
  `mafft`, `iqtree`, `mrbayes`, `igblast`, `blast`, `toytree`, `toyplot`, `matplotlib`, `seaborn`.

## Тесты

```bash
pytest tests/                                   # юнит-тесты шагов (часть требует внешних тулов)
python tests/test_analyze_mutations_unit.py     # быстрый регресс на разметку мутаций (без igblast)
python tests/test_fixtures.py                   # E2E на фикстурах (нужен MrBayes)
```

Подробности — `tests/README.md`, `tests/TEST_REPORT.md`.

## Известные ограничения

- В `report.json` поля `defining_mutations`, `defining_cdr`, `isotypes`, `days` пока не
  вычисляются (заглушки); `confident_both_models` заполняется (пересечение клад ML и Bayes).
- `mutation_stats.py` (dN/dS), `SHM/`, `tree_analytics.py`, `synthetic_data/validate_regions.py`
  запускаются отдельно от оркестратора, не входят в `run_pipeline.py`.
- `scripts/MACSE` требует ручной настройки путей и в оркестратор не встроен.

Каждый шаг снабжён своим README рядом со скриптом.
