# BioCad_bd — филогенетический конвейер по антителам (BCR repertoire)

Конвейер обработки B-клеточных рецепторов (BCR): от сырых прочтений до
эволюционных деревьев, уверенных клад и таблиц мутаций.

## Быстрый старт

```bash
# Весь пайплайн одним запуском:
python3 scripts/run_pipeline.py -k batch1

# С выбором выравнивателя и стратегии группировки:
python3 scripts/run_pipeline.py -k batch1 --aligner macse --grouping-strategy gene
```

`run_pipeline.py` — мастер-оркестратор. Он сам находит/создаёт conda-окружение
`biocad_bcr_pipeline_environment`, прогоняет все шаги и сохраняет конфиг.

## Конвейер

| # | Этап | Папка в `scripts/` | Вход | Выход |
|---|---|---|---|---|
| 1 | Фильтрация сиквенсов | `01_filter_sequences/` | `data/<key>/BCR_data.tsv` | `BCR_data_filtered.tsv` |
| 2 | Группировка по гермлайнам | `02_group_by_germlines/` | `BCR_data_filtered.tsv` | `grouped_by_germlines/vj/*.fasta` |
| 2b | Фильтрация групп по размеру | `filter_by_symbol_count/` | `vj/` | `vj_filtered/` |
| 2c | Верификация по аминокислотам | `verify_by_amino/` | `vj_filtered/` | `verify_by_amino/` (только перед MAFFT) |
| 3 | Множественное выравнивание (MSA) | `03_multiple_alignment/` | `vj_filtered/` или `verify_by_amino/` | `aligned_sequences/` |
| 4a | ML-деревья (IQ-TREE) | `04a_build_trees_iqtree/` | `aligned_sequences/` | `trees/<группа>/` |
| 4b | Байесовские деревья (MrBayes) | `04b_build_trees_mrbayes/` | `aligned_sequences/` | `mrbayes/` |
| 5 | Визуализация деревьев | `visualize_trees/` | `trees/` или `mrbayes/` | `trees_visualization/` |
| 6 | Уверенные клады | `05_clade_search/` | `trees/` + `mrbayes/` | `groups/report.json` + `groups/clades/` |
| 7 | Анализ мутаций | `06_analyze_mutations/` | fasta по кладам | `mutation_tables/` |

## Выравниватели

| Выравниватель | Флаг | Когда выбирать |
|---|---|---|
| MAFFT (по умолч.) | `--aligner mafft` | Быстрое выравнивание; требует `verify_by_amino` перед собой |
| MACSE | `--aligner macse` | Корректная работа с кодонами, frameshift, стоп-кодонами; verify не нужен (MACSE сам транслирует) |

## Стратегии группировки

| Стратегия | Пример ключа группы | Описание |
|---|---|---|
| `gene` (по умолч.) | `IGHV1-2_IGHJ4` | V-ген + J-ген, без аллелей |
| `allele` | `IGHV1-2\*01_IGHJ4\*01` | V-ген + J-ген с аллелями (как было изначально) |
| `v_only` | `IGHV1-2` | Только V-ген |

## Конфигурация

Все параметры пайплайна (пути, пороги, модели) хранятся в
`scripts/shared/config.py` (класс `PipelineConfig`). Можно передать JSON-конфиг:

```bash
python3 scripts/run_pipeline.py -k batch1 --config my_config.json
```

Ключи CLI переопределяют значения из JSON.

## Требования

- Python 3.11+, [Miniconda](https://docs.anaconda.com/miniconda/)
- Окружение `biocad_bcr_pipeline_environment` (создаётся автоматически из
  `environment.yml`): `mafft`, `macse`, `iqtree`, `mrbayes`, `pandas`,
  `biopython`, `toytree`, `toyplot`, `matplotlib`, `seaborn`

## Тестирование

```bash
python3 -m pytest tests/ -v
```

## Структура репозитория

```
data/<key>/BCR_data.tsv          — сырой вход (AIRR-подобный TSV)
data/references/                 — germline-справочники IMGT
results/<report>/                — результаты запуска (лог, конфиг, все шаги)
scripts/<этап>/                  — код каждого этапа + README
scripts/shared/                  — общие модули (config.py, utils.py)
tests/                           — тесты pytest
environment.yml                  — conda-окружение
```
