# Tests

## Структура

```
tests/
  ├── README.md                       # этот файл
  ├── TESTING.md                      # как запустить
  ├── TEST_REPORT.md                  # разбор находок и исправленных багов
  ├── fixtures/                       # СТАТИЧЕСКИЕ данные с заранее известными мутациями
  │   ├── two_clear_pairs_aligned.fasta      + .expected.json
  │   ├── no_shared_mutations_aligned.fasta  + .expected.json
  │   └── mixed_signal_aligned.fasta         + .expected.json
  ├── test_analyze_mutations_unit.py  # юнит на разметку мутаций (без igblast/mb)
  ├── test_fixtures.py                # E2E на фикстурах (нужен mb)
  ├── test_pipeline.py                # E2E-прогон mrbayes + clade_search на фикстурах (нужен mb)
  ├── test_*.py                       # юнит-тесты остальных шагов
  └── visualize_tree.py               # HTML/SVG-визуализация — ОТДЕЛЬНЫЙ файл, опционально
```

## Принцип: тесты и визуализация разделены

`test_fixtures.py` и `test_pipeline.py` не импортируют `visualize_tree.py`. Они
работают в консоли: запускают `04b_build_trees_mrbayes/build_trees_mrbayes.py`
и `05_clade_search/clade_search.py` как subprocess, читают
`.nex.con.tre`/`report.json`, печатают pass/fail.

## Fixtures — статические файлы

`tests/fixtures/*_aligned.fasta` — тот же формат, что даёт
`03_multiple_alignment` (`mafft --auto`, суффикс `_aligned.fasta`).
Последовательности из повторяющегося паттерна `ACGT` — любая мутация видна
глазами (напр. `...ACGTACTTACGT...` вместо `...ACGTACGTACGT...`). Рядом
`*.expected.json`: корневая последовательность, позиции/типы мутаций каждого
листа, ожидаемые уверенные клады с обоснованием. Без `random` в рантайме.

## Запуск

```bash
# быстрый регресс на разметку мутаций (чистый Python, без внешних тулов)
python tests/test_analyze_mutations_unit.py

# E2E на статических фикстурах (нужен mb)
python tests/test_fixtures.py

# E2E-прогон mrbayes + clade_search на фикстурах (нужен mb)
python tests/test_pipeline.py 3

# визуализация — отдельно, опционально
python tests/visualize_tree.py <group_key> mrbayes/ --report groups/report.json
```

## Требования

- Python 3.11, `biopython` (для скриптов, вызываемых тестами, и `visualize_tree.py`)
- бинарь `mb` (MrBayes): `conda install -c bioconda mrbayes` — для test_fixtures/test_pipeline
- для `test_analyze_mutations/` (smoke run_mutations.py) — `igblast` + BLAST-базы

Подробности находок — [TEST_REPORT.md](TEST_REPORT.md).
