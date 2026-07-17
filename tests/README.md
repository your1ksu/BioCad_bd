# Tests — открытость и достоверность pipeline Никиты

## Структура

```
tests/
  ├── README.md               # этот файл
  ├── TESTING.md               # как запустить
  ├── TEST_REPORT.md           # разбор находок, включая найденный и исправленный баг
  ├── fixtures/                 # СТАТИЧЕСКИЕ данные с заранее известными мутациями
  │   ├── two_clear_pairs_aligned.fasta
  │   ├── two_clear_pairs.expected.json
  │   ├── no_shared_mutations_aligned.fasta
  │   ├── no_shared_mutations.expected.json
  │   ├── mixed_signal_aligned.fasta
  │   └── mixed_signal.expected.json
  ├── test_fixtures.py          # тест на статических fixtures (консоль, без HTML)
  ├── test_pipeline.py          # E2E тест на реальных данных (консоль, без HTML)
  └── visualize_tree.py         # HTML/SVG-визуализация — ОТДЕЛЬНЫЙ файл, опционально
```

## Ключевой принцип: тесты и визуализация разделены

`test_fixtures.py` и `test_pipeline.py` не импортируют `visualize_tree.py` вообще.
Они работают в чистой консоли: запускают `04b_build_trees_mrbayes/build_trees_mrbayes.py` и
`05_clade_search/clade_search.py` как subprocess, читают получившиеся
`.nex.con.tre`/`report.json`, печатают pass/fail. `visualize_tree.py` — отдельный
скрипт, который по желанию рендерит HTML-филограмму из уже посчитанных
результатов, и без него тестирование работает нормально.

## Fixtures — не runtime-генератор, а статические файлы

`tests/fixtures/*_aligned.fasta` — точно формат, который производит
`Alina/03_multiple_alignment/multiple_alignment.py` (mafft --auto) и который читает `04b_build_trees_mrbayes/build_trees_mrbayes.py`.
Последовательности построены из повторяющегося паттерна `ACGT` — любая
мутация видна невооружённым глазом (напр. `...ACGTACTTACGT...` вместо
`...ACGTACGTACGT...`). Рядом — `*.expected.json`: корневая последовательность,
точные позиции и типы мутаций каждого листа, ожидаемые уверенные клады с
обоснованием. Никакого `random` в рантайме теста — файлы построены один раз
детерминированно и просто читаются тестом.

## Запуск

```bash
# Из корня проекта
python tests/test_fixtures.py                 # статические fixtures (~65 сек, 3 кейса)
python tests/test_pipeline.py 1                # E2E на реальных данных (~35 сек, 1 группа)

# Визуализация — отдельно, опционально, после того как выше уже отработало:
python tests/visualize_tree.py IGHV3-23_01_IGHJ3_01 04b_build_trees_mrbayes/ \
    --report 05_clade_search/report.json --out 04b_build_trees_mrbayes/IGHV3-23_01_IGHJ3_01.tree.html
```

## Требования

- Python 3.9+
- `biopython` (pip install biopython) — нужен и test_fixtures.py/test_pipeline.py
  (транзитивно, через вызываемые 04b_build_trees_mrbayes/05_clade_search скрипты), и visualize_tree.py
- Бинарь `mb` (MrBayes): conda install -c bioconda mrbayes
- Для test_pipeline.py — клон BIOCAD.bigchallenges на ветке main (реальные
  данные Ксюши/Алины)

## Связь с реальным pipeline (см. TEST_REPORT.md)

Формат fixtures подтверждён по факту прочтения реальных скриптов других
участников (ветки your1ksu/BioCad_bd): `ksu_branch/scripts/group_by_germlines`
→ `Alina/scripts/03_multiple_alignment/multiple_alignment.py` (mafft --auto, суффикс `_aligned.fasta`) →
наш `04b_build_trees_mrbayes/build_trees_mrbayes.py`. Параллельно — `Denis/scripts/04a_build_trees_iqtree`
(вход тот же `aligned_sequences/`, выход `trees/<группа>/<группа>.treefile`),
который `05_clade_search/clade_search.py` тоже умеет читать (`--iqtree-dir`).
