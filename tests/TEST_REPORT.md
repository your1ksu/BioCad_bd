# Отчёт о тестировании pipeline

## 1. Регрессия разметки мутаций: `is_silent` (найдена и исправлена)

### Симптом
`mutations_summary.tsv` выходил весь в нулях, а колонка `is_silent` в `mutations.tsv` —
всегда `yes`, хотя в данных есть несинонимные замены. Как следствие, `mutation_stats.py`
не мог посчитать dN/dS (все замены выглядели синонимными).

### Причина
В `06_analyze_mutations/analyze_mutations.py` `compute_mutations` клал `is_silent`
**строкой** `'yes'/'no'`, а потребители ждут bool:
- `write_mutations`: `'yes' if m['is_silent'] else 'no'` — обе строки истинны → всегда `yes`;
- `write_summary`: `if not m['is_silent']` — всегда False → счётчики не растут.

Старый тест `test_analyze_mutations` проверял только НАЛИЧИЕ файлов и этот баг не видел.

### Исправление
`compute_mutations` возвращает `is_silent` булевым (`'is_silent': is_silent`).
Добавлен юнит-тест `tests/test_analyze_mutations_unit.py` — без igblast/mb, проверяет
СОДЕРЖИМОЕ таблиц (синонимная → `yes`, несинонимная → `no`, summary считает несинонимные).
На «строковом» коде тест падает, на исправленном — проходит.

## 2. Потеря симметричной клады у MrBayes (найдена и исправлена ранее)

### Симптом
Fixture `two_clear_pairs` (4 таксона, 2 симметричные пары по 6 общих мутаций) — ожидались
обе клады `{seq1,seq2}` и `{seq3,seq4}` с posterior≈1.0, находилась только одна:

```
tree con_50_majrule = (1:0.0072,2:0.0074,(3:0.0073,4:0.0073)1.000:0.0627);
```

### Причина
MrBayes пишет консенсус НЕУКОРЕНЁННОГО дерева как rooted-newick с базовой политомией:
`seq1`,`seq2` висят на корне отдельными ветвями, а `(seq3,seq4)` обёрнуты в узел с posterior.
Для 4 таксонов это единственное внутреннее ребро — `{seq1,seq2}` обладает той же posterior,
что и `{seq3,seq4}` (две стороны одного разбиения), но наивный обход `get_nonterminals()`
видит только обёрнутую сторону.

### Исправление
В `05_clade_search/clade_search.py`:
- `_root_complement_supports` — достраивает комплементарную кладу из голых листьев корня;
- `_drop_complement_duplicates` — убирает бо́льшую сторону разбиения при РАЗНОМ размере,
  оставляя обе при РАВНОМ (как 2+2).

## 3. Фикстуры — статические данные с известными мутациями

`tests/fixtures/*_aligned.fasta` — статические файлы: последовательности из повторяющегося
паттерна `ACGT`, любая мутация видна глазами. Формат совпадает с выходом
`03_multiple_alignment` (`mafft --auto`, суффикс `_aligned.fasta`). Рядом `*.expected.json` —
ground truth: корень, позиции/ref/alt каждого листа, ожидаемые надёжные клады с обоснованием.

| Fixture | Таксонов | Структура | Ожидаемые клады |
|---|---|---|---|
| `two_clear_pairs` | 4 | 2 симметричные пары (по 6 общих мутаций) | `{seq1,seq2}`, `{seq3,seq4}` |
| `no_shared_mutations` | 5 | только приватные мутации (звезда) | нет (отрицательный контроль) |
| `mixed_signal` | 6 | 2 пары + 2 одиночки | `{seqA1,seqA2}`, `{seqB1,seqB2}` |

## 4. Файлы тестов

- `tests/test_analyze_mutations_unit.py` — юнит на разметку мутаций (§1). Чистый Python.
- `tests/test_fixtures.py` — E2E на фикстурах: `04b_build_trees_mrbayes/build_trees_mrbayes.py`
  → `05_clade_search/clade_search.py`, сверка с `*.expected.json`. Нужен `mb`.
- `tests/test_pipeline.py` — E2E-прогон mrbayes + clade_search на фикстурах (консоль). Нужен `mb`.
- `tests/test_analyze_mutations/` — smoke-тест `run_mutations.py` (наличие файлов). Нужен igblast.
- `tests/test_*` — юнит-тесты остальных шагов (filter, group, filter_by_symbol_count, MSA,
  iqtree, visualize, translate, SHM, tree_analytics).

> Историческая правка: `test_fixtures.py`/`test_pipeline.py` раньше ссылались на скрипты
> `mrbayes/run_mrbayes.py`, `clades/confident_clades_report.py` и внешний репозиторий
> `BIOCAD.bigchallenges` из старой раскладки — перенаправлены на актуальные нумерованные
> скрипты и локальные фикстуры.
