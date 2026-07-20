# Тестирование Pipeline

## 1. Юнит на разметку мутаций (быстро, без внешних тулов)

```bash
python tests/test_analyze_mutations_unit.py
```

Проверяет СОДЕРЖИМОЕ таблиц `analyze_mutations` (синонимная → `is_silent=yes`,
несинонимная → `no`, summary считает несинонимные). Ловит регрессию, при которой
`is_silent` записывался строкой и summary выходил нулевым. Чистый Python.

## 2. Тесты на статических fixtures с известным ответом

```bash
python tests/test_fixtures.py
```

Читает `tests/fixtures/*_aligned.fasta` (см. `*.expected.json` рядом), прогоняет
через `04b_build_trees_mrbayes/build_trees_mrbayes.py` →
`05_clade_search/clade_search.py`, сверяет найденные уверенные клады с ожидаемыми.

## 3. E2E-прогон на фикстурах

```bash
python tests/test_pipeline.py [N_GROUPS]     # напр. 3
```

Берёт N групп из `tests/fixtures/`, прогоняет mrbayes + clade_search, печатает
пути к результатам (`.nex.con.tre`, `report.json`). Визуализацию не генерирует.

## 4. Визуализация (отдельно, опционально)

```bash
python tests/visualize_tree.py <group_key> mrbayes/ --report groups/report.json
```

Не требуется для тестирования — намеренно вынесена в отдельный файл, который
`test_fixtures.py`/`test_pipeline.py` не импортируют.

## Требования

- Python 3.11, `biopython`
- бинарь `mb` (MrBayes): `conda install -c bioconda mrbayes` — для п.2–3

## Интерпретация

См. [TEST_REPORT.md](TEST_REPORT.md) — там разобраны находки: регрессия `is_silent`
в разметке мутаций и потеря симметричной клады из MrBayes-консенсуса (обе
исправлены и покрыты тестами).
