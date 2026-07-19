# Тестирование Pipeline

Два типа тестов, оба консольные (без HTML/GUI-зависимостей):

## 1. Тесты на статических fixtures с известным ответом

```bash
python tests/test_fixtures.py
```

Читает `tests/fixtures/*_aligned.fasta` (статические файлы с задокументированными
мутациями — см. `*.expected.json` рядом), прогоняет через реальный
`mrbayes/run_mrbayes.py` → `clades/confident_clades_report.py`, сверяет найденные
уверенные клады с ожидаемыми.

```
Пройдено: 3/3
  ✓ mixed_signal
  ✓ no_shared_mutations
  ✓ two_clear_pairs
```

## 2. E2E тест на реальных данных

```bash
python tests/test_pipeline.py [N_GROUPS]
```

Загружает N реальных групп из `BIOCAD.bigchallenges/anotherpipeline/`
(настоящие данные Ксюши/Алины), прогоняет через pipeline, печатает пути к
результатам. Визуализацию не генерирует — для неё см. ниже.

```bash
python tests/test_pipeline.py 1      # одна группа (~35 сек)
python tests/test_pipeline.py 3      # три группы (~100 сек)
```

## 3. Визуализация (отдельно, опционально)

```bash
python tests/visualize_tree.py <group_key> mrbayes/ \
    --report clades/report.json --out mrbayes/<group_key>.tree.html
```

Не требуется для тестирования — намеренно вынесена в отдельный файл, который
`test_fixtures.py`/`test_pipeline.py` не импортируют.

## Требования

- Python 3.9+, `biopython` (pip install biopython)
- Бинарь `mb` (MrBayes): conda install -c bioconda mrbayes
- Для test_pipeline.py — клон BIOCAD.bigchallenges на ветке main

## Интерпретация результатов

Смотрите [TEST_REPORT.md](TEST_REPORT.md) — там разобрана находка: реальный
баг в извлечении уверенных клад из MrBayes-консенсуса (часть клады терялась
при неукоренённом дереве с базовой политомией) и то, как он был исправлен и
покрыт тестом.
