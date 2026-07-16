# Тестирование Pipeline

Два типа тестов для проверки открытости и достоверности:

## 1. E2E тесты на реальных данных

Загружает реальные группы из `BIOCAD.bigchallenges/anotherpipeline/`, прогоняет через pipeline, генерирует визуализации.

```bash
python test_pipeline.py [N_GROUPS]
```

**Примеры:**
```bash
python test_pipeline.py 1      # одна группа (~32 сек)
python test_pipeline.py 3      # три группы (~90 сек)
```

**Выход:**
- `aligned_sequences/` — копии исходных FASTA
- `mrbayes/` — NEXUS, консенсус, лог MrBayes, `.tree.html` визуализации
- `groups/report.json` — уверенные клады всех групп

---

## 2. Синтетические тесты с известным ответом

Создаёт филогенетические деревья с контролируемой структурой, эволюционирует по ним последовательности, запускает pipeline, проверяет восстановление топологии.

```bash
python test_synthetic.py /path/to/venv/bin/python
```

**Пример (с тестовым окружением):**
```bash
VENV=/private/tmp/claude-501/.../scratchpad/venv
python test_synthetic.py "$VENV/bin/python"
```

Запускает 4 теста:
- ✓ **star_tree_5taxa** — политомия (0 клад ожидается, 0 найдено)
- ⚠️ **well_separated_pairs** — две пары с хорошим разделением
- ⚠️ **deep_nested_tree** — глубокая вложенная структура
- ⚠️ **uneven_branches** — неравномерные длины ветвей

**Выход:**
```
Прошло: 1/4
  ✓ star_tree_5taxa
  ❌ well_separated_pairs
  ❌ deep_nested_tree
  ❌ uneven_branches
```

Результаты интерпретируются в `TEST_REPORT.md`.

---

## Требования

- Python 3.9+
- `biopython` (pip install biopython)
- Бинарь `mb` (conda install -c bioconda mrbayes, опционально для фазы 2 MrBayes)
- Клон репозитория `BIOCAD.bigchallenges` с ветки `origin/main` (для test_pipeline.py)

---

## Интерпретация результатов

Смотрите [TEST_REPORT.md](TEST_REPORT.md) для подробного анализа:
- Почему некоторые синтетические тесты "не проходят" (это нормально)
- Что это говорит о достоверности pipeline
- Какие ограничения есть и когда их ожидать
