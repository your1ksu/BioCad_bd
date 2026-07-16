# Test

Полный набор тестов для проверки корректности и воспроизводимости pipelin (nexus/MrBayes/уверенные клады).

## Структура

- **test_pipeline.py** — E2E тесты на реальных данных из BIOCAD.bigchallenges
- **test_synthetic.py** — Синтетические тесты с контролируемой филогенетической структурой
- **TESTING.md** — Краткая инструкция по запуску
- **TEST_REPORT.md** — Полный анализ результатов и выводы о достоверности

## Запуск

```bash
# Из корневой папки проекта (там где mrbayes/, groups/, etc)
cd /path/to/BioCad_bd  # или /Users/nikitasyzdykov/Desktop/biocad/biocad

# E2E тест на реальных данных (одна группа, ~32 сек)
python tests/test_pipeline.py 1

# Синтетические тесты с известным ответом (все 4 теста, ~70 сек)
VENV=/path/to/venv/bin
python tests/test_synthetic.py "$VENV/python"
```

## Требования

- Python 3.9+
- biopython: `pip install biopython`
- MrBayes binary `mb` (опционально для фазы 2): `conda install -c bioconda mrbayes`
- Клон BIOCAD.bigchallenges на ветке main (для test_pipeline.py)

## Что тестируется

| Аспект | Тест | Статус |
|--------|------|--------|
| Функциональность | E2E на реальных данных | ✓ |
| Воспроизводимость | Повторные запуски | ✓ |
| Граничные случаи | Синтетические политомии | ✓ |
| Восстановление топологии | Синтетические деревья | ~ |
| Апостериорная поддержка | Консенсусное дерево | ✓ |

Подробнее: см. TEST_REPORT.md
