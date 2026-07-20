# MACSE MSA

Множественное выравнивание FASTA-файлов с помощью MACSE.

Скрипт ищет все FASTA-файлы во входной папке, запускает MACSE для каждого
файла и сохраняет:

- `*_aligned.fasta` — нуклеотидное выравнивание;
- `*_aligned_aa.fasta` — аминокислотное выравнивание;
- `manifest.tsv` — таблицу соответствия.

## Требования

- Python 3.11+
- MACSE в PATH (или путь через `--macse`)

```bash
conda install -c bioconda macse
```

## Запуск

```bash
python3 MACSEtry.py -i /path/to/input_fastas -o /path/to/output_dir

# С явным путём к MACSE:
python3 MACSEtry.py -i /path/to/input_fastas -o /path/to/output_dir --macse /opt/macse/bin/macse
```

## Аргументы

| Аргумент | Описание |
|---|---|
| `-i`, `--input` | Папка с входными FASTA-файлами |
| `-o`, `--output` | Папка для результатов |
| `--macse` | Путь к MACSE (по умолч. поиск в PATH) |
| `--threads` | Потоков на процесс (по умолч. 1) |
| `--workers` | Параллельных процессов (по умолч. число ядер CPU) |

## Использование в пайплайне

Через `run_pipeline.py`:

```bash
python3 scripts/run_pipeline.py -k batch1 --aligner macse
```

Или через `multiple_alignment.py`:

```bash
python3 scripts/03_multiple_alignment/multiple_alignment.py -i input -o output --aligner macse
```

## Важно

MACSE предназначен для кодирующих нуклеотидных последовательностей. Если
входные последовательности не являются coding DNA, имеют проблемную рамку
считывания или много стоп-кодонов, MACSE может завершиться с ошибкой.

## Тесты

```bash
python3 -m pytest tests/tests_MACSE/ -v
```

Тесты не запускают настоящий MACSE: внешний вызов подменяется фейком.