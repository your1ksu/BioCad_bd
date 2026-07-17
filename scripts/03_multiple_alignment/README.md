# MSA: выравнивание FASTA-файлов через MAFFT

`multiple_alignment.py` выполняет множественное выравнивание FASTA-файлов с помощью
MAFFT.

## Запуск

```bash
python3 multiple_alignment.py -i /path/to/input_fastas -o /path/to/output_dir
```

Пример:

```bash
python3 multiple_alignment.py \
  -i results/BCR/grouped_by_germlines/vj \
  -o aligned_sequences
```

Если `mafft` не находится автоматически, можно указать путь к нему:

```bash
python3 multiple_alignment.py \
  -i /path/to/input_dir \
  -o /path/to/output_dir \
  -m /path/to/mafft
```

## Входные данные

На вход подаётся папка с FASTA-файлами. `multiple_alignment.py` ищет файлы рекурсивно,
то есть также обрабатывает FASTA-файлы во всех подпапках.

Поддерживаемые расширения: `.fasta`, `.fa`, `.fna`, `.ffn`, `.faa`.

Каждый FASTA-файл выравнивается отдельно.

## Что делает

Скрипт:

- находит все FASTA-файлы во входной папке;
- запускает MAFFT для каждого файла отдельно;
- сохраняет выровненные последовательности в выходную папку;
- сохраняет структуру подпапок из входной папки;
- создаёт файл `manifest.tsv` со списком входных и выходных файлов.

## Результаты

Для каждого входного FASTA-файла создаётся выровненный файл с суффиксом
`_aligned`.

```text
group1.fasta -> group1_aligned.fasta
```

В выходной папке также появится `manifest.tsv` — соответствие между исходными
FASTA-файлами и результатами выравнивания.
