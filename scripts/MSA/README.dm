# MSA_final.py

Скрипт выполняет множественное выравнивание FASTA-файлов с помощью MAFFT.

## Входные данные

На вход подаётся папка с FASTA-файлами. Скрипт ищет файлы рекурсивно, то есть
также обрабатывает FASTA-файлы во всех подпапках.

Поддерживаемые расширения:

```text
.fasta, .fa, .fna, .ffn, .faa
```

Каждый FASTA-файл выравнивается отдельно.

## Запуск

Если `mafft` доступен в активном окружении:

```bash
python3 MSA_final.py --input /path/to/input_dir --output /path/to/output_dir
```

Короткий вариант:

```bash
python3 MSA_final.py -i /path/to/input_dir -o /path/to/output_dir
```

Если `mafft` не находится автоматически, укажите путь к нему:

```bash
python3 MSA_final.py \
  -i /path/to/input_dir \
  -o /path/to/output_dir \
  -m /path/to/mafft
```

Например:

```bash
python3 MSA_final.py \
  -i /home/hellstrom/Загрузки/grouped_by_germlines \
  -o /home/hellstrom/Загрузки/aligned_sequences \
  -m /home/hellstrom/.conda/envs/bio_env/bin/mafft
```

## Что делает скрипт

Скрипт:

- находит все FASTA-файлы во входной папке;
- запускает MAFFT для каждого файла отдельно;
- сохраняет выровненные последовательности в выходную папку;
- сохраняет структуру подпапок из входной папки;
- создаёт файл `manifest.tsv` со списком входных и выходных файлов.

## Результаты

Для каждого входного файла создаётся файл с суффиксом `_aligned`.

Пример:

```text
group1.fasta -> group1_aligned.fasta
```

В выходной папке также появится:

```text
manifest.tsv
```

В нём записано соответствие между исходными FASTA-файлами и результатами
выравнивания.
