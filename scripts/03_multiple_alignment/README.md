# MSA: выравнивание FASTA-файлов через MAFFT

`multiple_alignment.py` выполняет множественное выравнивание FASTA-файлов с помощью
MAFFT.

Для удобного запуска есть shell-скрипт `multiple_alignment.sh`, который сам создаёт conda-
окружение и устанавливает MAFFT.

## Быстрый запуск через multiple_alignment.sh

Рекомендуемый способ запуска:

```bash
cd /home/hellstrom/Документы/MACSE
./multiple_alignment.sh -i /path/to/input_fastas -o /path/to/output_dir
```

Пример:

```bash
cd /home/hellstrom/Документы/MACSE
./multiple_alignment.sh \
  -i /home/hellstrom/Загрузки/grouped_by_germlines \
  -o /home/hellstrom/Загрузки/aligned_sequences
```

Аргументы:

- `-i` или `--input` — папка с входными FASTA-файлами;
- `-o` или `--output` — папка, куда будут сохранены результаты.

## Что делает multiple_alignment.sh

`multiple_alignment.sh`:

- находит `multiple_alignment.py`;
- проверяет, установлен ли `conda`;
- если `conda` не найден, скачивает и устанавливает локальный Miniforge;
- создаёт conda-окружение `msa_final_env`;
- устанавливает в окружение `python=3.12` и `mafft`;
- запускает `multiple_alignment.py` с переданными аргументами.

На устройстве должен быть установлен только Python 3. Для первого запуска нужен
интернет, чтобы скачать Miniforge и MAFFT.

При повторном запуске окружение не создаётся заново, если оно уже существует.
Скрипт использует готовое окружение `msa_final_env`.

## Входные данные

На вход подаётся папка с FASTA-файлами. `multiple_alignment.py` ищет файлы рекурсивно,
то есть также обрабатывает FASTA-файлы во всех подпапках.

Поддерживаемые расширения:

```text
.fasta, .fa, .fna, .ffn, .faa
```

Каждый FASTA-файл выравнивается отдельно.

## Что делает multiple_alignment.py

Скрипт:

- находит все FASTA-файлы во входной папке;
- запускает MAFFT для каждого файла отдельно;
- сохраняет выровненные последовательности в выходную папку;
- сохраняет структуру подпапок из входной папки;
- создаёт файл `manifest.tsv` со списком входных и выходных файлов.

## Результаты

Для каждого входного FASTA-файла создаётся выровненный файл с суффиксом
`_aligned`.

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

## Ручной запуск multiple_alignment.py

Если `mafft` уже доступен в активном окружении:

```bash
python3 multiple_alignment.py --input /path/to/input_dir --output /path/to/output_dir
```

Короткий вариант:

```bash
python3 multiple_alignment.py -i /path/to/input_dir -o /path/to/output_dir
```

Если `mafft` не находится автоматически, можно указать путь к нему:

```bash
python3 multiple_alignment.py \
  -i /path/to/input_dir \
  -o /path/to/output_dir \
  -m /path/to/mafft
```

Пример:

```bash
python3 multiple_alignment.py \
  -i /home/hellstrom/Загрузки/grouped_by_germlines \
  -o /home/hellstrom/Загрузки/aligned_sequences \
  -m /home/hellstrom/.conda/envs/bio_env/bin/mafft
```

## Дополнительные настройки multiple_alignment.sh

Можно изменить имя conda-окружения:

```bash
MSA_ENV_NAME=my_msa_env ./multiple_alignment.sh -i input_fastas -o aligned_output
```

Можно указать другой путь к `multiple_alignment.py`:

```bash
MSA_PY=/path/to/multiple_alignment.py ./multiple_alignment.sh -i input_fastas -o aligned_output
```

Можно указать другую папку для локальной установки Miniforge:

```bash
MSA_MINIFORGE_DIR=/path/to/miniforge ./multiple_alignment.sh -i input_fastas -o aligned_output
```
