# MSA.sh

`MSA.sh` запускает скрипт `MSA_final.py` и автоматически подготавливает всё
необходимое окружение для работы с MAFFT.

## Что делает скрипт

Скрипт:

- находит `MSA_final.py`;
- проверяет, установлен ли `conda`;
- если `conda` не найден, скачивает и устанавливает локальный Miniforge;
- создаёт conda-окружение `msa_final_env`;
- устанавливает в окружение `python=3.12` и `mafft`;
- запускает `MSA_final.py` с переданными аргументами.

## Требования

На устройстве должен быть установлен только Python 3.

Остальное скрипт скачает и установит сам:

- Miniforge/conda, если conda отсутствует;
- MAFFT;
- отдельное conda-окружение для запуска.

Для первого запуска нужен интернет.

## Запуск

Перейдите в папку со скриптом:

```bash
cd /home/hellstrom/Документы/MACSE
```

Запустите выравнивание:

```bash
./MSA.sh -i /path/to/input_fastas -o /path/to/output_dir
```

Где:

- `-i` или `--input` — папка с входными FASTA-файлами;
- `-o` или `--output` — папка, куда будут сохранены результаты.

Пример:

```bash
./MSA.sh \
  -i /home/hellstrom/Загрузки/grouped_by_germlines \
  -o /home/hellstrom/Загрузки/aligned_sequences
```

## Входные файлы

`MSA_final.py` ищет FASTA-файлы рекурсивно во входной папке.

Поддерживаемые расширения:

```text
.fasta, .fa, .fna, .ffn, .faa
```

## Результат

Для каждого входного FASTA-файла создаётся выровненный файл с суффиксом
`_aligned`.

Пример:

```text
group1.fasta -> group1_aligned.fasta
```

Также в выходной папке создаётся файл:

```text
manifest.tsv
```

В нём записано соответствие между исходными FASTA-файлами и результатами
выравнивания.

## Повторный запуск

При повторном запуске окружение не создаётся заново, если оно уже существует.
Скрипт просто использует готовое окружение `msa_final_env` и запускает
`MSA_final.py`.

## Дополнительные настройки

Можно изменить имя conda-окружения:

```bash
MSA_ENV_NAME=my_msa_env ./MSA.sh -i input_fastas -o aligned_output
```

Можно указать другой путь к `MSA_final.py`:

```bash
MSA_PY=/path/to/MSA_final.py ./MSA.sh -i input_fastas -o aligned_output
```

Можно указать другую папку для локальной установки Miniforge:

```bash
MSA_MINIFORGE_DIR=/path/to/miniforge ./MSA.sh -i input_fastas -o aligned_output
```
