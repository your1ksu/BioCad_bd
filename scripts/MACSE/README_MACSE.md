# MACSE MSA

Короткий модуль пайплайна для множественного выравнивания FASTA-файлов с помощью MACSE.

Скрипт ищет все FASTA-файлы во входной папке, запускает MACSE для каждого файла и сохраняет:

- `*_aligned.fasta` - нуклеотидное выравнивание;
- `*_aligned_aa.fasta` - аминокислотное выравнивание;
- `manifest.tsv` - таблицу соответствия входных и выходных файлов.

## Требования

- Linux/macOS shell;
- `python3`;
- `conda` или доступ к интернету для автоматической установки Miniforge;
- MACSE устанавливается автоматически в conda-окружение `macse_env`.

## Запуск

Рекомендуемый запуск через shell-обертку:

```bash
./run_macse.sh \
  -i /path/to/input_fastas \
  -o /path/to/output_dir
```

Пример для V-генов:

```bash
./run_macse.sh \
  -i /home/hellstrom/Загрузки/grouped_by_germlines/v \
  -o /home/hellstrom/Документы/MACSE/v
```

Можно задать имя conda-окружения:

```bash
./run_macse.sh \
  -i /path/to/input_fastas \
  -o /path/to/output_dir \
  -e my_macse_env
```

`MSA.sh` делает то же самое и просто вызывает `run_macse.sh`:

```bash
./MSA.sh -i /path/to/input_fastas -o /path/to/output_dir
```

## Что делает run_macse.sh

1. Ищет `conda`.
2. Если `conda` нет, ставит локальный Miniforge в `.miniforge`.
3. Создает conda-окружение, если его еще нет.
4. Устанавливает `macse`.
5. Запускает `MACSEtry.py` внутри этого окружения.

## Запуск Python-скрипта напрямую

Если окружение уже активировано и `macse` доступен:

```bash
python3 MACSEtry.py \
  --input /path/to/input_fastas \
  --output /path/to/output_dir
```

Если нужно указать MACSE явно:

```bash
python3 MACSEtry.py \
  --input /path/to/input_fastas \
  --output /path/to/output_dir \
  --macse /path/to/macse
```

## Тесты

```bash
python3 -m unittest discover -s tests_MACSE -p 'test_*.py'
```

Тесты не запускают настоящий MACSE: внешний вызов подменяется фейковой функцией, чтобы проверить логику скрипта быстро и без зависимости от окружения.

## Важно

MACSE предназначен для кодирующих нуклеотидных последовательностей. Если входные последовательности не являются coding DNA, имеют проблемную рамку считывания или много стоп-кодонов, MACSE может завершиться с ошибкой.
