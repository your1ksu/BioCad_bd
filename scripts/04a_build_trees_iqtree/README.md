# build_trees_iqtree.sh

Скрипт строит филогенетические деревья по FASTA-выравниваниям через IQ-TREE.

## Требования

- [Miniconda](https://docs.anaconda.com/miniconda/) (или Anaconda)

## Использование

```bash
./build_trees_iqtree.sh                           # input: aligned_sequences, output: trees
./build_trees_iqtree.sh aligned_sequences         # явный input, output: trees
./build_trees_iqtree.sh aligned_sequences my_trees # явные input и output
```

Скрипт сам создаст conda-окружение `trees_building_env`, установит `iqtree` и запустит его для всех `.fa/.fasta/.fas/.aln` файлов из указанной папки.

## Что делает

- Читает каждый FASTA-файл из указанной директории
- Определяет число CPU (`nproc` / `sysctl -n hw.ncpu`)
- Запускает IQ-TREE с параметрами:

```
-m MFP       ModelFinder Plus (подбор модели)
-B 1000      UFBoot2 (1000 реплик)
-T <NPROC>   все доступные ядра
-redo        перезапись старых файлов
```

- Результаты сохраняет в `<output_dir>/<имя_файла>/`

## Вход/выход

```
aligned_sequences/  →  trees/
  file1.fasta           file1/
  file2.fasta           file2/
```

Вся работа происходит внутри conda-окружения `trees_building_env`. Окружение создаётся один раз, при повторных запусках переиспользуется.

## Тесты

```bash
python3 -m pytest tests/test_build_trees_iqtree/test_build_trees_iqtree.py -v
```

Тесты:
- `test_successful_run` — прогон на синтетических данных (`test_build_iqtree_data/` с 2 файлами по 5 последовательностей), проверяет создание `test_result/` с подпапками для каждого fasta
- `test_nonexistent_input_dir` — проверка ошибки при несуществующей входной папке

Тестовые данные в `test_build_iqtree_data/`:
- `family1.fasta` — 5 последовательностей
- `family2.fasta` — 5 последовательностей