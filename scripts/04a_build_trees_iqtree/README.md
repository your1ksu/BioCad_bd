# build_trees_iqtree.py

Скрипт строит филогенетические деревья по FASTA-выравниваниям через IQ-TREE.

## Использование

```bash
python3 build_trees_iqtree.py -i aligned_sequences -o trees
```

## Что делает

- Читает каждый FASTA-файл из указанной директории
- Определяет число CPU
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

Требует `iqtree` в `$PATH` (можно установить через `conda install -c bioconda iqtree`).

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