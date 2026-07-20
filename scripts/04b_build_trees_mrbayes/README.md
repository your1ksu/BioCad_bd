# 04b_build_trees_mrbayes — байесовские деревья (MrBayes)

Строит по одному байесовскому дереву на каждую группу выравнивания.

## Вход / выход

| | формат | путь по умолчанию |
|---|---|---|
| вход | `*.fasta`/`*.aln` (по одному на группу V+J, выровненные) | `aligned_sequences/` |
| выход | NEXUS-задание | `mrbayes/<группа>.nex` |
| выход | консенсусное дерево | `mrbayes/<группа>.nex.con.tre` |
| выход | карта имён | `mrbayes/<группа>.names.tsv` (`safe_id → исходный id`) |

`.names.tsv` нужен, потому что MrBayes не принимает `-`/спецсимволы в именах
таксонов (обычные в 10x-баркодах, напр. `GTTTCTATCATTATCC-1_contig_1`): таксоны
временно переименовываются в `T0000…`. По этому файлу
`../05_clade_search/clade_search.py` возвращает исходные id в отчёт.

## Запуск

```bash
python build_trees_mrbayes.py aligned_sequences --out mrbayes
python build_trees_mrbayes.py aligned_sequences --out mrbayes --nexus-only   # только .nex
python build_trees_mrbayes.py aligned_sequences --out mrbayes --outgroup <seq_id>
```

Особенности реализации:
- **Изоляция cwd**: каждый `mb` запускается в своей подпапке `_work/<группа>/`,
  поэтому параллельные процессы не мешают друг другу общими файлами.
- **stoprule**: MCMC останавливается по достижении сходимости
  (avg split freq < `--stopval`, по умолчанию 0.01); `--mb-ngen` — верхний предел.
- **Параллелизм** по группам на всех ядрах (`--workers`).

## CPU / GPU

По умолчанию — CPU (параллельно по группам). Для GPU (BEAGLE-CUDA) укажите бинарь:

```bash
python build_trees_mrbayes.py aligned_sequences --out mrbayes \
    --gpu-mb-bin /path/to/mrbayes-gpu/bin/mb --gpu-min-taxa 60
```

Группы с числом таксонов ≥ `--gpu-min-taxa` считаются на GPU последовательно
(одна карта), мелкие — параллельно на CPU. Без `--gpu-mb-bin` всё идёт на CPU.
GPU-бинарь MrBayes с BEAGLE-CUDA собирается отдельно (в conda его нет).

Требует `biopython` (парсинг `.nex.con.tre`) и бинарь `mb`
(`conda install -c bioconda mrbayes`).

## Дальше по конвейеру

`mrbayes/<группа>.nex.con.tre` → [../05_clade_search/clade_search.py](../05_clade_search/clade_search.py)
(уверенные клады по posterior ≥ 0.95).
