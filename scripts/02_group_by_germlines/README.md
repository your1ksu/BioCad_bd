# group_by_germlines.py

Группировка отфильтрованных BCR-последовательностей по наиболее похожим germline-генам (V, D, J).

## Вход и выход

| | Путь | Что это |
|---|---|---|
| вход (`-i`) | `<out>/BCR_data_filtered.tsv` | результат предыдущего шага (`filter_sequences.py`) |
| выход (`-o`) | `<out>/grouped_by_germlines/` | выходные `.fasta` (подпапки `v/`, `d/`, `j/`, `vj/`) |
| референсы (`-r`) | `data/references/` | IMGT-гермлайновые справочники |

## Задача

Скрипт берёт отфильтрованный файл (результат `filter_sequences.py`) и для каждой последовательности находит наиболее похожий V-, D- (только для локуса `IGH`) и J-ген методом попарного локального выравнивания (Biopython `PairwiseAligner`).

**Результат:** набор `.fasta`-файлов, разложенных по подпапкам:
- `v/` — сгруппировано по V-гену
- `d/` — сгруппировано по D-гену (только тяжёлая цепь)
- `j/` — сгруппировано по J-гену
- `vj/` — сгруппировано по паре V+J

## Запуск

```bash
python3 group_by_germlines.py \
  -i results/BCR/BCR_data_filtered.tsv \
  -o results/BCR/grouped_by_germlines \
  -r data/references
```