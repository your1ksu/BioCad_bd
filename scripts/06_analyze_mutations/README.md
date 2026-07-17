# run_mutations.py

Анализирует антигенные рецепторы (антитела): находит мутации относительно гермлайна, аннотирует по FR/CDR регионам (IMGT нумерация), выводит таблицу мутаций.

## Использование

```bash
python3 run_mutations.py -i fasta_from_clades -o mutation_tables -r data/references
```

Входная директория должна содержать `.fasta`/`.fa`/`.fas` файлы с нуклеотидными последовательностями антител.

## Что делает

1. Форматирует BLAST-базы из референсных FASTA
2. Находит все FASTA-файлы во входной директории
3. Для каждого файла запускает IgBLAST с гермлайновыми базами V/D/J (human, IMGT)
4. Парсит выход IgBLAST: V/D/J гены, тип цепи, границы доменов FR1-CDR3, FR4
5. Вычисляет аминокислотные мутации из выравнивания нуклеотидов
6. Сохраняет результаты в подпапке `output_dir/<basename>/`:
   - `mutations.tsv` — таблица мутаций (по одной строке на мутацию)
   - `mutations_summary.tsv` — сводка по последовательности (счётчики по регионам)

## Выходные файлы

```
mutation_tables/
├── clade1/
│   ├── mutations.tsv
│   └── mutations_summary.tsv
├── clade2/
│   ├── mutations.tsv
│   └── mutations_summary.tsv
└── ...
```

### mutations.tsv

| Колонка | Описание |
|---|---|
| sequence_id | ID последовательности из FASTA |
| chain_type | Тип цепи (IGH/IGK/IGL) |
| v_gene | V-ген |
| d_gene | D-ген |
| j_gene | J-ген |
| region | Регион (FR1-IMGT, CDR1-IMGT, ..., FR4-IMGT) |
| aa_position | Позиция АК в регионе |
| imgt_position | IMGT позиция |
| ref_aa | Референсная АК |
| mut_aa | АК в последовательности |
| ref_codon | Референсный кодон |
| mut_codon | Мутированный кодон |
| query_nt_pos | Позиция в query последовательности |
| is_silent | yes/no (синонимная мутация) |

### mutations_summary.tsv

| Колонка | Описание |
|---|---|
| sequence_id | ID последовательности |
| chain_type | Тип цепи |
| v_gene, d_gene, j_gene | Гены |
| fr1_muts ... fr4_muts | Число несинонимных мутаций по регионам |
| total_muts | Всего несинонимных мутаций |

## Референсные базы

Ожидаются в папке, переданной аргументом `-r` (например `data/references/`):
- `all_V.fasta` — V-сегменты (IGHV, IGKV, IGLV)
- `all_D.fasta` — D-сегменты (IGH, IGK, IGL)
- `all_J.fasta` — J-сегменты
- `human_gl.aux` — вспомогательный файл IgBLAST

BLAST-базы (`.nsq`, `.nin`, `.nhr`...) создаются автоматически при первом запуске.

## Запуск тестов

```bash
pytest tests/test_analyze_mutations/test_analyze_mutations.py -v
```