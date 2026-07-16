# BioCad_bd — филогенетический конвейер по антителам (BCR repertoire)

Конвейер обработки B-клеточных рецепторов (BCR): от сырых прочтений до
эволюционных деревьев, уверенных клад и таблиц мутаций. Каждый этап —
самостоятельный скрипт одного участника команды, соединённые в единую цепочку
через файлы на диске (fasta → fasta → nexus/newick → json/html/tsv).

## Конвейер

| # | Этап | Автор | Папка в `scripts/` | Вход | Выход |
|---|---|---|---|---|---|
| 1 | Фильтрация сиквенсов | Ксюша | `filtered/` | `data/<key>/BCR_data.tsv` | `results/<key>/BCR_data_filtered.tsv` |
| 2 | Группировка по гермлайнам | Ксюша | `group_by_germlines/` | `BCR_data_filtered.tsv` | `results/<key>/grouped_by_germlines/{v,d,j,vj}/*.fasta` |
| 3 | Множественное выравнивание (MSA) | Алина | `MSA/` | `grouped_by_germlines/vj/*.fasta` | `aligned_sequences/*_aligned.fasta` + `manifest.tsv` |
| 4 | ML-деревья (IQ-TREE) | Денис | `build_trees_iqtree/` | `aligned_sequences/` | `trees/<группа>/<группа>.treefile` (UFBoot+SH-aLRT) |
| 5a | NEXUS + байесовское дерево (MrBayes) | Никита | `mrbayes/` | `aligned_sequences/*_aligned.fasta` | `mrbayes/<группа>.nex`, `.nex.con.tre`, `.mb.log` |
| 5b | Уверенные клады | Никита | `groups/` | `mrbayes/*.nex.con.tre` и/или `trees/*/*.treefile` | `groups/report.json` |
| 6 | Визуализация деревьев | Денис | `visualize_trees/` | `trees/` или `mrbayes/` | `<out>/<группа>/<группа>.html` (toytree, интерактивный SVG) |
| 7 | Анализ мутаций | Денис | `analyze_mutations/` | fasta по кладам (IgBLAST) | `<out>/<клада>/mutations.tsv`, `mutations_summary.tsv` |

Шаги 1–2 (Ксюша) работают по единому ключу через `scripts/paths.py`
(`data/<key>/…` → `results/<key>/…`); шаги 3–7 принимают входную/выходную
папку явными аргументами командной строки — см. «Известные несоответствия»
ниже.

## Как прогнать весь конвейер (из корня репозитория)

```bash
# 1–2: фильтрация + группировка (ключ = имя подпапки в data/ и results/)
python3 scripts/filtered/filter_and_save.py --key BCR
python3 scripts/group_by_germlines/group_by_germlines.py --key BCR

# 3: выравнивание (MAFFT)
python3 scripts/MSA/MSA_final.py \
  -i results/BCR/grouped_by_germlines/vj \
  -o aligned_sequences

# 4: ML-деревья (IQ-TREE) — опционально, независимо от шага 5
./scripts/build_trees_iqtree/build_trees_iqtree.sh aligned_sequences trees

# 5a-b: байесовское дерево + уверенные клады (Никита)
python3 scripts/mrbayes/run_mrbayes.py aligned_sequences --out mrbayes
python3 scripts/groups/confident_clades_report.py \
  --mrbayes-dir mrbayes --iqtree-dir trees --out groups/report.json

# 6: визуализация деревьев
./scripts/visualize_trees/visualize_trees.sh trees trees_visualization
./scripts/visualize_trees/visualize_trees.sh mrbayes mrbayes_visualization

# 7: анализ мутаций (нужны fasta-файлы по кладам — см. «Известные пробелы»)
./scripts/analyze_mutations/analyze_mutations.sh fasta_from_clades mutation_tables
```

## Структура репозитория

```
data/<key>/BCR_data.tsv        — сырой вход (AIRR-подобный TSV)
data/IMGT/<вид>/{IG,TR}/*.fasta — germline-справочники IMGT (используется Homo_sapiens/IG)
results/<key>/…                 — результаты шагов 1-2 (Ксюша)
scripts/<этап>/                 — код каждого этапа + README на этапе
tests/                          — тесты (pytest для большинства этапов,
                                   консольные для mrbayes/groups — см. tests/README.md)
misc/                           — пусто (зарезервировано)
```

## Требования

- Python 3.9+, [Miniconda](https://docs.anaconda.com/miniconda/)
- Python-пакеты: `pandas`, `biopython`
- Внешние инструменты (ставятся conda-окружениями автоматически там, где
  скрипты это умеют — `trees_building_env` для IQ-TREE/визуализации):
  `mafft`, `iqtree`, `mb` (MrBayes, `conda install -c bioconda mrbayes`),
  `igblastn` + BLAST+ (для анализа мутаций), `toytree`/`toyplot` (визуализация)

## Тестирование

Большинство этапов покрыто `pytest` (`tests/test_*`). Этап Никиты
(`mrbayes/` + `groups/`) тестируется отдельными консольными скриптами без
pytest — см. [tests/README.md](tests/README.md) и [tests/TEST_REPORT.md](tests/TEST_REPORT.md)
(там же разобран найденный и исправленный баг в извлечении уверенных клад).

## Известные пробелы и несоответствия

- **`biocode/` не опубликован.** `scripts/mrbayes/run_mrbayes.py` и
  `scripts/groups/confident_clades_report.py` импортируют пакет `biocode`
  (вендорен из `BIOCAD.bigchallenges@main`) из родительской директории — без
  него оба скрипта падают на импорте. См. `scripts/README.md`.
- **Два разных соглашения о путях.** Шаги 1–2 (Ксюша) используют
  централизованный `scripts/paths.py` с ключом (`data/<key>` ↔
  `results/<key>`); шаги 3–7 — обычные `--input`/`--output`/позиционные
  аргументы с папками по умолчанию в текущей директории
  (`aligned_sequences`, `trees`, `mrbayes`, `groups`). Единого раннера,
  который бы прокидывал результат одного шага во вход следующего
  автоматически, нет — команды выше нужно запускать по очереди вручную.
- **Шаг 7 (анализ мутаций) ожидает fasta по кладам**, а `groups/report.json`
  (выход шага 5b) отдаёт клады как списки id внутри JSON, а не как отдельные
  fasta-файлы — недостающее звено: скрипт, который бы разбивал
  `aligned_sequences/*.fasta` на `fasta_from_clades/<клада>.fasta` по данным
  `report.json`, пока не написан.
- `.gitignore` и `gitignore.txt` — два файла с частично разным содержимым;
  активный (учитываемый git) — `.gitignore`, `gitignore.txt` — вероятно,
  забытый черновик.

## Участники

| Участник | Этапы |
|---|---|
| Ксюша | фильтрация, группировка по гермлайнам |
| Алина | множественное выравнивание (MSA) |
| Денис | ML-деревья (IQ-TREE), визуализация деревьев, анализ мутаций |
| Никита | NEXUS/MrBayes (байесовское дерево), уверенные клады |
