# nexus, MrBayes, уверенные клады

| задача | вход | выход | папка |
|---|---|---|---|
| nexus | fasta | nexus | `mrbayes/` |
| MRBayes | nexus | nexus (дерево) | `mrbayes/` |
| уверенные клады | nexus | json | `groups/` |

Обе папки — тонкие файловые обёртки над продакшн-движком `biocode`
(`biocode/trees/mrbayes.py`, `biocode/clades.py`), адаптированные под входной/
выходной формат остальных участников конвейера (Ксюша → Алина → [mrbayes/groups] → Денис).

- [mrbayes/readme.md](mrbayes/readme.md) — генерация NEXUS (GTR+I+G) и запуск MrBayes:
  `aligned_sequences/*.fasta` → `mrbayes/<группа>.nex`, `.nex.con.tre`, `.mb.log`.
- [groups/readme.md](groups/readme.md) — уверенные клады: `mrbayes/*.nex.con.tre`
  (posterior ≥ 0.95) и/или ML-дерево IQ-TREE Дениса (UFBoot ≥ 95, aLRT ≥ 80)
  → единый `groups/report.json`.

## Зависимость: пакет `biocode`

Оба скрипта импортируют `biocode` из родительской директории
(`sys.path.insert(..., parent.parent)` → `from biocode... import ...`). В этой
ветке папка `biocode/` **не публикуется** (пушились только `mrbayes/`, `groups/`
и этот readme) — при выкладке в отдельный репозиторий её нужно либо скопировать
рядом самостоятельно (вендорена без изменений из `BIOCAD.bigchallenges@main`,
модули `model.py`, `config.py`, `tools.py`, `errors.py`, `logging_.py`,
`clades.py`, `trees/`), либо поправить импорт под фактическое расположение
пакета в целевом репозитории. Без неё оба скрипта падают на импорте.

Также требуется Python-пакет `biopython` (парсинг NEXUS/Newick) и `mb` (`conda install -c bioconda mrbayes`).
