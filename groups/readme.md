# groups

Уверенные клады (задача Никиты «уверенные клады» из [../BioCad_repo.md](../BioCad_repo.md)):
клады с надёжной поддержкой, на которые можно опираться в выводах. Обёртка над
продакшн-движком [`../biocode/clades.py`](../biocode/clades.py) (вендорено без
изменений из `BIOCAD.bigchallenges@main`).

## Вход / выход

| | формат | путь по умолчанию |
|---|---|---|
| вход | nexus | `mrbayes/<группа>.nex.con.tre` (+ `.names.tsv`) — выход [../mrbayes/run_mrbayes.py](../mrbayes/run_mrbayes.py) |
| вход (опционально) | newick с aLRT/UFBoot в метках узлов | `trees/<группа>/<группа>.treefile` — выход `anotherpipeline/build_trees/build_trees.sh` (Денис, IQ-TREE) |
| выход | json | `groups/report.json` |

Два источника, любой можно включить/выключить, оба пишутся в один `report.json`:

- **mrbayes** (по умолчанию включён, `--mrbayes-dir mrbayes`) — критерий
  **posterior ≥ 0.95** на консенсусном дереве MrBayes.
- **iqtree** (по умолчанию выключен, включается `--iqtree-dir trees`) — критерий
  **UFBoot ≥ 95 И aLRT ≥ 80** на ML-дереве IQ-TREE; здесь напрямую вызывается
  `biocode.clades.confident_clades` без каких-либо изменений алгоритма.

## Запуск

```bash
python confident_clades_report.py                                   # mrbayes/ → report.json
python confident_clades_report.py --iqtree-dir trees                # + ML-путь Дениса
python confident_clades_report.py --posterior-min 0.9 --ufboot-min 90 --alrt-min 70
```

## Формат report.json

```json
{
  "<группа>": {
    "mrbayes": {"threshold": {"posterior_min": 0.95}, "clades": [ {...} ]},
    "iqtree":  {"threshold": {"ufboot_min": 95, "alrt_min": 80}, "clades": [ {...} ]}
  }
}
```

Каждая клада: `clade`, `size`, `leaves` (исходные id таксонов), `ufboot`, `alrt`,
`posterior`, `defining_mutations`, `defining_cdr`, `isotypes`, `days`,
`confident_both_models` — единая схема для обоих источников (поля, не относящиеся
к источнику, — `null`).

Требует Python-пакет `biopython`.
