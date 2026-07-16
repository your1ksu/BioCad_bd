# mrbayes

Байесовские деревья по группам

## Вход / выход

| | формат | путь по умолчанию |
|---|---|---|
| вход | несколько `.fasta` (по одному на группу V+J, уже выровненные) | `aligned_sequences/` — тот же вход, что берёт `anotherpipeline/build_trees/build_trees.sh` (Денис, IQ-TREE) |
| выход | nexus | `mrbayes/<группа>.nex` |
| выход | nexus (дерево) | `mrbayes/<группа>.nex.con.tre`, `mrbayes/<группа>.mb.log` |

Дополнительно пишется `mrbayes/<группа>.names.tsv` (`safe_id → исходный id`) —
MrBayes не принимает `-`/спецсимволы в именах таксонов, которые обычны в
10x-баркодах (`GTTTCTATCATTATCC-1_contig_1`), поэтому таксоны временно
переименовываются в `T0001…`. Этот файл нужен `../groups/confident_clades_report.py`,
чтобы вернуть исходные id в отчёте.

## Запуск

```bash
python run_mrbayes.py                       # aligned_sequences/ → mrbayes/
python run_mrbayes.py путь/до/fasta --out mrbayes
python run_mrbayes.py --nexus-only           # только сгенерировать .nex, не запускать mb
python run_mrbayes.py --outgroup <seq_id>    # если в fasta есть germline-предок
```

Скрипт состоит из двух фаз:
1. **nexus** — MSA → `.nex` (DATA-блок + MRBAYES-блок: GTR+I+G, MCMC). Не требует
   бинаря `mb`, всегда отрабатывает.
2. **MRBayes** — запускает `mb` на только что написанном `.nex`, парсит
   `.nex.con.tre` (апостериорные вероятности клад). Если `mb` не найден в `$PATH`
   (`conda install -c bioconda mrbayes`), фаза 2 пропускается с понятным
   сообщением, а `.nex` остаётся на диске — можно прогнать `mb` вручную или
   перезапустить скрипт позже.

Требует Python-пакет `biopython` (парсинг `.nex.con.tre`).

## Дальше по конвейеру

`mrbayes/<группа>.nex.con.tre` — вход для [../groups/confident_clades_report.py](../groups/confident_clades_report.py)
(поиск уверенных клад по апостериорной поддержке, posterior ≥ 0.95).
