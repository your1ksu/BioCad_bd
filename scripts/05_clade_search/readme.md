# 05_clade_search

Поиск **уверенных клад** — под-линий клона с надёжной поддержкой ветви, на
которые можно опираться в выводах. Читает деревья предыдущего шага, отбирает
клады по порогу поддержки и пишет единый `report.json`. По желанию сохраняет
FASTA каждой клады (вход для шага мутаций).

## Вход / выход

| | формат | путь |
|---|---|---|
| вход (Байес) | NEXUS-консенсус + карта имён | `mrbayes/<группа>.nex.con.tre` (+ `<группа>.names.tsv`) — выход `04b_build_trees_mrbayes` |
| вход (ML, опц.) | newick с `aLRT/UFboot` в метках узлов | `trees/<группа>/<группа>.treefile` — выход `04a_build_trees_iqtree` |
| выход | JSON | `report.json` |
| выход (опц.) | FASTA по кладам | `--clades-fasta-dir` (вход шага `06_analyze_mutations`) |

Два источника, любой включается/выключается, оба пишутся в один `report.json`:

- **mrbayes** (по умолчанию, `--mrbayes-dir`) — критерий **posterior ≥ 0.95** на
  консенсусном дереве MrBayes.
- **iqtree** (опц., `--iqtree-dir`) — критерий **UFBoot ≥ 95 И aLRT ≥ 80** на
  ML-дереве IQ-TREE.

## Запуск

```bash
python clade_search.py --mrbayes-dir mrbayes --out report.json
python clade_search.py --mrbayes-dir mrbayes --iqtree-dir trees --out report.json
python clade_search.py --mrbayes-dir mrbayes \
    --aligned-dir aligned_sequences --clades-fasta-dir clades   # + FASTA клад
```

Пороги настраиваются: `--posterior-min`, `--ufboot-min`, `--alrt-min`.
Требует Python-пакет `biopython`.

## Формат `report.json`

```json
{
  "<группа>": {
    "mrbayes": {"threshold": {"posterior_min": 0.95}, "clades": [ {клада}, ... ]},
    "iqtree":  {"threshold": {"ufboot_min": 95, "alrt_min": 80}, "clades": [ {клада}, ... ]}
  }
}
```

Схема одной клады (порядок ключей стабилен):

| поле | смысл |
|---|---|
| `clade` | подпись клады (внутреннее имя узла или `id\|id\|…`) |
| `size` | число листьев в кладе |
| `leaves` | список исходных id таксонов |
| **`depth`** | **число внутренних узлов внутри клады** (не считая сам узел-предок). Пара листьев от общего предка → `0`; каждый уровень ветвления внутри → `+1` |
| **`ancestor_to_leaves`** | **дистанция от предка клады до её листьев** в длинах ветвей (замен на сайт): `{"max", "mean", "min"}` |
| `posterior` | апостериорная вероятность (источник `mrbayes`) |
| `ufboot`, `alrt` | опоры ML-ветви (источник `iqtree`) |
| `defining_mutations`, `defining_cdr`, `isotypes`, `days` | зарезервированы (заполняются на других шагах) |
| `confident_both_models` | согласие ML↔Bayes (по умолчанию `false`) |

### Поля `depth` и `ancestor_to_leaves` — для аналитики деревьев

Следующий шаг `tree_analytics/tree_analytics.py` берёт из каждой клады `mrbayes`
уже готовые (не пересчитывает) `size`, `depth` и `ancestor_to_leaves.mean` и
строит по ним распределения. Без этих полей его панели глубины и дистанции
остаются пустыми, поэтому они рассчитываются здесь и обязательны в выходе.

- **`depth`** считается по РЕАЛЬНОМУ дереву, а не как `size − 1`: majority-rule
  консенсус MrBayes бывает политомным, и формула по размеру завышала бы глубину.
- **`ancestor_to_leaves.mean`** — средняя генетическая дистанция (subs/site) от
  предка клады до листьев; `max`/`min` дают разброс.

## Устойчивость к особенностям деревьев MrBayes

- **Комплемент корня.** MrBayes пишет неукоренённое дерево как rooted-newick с
  базовой политомией — часть клады может висеть «голыми» листьями на корне вместо
  отдельного узла. Такая сторона разбиения достраивается с той же posterior
  (`_root_complement_supports`), иначе симметричная клада терялась бы.
- **Дедупликация сторон ребра.** Если клада и её дополнение обе присутствуют и
  разного размера — большая сторона это «всё, кроме меньшей клады»; она убирается
  (`_drop_complement_duplicates`), при равном размере обе стороны сохраняются.
