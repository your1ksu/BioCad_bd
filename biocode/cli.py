"""CLI: ``python -m biocode <команда>``.

Команды:
  run   — полный ML-прогон по AIRR-TSV (см. pipeline.run)
"""
from __future__ import annotations

import argparse
import sys

from .config import RunConfig
from .errors import BiocodeError
from .logging_ import get_logger

log = get_logger("cli")


def _add_run(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("run", help="полный ML-прогон по AIRR-TSV")
    p.add_argument("--input", required=True, help="AIRR-TSV (напр. EDU/data/BCR_data.tsv)")
    p.add_argument("--out", default="EDU/results/biocode_run", help="каталог прогона")
    p.add_argument("--working-seq", default="sequence_alignment")
    p.add_argument("--group-by", default="v_j_gene")
    p.add_argument("--locus", default="IGH")
    p.add_argument("--min-group-size", type=int, default=4)
    p.add_argument("--limit-groups", type=int, default=0, help="0 = без лимита")
    p.add_argument("--mafft-mode", default="linsi")
    p.add_argument("--ufboot", type=int, default=1000)
    p.add_argument("--alrt", type=int, default=1000)
    p.add_argument("--no-asr", action="store_true", help="отключить реконструкцию предков")
    p.add_argument("--bayes", action="store_true", help="также строить байесовское дерево (MrBayes) и сравнивать с ML")
    p.add_argument("--mb-ngen", type=int, default=200000, help="длина цепи MCMC для MrBayes")
    p.add_argument("--no-plots", action="store_true", help="не рисовать дерево/FR-CDR PNG")
    p.add_argument("--threads", default="AUTO")
    p.add_argument("--jobs", type=int, default=1, help="групп в параллель")
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--no-resume", action="store_true")
    p.add_argument("--log-level", default="INFO")


def _add_gui(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("gui", help="структурированный браузер результатов (UGENE-стиль)")
    p.add_argument("--run", required=True, help="каталог прогона <out>/run (с manifest.json)")
    p.add_argument("--port", type=int, default=8766)


def _cfg_from_run(a: argparse.Namespace) -> RunConfig:
    return RunConfig(
        input=a.input, out=a.out, working_seq=a.working_seq, group_by=a.group_by,
        locus=a.locus, min_group_size=a.min_group_size, limit_groups=a.limit_groups,
        mafft_mode=a.mafft_mode, ufboot=a.ufboot, alrt=a.alrt, asr=not a.no_asr,
        run_bayes=a.bayes, mb_ngen=a.mb_ngen, make_plots=not a.no_plots,
        threads=a.threads, jobs=a.jobs, seed=a.seed, resume=not a.no_resume,
        log_level=a.log_level,
    ).validate()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="biocode",
                                 description="BIOCODE — филогенетика и анализ мутаций антител")
    sub = ap.add_subparsers(dest="cmd", required=True)
    _add_run(sub)
    _add_gui(sub)
    args = ap.parse_args(argv)

    try:
        if args.cmd == "run":
            from . import pipeline
            totals = pipeline.run(_cfg_from_run(args))
            print(f"\nГОТОВО. Группы: ok={totals['groups_ok']} "
                  f"skip={totals['groups_skipped']} fail={totals['groups_failed']} | "
                  f"мутаций={totals['mutations']} | кандидатов={totals['candidates']} | "
                  f"уверенных клад={totals['confident_clades']}")
            return 0
        if args.cmd == "gui":
            from . import gui
            gui.launch(args.run, args.port)
            return 0
    except BiocodeError as e:
        log.error("%s", e)
        return 2
    return 1


if __name__ == "__main__":
    sys.exit(main())
