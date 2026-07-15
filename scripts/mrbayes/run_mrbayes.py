#!/usr/bin/env python3
"""ШАГ «Никиты» (nexus + MRBayes) — байесовские деревья по группам.

Обёртка над продакшн-движком ``biocode.trees.mrbayes`` (см. ../biocode/,
вендорено без изменений из BIOCAD.bigchallenges@main), адаптированная под
файловый контракт остальных участников конвейера (см. ../BioCad_repo.md —
это ДВЕ отдельные строки таблицы: «nexus» (fasta→nexus) и «MRBayes»
(nexus→nexus-дерево), выполняются здесь как две фазы одного скрипта):

  вход:  aligned_sequences/*.fasta — множественные выравнивания по группам
         V+J (выход шага Алины «MSA»; по умолчанию берём ту же папку, что и
         anotherpipeline/build_trees/build_trees.sh у Дениса — единый вход
         для ML- и Bayes-путей)
  выход: mrbayes/<группа>.nex          — NEXUS (DATA + MRBAYES-блок, GTR+I+G)
         mrbayes/<группа>.nex.con.tre  — консенсусное дерево (фаза 2, нужен mb)
         mrbayes/<группа>.mb.log       — лог MrBayes (фаза 2)
         mrbayes/<группа>.names.tsv    — safe_id → исходный id (см. ниже)

Фаза 1 (генерация .nex) не требует бинаря mb и всегда отрабатывает — так
шаг «nexus» из таблицы остаётся самостоятельным даже если MrBayes ещё не
установлен. Фаза 2 запускает ``mb`` на уже написанном .nex, если бинарь
найден; иначе печатает подсказку и не падает (аналогично поведению
``biocode.pipeline._bayes_tree`` — Bayes-путь опционален).

MrBayes капризен к именам таксонов (не любит '-'/спецсимволы, обычные в
10x-баркодах вида ``GTTTCTATCATTATCC-1_contig_1``), поэтому таксоны
переименовываются в безопасные ``T0001``... Генерация NEXUS и чтение
консенсуса — два независимых запуска (может быть даже на разных машинах,
раз MrBayes долгий), поэтому маппинг сохраняется рядом в
``<группа>.names.tsv``, чтобы ``groups/confident_clades_report.py`` мог
вернуть исходные id.

Требует бинарь ``mb`` для фазы 2: conda install -c bioconda mrbayes.
Требует Python-пакет biopython (парсинг .nex.con.tre в фазе 2).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from biocode import tools
from biocode.config import RunConfig
from biocode.errors import ToolNotFoundError, ToolRunError
from biocode.trees import mrbayes

FASTA_EXTS = (".fa", ".fasta", ".fas", ".aln")


def read_fasta(path: Path) -> dict[str, str]:
    """FASTA → {id: seq}. Гэпы '-' сохраняются как есть — это уже MSA Алины."""
    seqs: dict[str, str] = {}
    cur = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].strip().split()[0]
            seqs[cur] = ""
        elif cur is not None:
            seqs[cur] += line.strip()
    return seqs


def group_key_of(fasta_path: Path) -> str:
    name = fasta_path.stem
    if name.endswith("_aligned"):
        name = name[: -len("_aligned")]
    return name


def write_nexus(key: str, msa: dict[str, str], outgroup_id: str | None,
                cfg: RunConfig, out_dir: Path) -> tuple[Path, dict[str, str]]:
    """Фаза 1: MSA → <группа>.nex (+ .names.tsv). Не требует бинаря mb."""
    fwd, rev = mrbayes.safe_names(list(msa))
    names_tsv = out_dir / f"{key}.names.tsv"
    names_tsv.write_text(
        "\n".join(f"{safe}\t{orig}" for orig, safe in fwd.items()) + "\n",
        encoding="utf-8")
    outg = fwd.get(outgroup_id) if outgroup_id else None
    nexus_text = mrbayes.build_nexus(msa, fwd, outg, cfg)
    nex_path = out_dir / f"{key}.nex"
    nex_path.write_text(nexus_text, encoding="utf-8")
    return nex_path, rev


def run_mb(key: str, nex_path: Path, rev: dict[str, str], cfg: RunConfig,
          out_dir: Path) -> str:
    """Фаза 2: запустить mb на уже готовом .nex. Возвращает статус-строку."""
    try:
        binpath = tools.find_tool("mb", cfg.mrbayes_bin)
    except ToolNotFoundError as e:
        print(f"[{key}] mb не найден — nexus сохранён ({nex_path.name}), "
              f"MrBayes-фаза пропущена: {e}", file=sys.stderr)
        return "nexus_only"

    t0 = time.time()
    try:
        proc = tools.run([str(binpath), nex_path.name], cwd=out_dir,
                         timeout=cfg.timeout_s, tool="mb")
    except ToolRunError as e:
        print(f"[{key}] ОШИБКА mb — {e}", file=sys.stderr)
        return "failed"
    (out_dir / f"{key}.mb.log").write_text(proc.stdout or "", encoding="utf-8")
    runtime = time.time() - t0

    newick, supports = mrbayes.parse_con_tre(out_dir / f"{key}.nex.con.tre", rev)
    avg_std, converged = mrbayes.convergence(proc.stdout or "")
    print(f"[{key}] готово: {out_dir / (key + '.nex.con.tre')} "
          f"(runtime={runtime:.1f}s, клад с posterior={len(supports)}, "
          f"сошлось={'да' if converged else 'нет'} avg_std={avg_std})")
    return "ok"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("input_dir", nargs="?", default="aligned_sequences",
                    help="папка с *.fasta (одно выравнивание = одна группа V+J); "
                         "по умолчанию 'aligned_sequences'")
    ap.add_argument("--out", default="mrbayes", help="выходная папка (по умолчанию 'mrbayes')")
    ap.add_argument("--outgroup", default=None,
                    help="id таксона-outgroup (germline-предок), если он есть в fasta; "
                         "по умолчанию не укореняем")
    ap.add_argument("--mb-ngen", type=int, default=200_000, help="длина цепи MCMC")
    ap.add_argument("--mb-burnin-frac", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--mb-bin", default="", help="явный путь к бинарю mb (по умолчанию — автопоиск в $PATH)")
    ap.add_argument("--timeout-s", type=int, default=3600)
    ap.add_argument("--nexus-only", action="store_true",
                    help="только фаза 1 (сгенерировать .nex), не пытаться запускать mb")
    args = ap.parse_args(argv)

    in_dir = Path(args.input_dir)
    out_dir = Path(args.out)
    if not in_dir.is_dir():
        print(f"Ошибка: папка не найдена: {in_dir}", file=sys.stderr)
        return 1

    fasta_files = sorted(p for p in in_dir.iterdir() if p.suffix.lower() in FASTA_EXTS)
    if not fasta_files:
        print(f"В {in_dir} не найдено *.fa/*.fasta/*.fas/*.aln", file=sys.stderr)
        return 0

    cfg = RunConfig(mb_ngen=args.mb_ngen, mb_burnin_frac=args.mb_burnin_frac,
                    seed=args.seed, mrbayes_bin=args.mb_bin, timeout_s=args.timeout_s)

    out_dir.mkdir(parents=True, exist_ok=True)
    counts = {"ok": 0, "nexus_only": 0, "failed": 0, "skipped": 0}
    for fasta in fasta_files:
        key = group_key_of(fasta)
        msa = read_fasta(fasta)
        if len(msa) < 3:
            print(f"[{key}] пропущено: меньше 3 последовательностей ({len(msa)})")
            counts["skipped"] += 1
            continue

        nex_path, rev = write_nexus(key, msa, args.outgroup, cfg, out_dir)
        print(f"[{key}] nexus готов: {nex_path}")

        status = "nexus_only" if args.nexus_only else run_mb(key, nex_path, rev, cfg, out_dir)
        counts[status] += 1

    print(f"\nИтого: ok={counts['ok']} только-nexus={counts['nexus_only']} "
          f"ошибок={counts['failed']} пропущено={counts['skipped']} из {len(fasta_files)} групп")
    return 0 if counts["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
