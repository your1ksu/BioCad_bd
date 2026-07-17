"""Множественное выравнивание группы (MAFFT).

Запуск MAFFT без ``shell=True`` (кросс-платформенность). На вход — ungapped
рабочие последовательности группы + outgroup; на выход — MSA (dict id→строка)
и файл. FR/CDR-метки проецируются на колонки MSA через annotate.project_onto_msa.
"""
from __future__ import annotations

from pathlib import Path

from . import annotate, io, tools
from .config import RunConfig
from .logging_ import get_logger
from .model import Group

log = get_logger("align")

_MODE_ARGS = {
    "linsi": ["--maxiterate", "1000", "--localpair"],
    "ginsi": ["--maxiterate", "1000", "--globalpair"],
    "einsi": ["--maxiterate", "1000", "--genafpair"],
    "auto": ["--auto"],
}


def align_group(group: Group, out_dir: Path, cfg: RunConfig) -> dict[str, str]:
    """Выровнять группу (+ outgroup). Пишет input.fasta и align.fasta, возвращает MSA."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    recs = group.all_records()                       # записи + outgroup
    in_fa = io.write_fasta(recs, out_dir / "input.fasta")
    out_fa = out_dir / "align.fasta"

    mafft = tools.find_tool("mafft", cfg.mafft_bin)
    mode_args = _MODE_ARGS.get(cfg.mafft_mode, _MODE_ARGS["auto"])
    cmd = [str(mafft), *mode_args, "--quiet", "--thread", "-1", str(in_fa)]
    log.info("MAFFT (%s) группа %s: %d последовательностей", cfg.mafft_mode, group.key, len(recs))
    proc = tools.run(cmd, timeout=cfg.timeout_s, tool="mafft")

    out_fa.write_text(proc.stdout, encoding="utf-8")
    (out_dir / "align.mafft.log").write_text(proc.stderr or "", encoding="utf-8")
    msa = io.read_fasta(out_fa)
    log.info("MAFFT готово: %d строк × %d колонок", len(msa),
             max((len(s) for s in msa.values()), default=0))
    return msa


def region_track(group: Group, msa: dict[str, str], out_dir: Path) -> list[dict]:
    """FR/CDR-трек по колонкам MSA (из per-record меток) + запись align.regions.tsv."""
    labeled = {r.id: (r.meta.get("labels") or [])
               for r in group.records if r.meta.get("labels")}
    track = annotate.project_onto_msa(labeled, msa)
    io.write_tsv([{"column": t["column"], "region": t["region"],
                   "purity": t["purity"], "coverage": t["coverage"]} for t in track],
                 Path(out_dir) / "align.regions.tsv",
                 columns=["column", "region", "purity", "coverage"])
    return track