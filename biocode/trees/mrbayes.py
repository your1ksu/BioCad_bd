"""Байесовские деревья через MrBayes.

Генерируем NEXUS (DATA + MRBAYES-блок с GTR+I+G, outgroup, MCMC), запускаем ``mb``,
парсим консенсус ``.con.tre`` (апостериорные вероятности клад) и проверяем
сходимость (средний стандартный разброс частот разбиений < 0.01).

Имена таксонов санитизируются (``T0001``…) — MrBayes капризен к '-'/спецсимволам;
после парсинга имена возвращаются к исходным id (для сопоставления с ML-деревом).

MrBayes не установлен по умолчанию → ``conda install -c bioconda mrbayes``. Обёртка
написана prod-качества и покрыта тестами на фикстурах; live-прогон включается при
наличии бинаря.
"""
from __future__ import annotations

import re
import time
from io import StringIO
from pathlib import Path

from .. import tools
from ..config import RunConfig
from ..logging_ import get_logger
from ..model import Group, TreeResult

log = get_logger("mrbayes")

_SPLIT_RE = re.compile(r"Average standard deviation of split frequencies:\s*([0-9.]+)")


def safe_names(ids: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    """id → безопасное имя (T0001…) и обратный маппинг."""
    fwd = {sid: f"T{i:04d}" for i, sid in enumerate(ids)}
    rev = {v: k for k, v in fwd.items()}
    return fwd, rev


def build_nexus(msa: dict[str, str], fwd: dict[str, str],
                outgroup_safe: str | None, cfg: RunConfig) -> str:
    """Собрать NEXUS для MrBayes из MSA."""
    ids = list(msa)
    width = len(next(iter(msa.values())))
    ngen = cfg.mb_ngen
    samplefreq = max(1, ngen // 1000)
    diagnfreq = max(1, ngen // 10)
    printfreq = max(1, ngen // 10)

    lines = ["#NEXUS", "", "begin data;",
             f"  dimensions ntax={len(ids)} nchar={width};",
             "  format datatype=DNA gap=- missing=? interleave=no;",
             "  matrix"]
    for sid in ids:
        lines.append(f"  {fwd[sid]}  {msa[sid]}")
    lines += ["  ;", "end;", "", "begin mrbayes;",
              # seed/swapseed через set (в mcmc параметр seed объявлен deprecated)
              f"  set autoclose=yes nowarn=yes seed={cfg.seed} swapseed={cfg.seed};",
              "  lset nst=6 rates=invgamma;                 [ GTR+I+G ]"]
    if outgroup_safe:
        lines.append(f"  outgroup {outgroup_safe};")
    lines += [
        f"  mcmc ngen={ngen} samplefreq={samplefreq} nchains=4 nruns=2 "
        f"printfreq={printfreq} diagnfreq={diagnfreq} "
        f"starttree=random savebrlens=yes;",
        f"  sumt burninfrac={cfg.mb_burnin_frac} conformat=simple;",
        f"  sump burninfrac={cfg.mb_burnin_frac};",
        "end;", ""]
    return "\n".join(lines)


def parse_con_tre(path: str | Path, rev: dict[str, str]) -> tuple[str, dict[str, dict]]:
    """``.con.tre`` (NEXUS) → (newick с исходными id, supports {clade_sig: {posterior}})."""
    from Bio import Phylo

    path = Path(path)
    if not path.is_file():
        return "", {}
    tree = next(Phylo.parse(str(path), "nexus"), None)
    if tree is None:
        return "", {}
    for tip in tree.get_terminals():
        tip.name = rev.get(tip.name, tip.name)

    all_leaves = {t.name for t in tree.get_terminals()}
    supports: dict[str, dict] = {}
    for clade in tree.get_nonterminals():
        leaves = frozenset(t.name for t in clade.get_terminals())
        if len(leaves) < 2 or len(leaves) >= len(all_leaves):
            continue
        if clade.confidence is not None:
            pp = float(clade.confidence)
            pp = pp / 100.0 if pp > 1.0 else pp        # проценты → доля
            supports["|".join(sorted(leaves))] = {"posterior": round(pp, 4)}
    out = StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue().strip(), supports


def convergence(stdout: str) -> tuple[float | None, bool]:
    """Средний станд. разброс частот разбиений из вывода mb; сошлось если < 0.01."""
    vals = _SPLIT_RE.findall(stdout or "")
    if not vals:
        return None, False
    last = float(vals[-1])
    return last, last < 0.01


def run_mrbayes(group: Group, msa: dict[str, str], out_dir: Path,
                cfg: RunConfig) -> TreeResult:
    """Построить байесовское дерево группы. Пишет вывод MrBayes в out_dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    binpath = tools.find_tool("mb", cfg.mrbayes_bin)

    fwd, rev = safe_names(list(msa))
    outg = fwd.get(group.outgroup.id) if group.outgroup else None
    nexus = build_nexus(msa, fwd, outg, cfg)
    nex_path = out_dir / f"{group.key}.nex"
    nex_path.write_text(nexus, encoding="utf-8")

    log.info("MrBayes группа %s: ngen=%d, %d таксонов", group.key, cfg.mb_ngen, len(msa))
    t0 = time.time()
    proc = tools.run([str(binpath), nex_path.name], cwd=out_dir,
                     timeout=cfg.timeout_s, tool="mb")
    (out_dir / f"{group.key}.mb.log").write_text(proc.stdout or "", encoding="utf-8")
    runtime = time.time() - t0

    newick, supports = parse_con_tre(out_dir / f"{group.key}.nex.con.tre", rev)
    avg_std, converged = convergence(proc.stdout or "")
    if not converged:
        log.warning("MrBayes группа %s: не сошлось (avg std split freq=%s ≥ 0.01) — "
                    "увеличьте mb_ngen", group.key, avg_std)

    return TreeResult(method="mrbayes", newick=newick, model="GTR+I+G",
                      supports=supports, tool_version=tools.tool_version(binpath, "mb"),
                      log_path=str(out_dir / f"{group.key}.mb.log"), runtime_s=runtime,
                      outgroup=group.outgroup.id if group.outgroup else None)
