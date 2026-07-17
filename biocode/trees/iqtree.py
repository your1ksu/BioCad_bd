"""ML-деревья через IQ-TREE.

Запуск: ModelFinder (`-m MFP`) + UFBoot (`-B`) + SH-aLRT (`-alrt`) +
реконструкция предков (`--asr`), укоренение в germline-outgroup (`-o`).
Парсинг: модель и лог-правдоподобие из ``.iqtree``, дерево из ``.treefile``,
опоры ветвей (UFBoot/aLRT) по биразбиениям, предки — в :mod:`biocode.ancestral`.
"""
from __future__ import annotations

import time
from pathlib import Path

from .. import tools
from ..config import RunConfig
from ..errors import ToolRunError
from ..logging_ import get_logger
from ..model import Group, TreeResult

log = get_logger("iqtree")


def _parse_model(iqtree_file: Path) -> tuple[str | None, float | None]:
    """Модель (по BIC) и лог-правдоподобие из <prefix>.iqtree."""
    model, logl = None, None
    if not iqtree_file.is_file():
        return model, logl
    for line in iqtree_file.read_text(errors="ignore").splitlines():
        s = line.strip()
        if s.startswith("Best-fit model according to BIC:"):
            model = s.split(":", 1)[1].strip()
        elif s.startswith("Log-likelihood of the tree:"):
            try:
                logl = float(s.split(":", 1)[1].split()[0])
            except (ValueError, IndexError):
                pass
    return model, logl


def parse_supports(newick: str) -> dict[str, dict]:
    """Опоры ветвей по биразбиениям: подпись клады (frozenset листьев) → {alrt,ufboot}.

    IQ-TREE с -alrt и -B пишет метки внутренних узлов как 'SH-aLRT/UFboot' (напр.
    '98.5/100'). Ключ — стабильная подпись клады (отсортированные имена листьев),
    чтобы сопоставлять поддержку между деревьями (ML↔Bayes) в compare.py.
    """
    from Bio import Phylo
    from io import StringIO

    tree = Phylo.read(StringIO(newick), "newick")
    all_leaves = {t.name for t in tree.get_terminals()}
    supports: dict[str, dict] = {}
    for clade in tree.get_nonterminals():
        leaves = frozenset(t.name for t in clade.get_terminals())
        if len(leaves) < 2 or len(leaves) >= len(all_leaves):
            continue
        raw = clade.name
        if raw is None and clade.confidence is not None:
            raw = str(clade.confidence)
        if not raw:
            continue
        # метка внутреннего узла IQ-TREE: 'NodeName/aLRT/UFboot' (напр. 'Node2/100/100')
        # либо 'aLRT/UFboot', либо просто 'UFboot' — берём хвостовые числовые токены.
        nums: list[float] = []
        for tok in str(raw).split("/"):
            try:
                nums.append(float(tok))
            except ValueError:
                pass
        sig = "|".join(sorted(leaves))
        if len(nums) >= 2:
            supports[sig] = {"alrt": nums[-2], "ufboot": nums[-1]}
        elif len(nums) == 1:
            supports[sig] = {"ufboot": nums[0]}
    return supports


def load_result(out_dir: Path, group: Group, cfg: RunConfig) -> TreeResult | None:
    """Собрать TreeResult из уже посчитанных файлов IQ-TREE (для resume), без запуска."""
    prefix = Path(out_dir) / group.key
    treefile = Path(f"{prefix}.treefile")
    if not treefile.is_file():
        return None
    newick = treefile.read_text().strip()
    if not newick:
        return None
    model, _ = _parse_model(Path(f"{prefix}.iqtree"))
    ancestral = None
    state_file = Path(f"{prefix}.state")
    if cfg.asr and state_file.is_file():
        from ..ancestral import parse_state_file
        ancestral = parse_state_file(state_file)
    return TreeResult(method="iqtree", newick=newick, model=model,
                      supports=parse_supports(newick), ancestral=ancestral,
                      log_path=str(prefix) + ".log",
                      outgroup=group.outgroup.id if group.outgroup else None)


def run_iqtree(group: Group, align_fasta: Path, out_dir: Path, cfg: RunConfig) -> TreeResult:
    """Построить ML-дерево группы. Пишет вывод IQ-TREE в out_dir, возвращает TreeResult."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / group.key
    binpath = tools.find_tool("iqtree", cfg.iqtree_bin)

    cmd = [str(binpath), "-s", str(align_fasta),
           "-m", cfg.iqtree_model,
           "-B", str(cfg.ufboot), "-alrt", str(cfg.alrt),
           "-T", str(cfg.threads), "--seed", str(cfg.seed),
           "--prefix", str(prefix), "-redo"]
    if cfg.asr:
        cmd += ["--ancestral"]        # IQ-TREE 2/3: реконструкция предков (было -asr в 1.x)
    if group.outgroup is not None:
        cmd += ["-o", group.outgroup.id]

    log.info("IQ-TREE группа %s (model=%s, B=%d, alrt=%d, asr=%s)",
             group.key, cfg.iqtree_model, cfg.ufboot, cfg.alrt, cfg.asr)
    t0 = time.time()
    try:
        tools.run(cmd, timeout=cfg.timeout_s, tool="iqtree",
                  log_path=out_dir / f"{group.key}.run.log")
    except ToolRunError:
        # хвост родного лога IQ-TREE помогает диагностике
        native = prefix.with_suffix(".log")
        if native.is_file():
            log.error("IQ-TREE лог:\n%s", native.read_text(errors="ignore")[-1500:])
        raise
    runtime = time.time() - t0

    treefile = Path(f"{prefix}.treefile")
    newick = treefile.read_text().strip() if treefile.is_file() else ""
    model, logl = _parse_model(Path(f"{prefix}.iqtree"))

    ancestral = None
    state_file = Path(f"{prefix}.state")
    if cfg.asr and state_file.is_file():
        from ..ancestral import parse_state_file
        ancestral = parse_state_file(state_file)

    supports = parse_supports(newick) if newick else {}
    ver = tools.tool_version(binpath, "iqtree")
    log.info("IQ-TREE готово за %.1fs: model=%s, узлов с опорой=%d, предков=%s",
             runtime, model, len(supports), len(ancestral) if ancestral else 0)

    return TreeResult(method="iqtree", newick=newick, model=model, supports=supports,
                      ancestral=ancestral, tool_version=ver,
                      log_path=str(prefix) + ".log", runtime_s=runtime,
                      outgroup=group.outgroup.id if group.outgroup else None)
