"""Оркестратор end-to-end (блок B13).

Прогоняет весь ML-путь по всем группам входа и собирает выходы прогона:
загрузка → валидация → группировка → (на группу) MSA → IQ-TREE → мутации →
ранжирование → уверенные клады → артефакты группы; затем объединённая таблица
мутаций, глобальный список кандидатов, манифест и HTML-сводка.

Устойчивость: сбой одной группы не роняет прогон (статус ``failed`` + причина в
манифесте). Resume: при наличии готовых align.fasta / IQ-TREE-файлов тяжёлые шаги
(MAFFT/IQ-TREE) пропускаются, лёгкие (мутации/ранжирование) пересчитываются.
"""
from __future__ import annotations

from pathlib import Path

from . import align, grouping, io, report, tools, validate
from .clades import confident_clades
from .config import RunConfig
from .errors import BiocodeError, ToolNotFoundError
from .logging_ import configure, get_logger
from .model import Group, SequenceRecord, TreeResult
from .mutations import extract_mutations
from .ranking import rank_mutations
from .trees import iqtree, mrbayes
from .trees.compare import compare_trees, posterior_by_signature

log = get_logger("pipeline")


def _bayes_tree(group: Group, msa: dict[str, str], gdir: Path,
                cfg: RunConfig) -> TreeResult | None:
    """MrBayes-дерево группы (или None, если Bayes выключен / бинарь не установлен)."""
    if not cfg.run_bayes:
        return None
    mb_dir = gdir / "mrbayes"
    con = mb_dir / f"{group.key}.nex.con.tre"
    try:
        if cfg.resume and con.is_file():
            _, rev = mrbayes.safe_names(list(msa))
            newick, sup = mrbayes.parse_con_tre(con, rev)
            if newick:
                return TreeResult(method="mrbayes", newick=newick, model="GTR+I+G",
                                  supports=sup,
                                  outgroup=group.outgroup.id if group.outgroup else None)
        return mrbayes.run_mrbayes(group, msa, mb_dir, cfg)
    except ToolNotFoundError as e:
        log.warning("группа %s: байесовский путь пропущен — %s", group.key, e)
        return None


def process_group(group: Group, cfg: RunConfig, groups_root: Path,
                  records_by_id: dict[str, SequenceRecord]) -> tuple[dict, list, list]:
    """Обработать одну группу. Возвращает (статус, мутации, кандидаты)."""
    gdir = groups_root / group.key
    gdir.mkdir(parents=True, exist_ok=True)
    st: dict = {"key": group.key, "size": group.size,
                "v_gene": group.v_gene, "j_gene": group.j_gene}
    try:
        align_fa = gdir / "align.fasta"
        if cfg.resume and align_fa.is_file():
            msa = io.read_fasta(align_fa)
            log.debug("resume: %s — MSA из кэша", group.key)
        else:
            msa = align.align_group(group, gdir, cfg)

        if validate.has_errors(validate.validate_alignment(msa)):
            st.update(status="failed", reason="длины в выравнивании не совпадают")
            return st, [], []
        track = align.region_track(group, msa, gdir)

        iq_dir = gdir / "iqtree"
        tr = iqtree.load_result(iq_dir, group, cfg) if cfg.resume else None
        if tr is None:
            tr = iqtree.run_iqtree(group, align_fa, iq_dir, cfg)

        muts = extract_mutations(group.key, tr.newick, tr.ancestral or {}, msa, track)
        cands = rank_mutations(muts, cfg.weights)

        # опциональный байесовский путь (MrBayes) + сравнение ML↔Bayes
        posterior = None
        btr = _bayes_tree(group, msa, gdir, cfg)
        if btr is not None:
            cmp = compare_trees(tr, btr)
            io.write_json(cmp, gdir / "compare.json")
            posterior = posterior_by_signature(btr)
            st["rf_norm"] = cmp["robinson_foulds"]["normalized_rf"]
            st["shared_clades"] = cmp["shared_clades"]

        clds = confident_clades(tr, muts, records_by_id, posterior=posterior)

        if cfg.make_plots:
            try:
                from . import viz
                viz.plot_tree(tr, records_by_id, muts, clds, gdir / "tree.png")
                viz.plot_region_track(track, gdir / "regions.png")
            except Exception as e:                       # графики вторичны — не валим группу
                log.warning("группа %s: график не построен (%s)", group.key, e)

        report.write_mutations(muts, gdir / "mutations.tsv")
        report.write_candidates(cands, gdir / "candidates.tsv")
        report.write_group_report(gdir / "report.json", group_key=group.key, tree=tr,
                                  n_records=group.size, n_mutations=len(muts),
                                  candidates=cands, clades=clds)
        st.update(status="ok", model=tr.model, n_mutations=len(muts),
                  n_candidates=len(cands), n_confident_clades=len(clds),
                  bayes=btr is not None)
        return st, muts, cands
    except BiocodeError as e:
        log.error("группа %s провалилась: %s", group.key, str(e).splitlines()[0])
        st.update(status="failed", reason=str(e).splitlines()[0])
        return st, [], []


def _tool_versions(cfg: RunConfig) -> dict:
    versions = {}
    for name, hint in (("mafft", cfg.mafft_bin), ("iqtree", cfg.iqtree_bin)):
        try:
            b = tools.find_tool(name, hint)
            versions[name] = tools.tool_version(b, name)
        except ToolNotFoundError as e:
            versions[name] = f"НЕ НАЙДЕН: {e}"
    return versions


def run(cfg: RunConfig) -> dict:
    """Полный прогон. Возвращает словарь totals."""
    cfg.validate()
    run_dir = cfg.run_dir
    run_dir.mkdir(parents=True, exist_ok=True)
    configure(cfg.log_level, run_dir / "run.log")
    groups_root = run_dir / "groups"

    versions = _tool_versions(cfg)
    log.info("инструменты: %s", versions)

    records = io.read_airr_tsv(cfg.input, cfg)

    # дедупликация по id (реальные BCR-данные иногда содержат повторы contig-id):
    # предупреждаем и оставляем первый — не роняем прогон из-за квирка данных.
    seen: set[str] = set()
    uniq: list[SequenceRecord] = []
    n_dup = 0
    for r in records:
        if r.id in seen:
            n_dup += 1
            continue
        seen.add(r.id)
        uniq.append(r)
    if n_dup:
        log.warning("удалено %d записей с повторяющимся sequence_id (оставлен первый)", n_dup)
    records = uniq

    issues = validate.validate_records(records, min_group_size=cfg.min_group_size)
    errs = [i for i in issues if i.level == "error"]
    warns = [i for i in issues if i.level == "warn"]
    if errs:
        for i in errs[:20]:
            log.error("%s", i)
        raise BiocodeError(f"вход не прошёл валидацию ({len(errs)} ошибок) — см. лог")
    log.info("валидация входа: %d предупреждений, ошибок нет", len(warns))

    groups = grouping.group_records(records, cfg)
    ok_groups, skip_groups = grouping.analyzable(groups, cfg)
    records_by_id = {r.id: r for r in records}

    statuses: list[dict] = [
        {"key": g.key, "size": g.size, "status": "skipped",
         "reason": f"size<{cfg.min_group_size}"} for g in skip_groups]

    all_muts: list = []
    all_cands: list = []

    def work(g: Group):
        return process_group(g, cfg, groups_root, records_by_id)

    if cfg.jobs > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=cfg.jobs) as ex:
            results = list(ex.map(work, ok_groups))
    else:
        results = [work(g) for g in ok_groups]

    for stt, muts, cands in results:
        statuses.append(stt)
        all_muts.extend(muts)
        all_cands.extend(cands)

    all_cands.sort(key=lambda c: (c.score, c.n_branches, c.max_support), reverse=True)
    totals = {
        "groups_total": len(groups),
        "groups_ok": sum(1 for s in statuses if s.get("status") == "ok"),
        "groups_skipped": sum(1 for s in statuses if s.get("status") == "skipped"),
        "groups_failed": sum(1 for s in statuses if s.get("status") == "failed"),
        "mutations": len(all_muts),
        "candidates": len(all_cands),
        "confident_clades": sum(s.get("n_confident_clades", 0) for s in statuses),
    }

    report.write_mutations(all_muts, run_dir / "mutations.tsv")
    report.write_candidates(all_cands, run_dir / "candidates.tsv")
    report.write_manifest(run_dir / "manifest.json", config=cfg.to_dict(),
                          tool_versions=versions, group_statuses=statuses, totals=totals)
    report.write_summary_html(run_dir / "summary.html", totals=totals,
                              group_statuses=statuses, top_candidates=all_cands)

    log.info("ГОТОВО: групп ok=%d skip=%d fail=%d | мутаций=%d | кандидатов=%d | клад=%d",
             totals["groups_ok"], totals["groups_skipped"], totals["groups_failed"],
             totals["mutations"], totals["candidates"], totals["confident_clades"])
    log.info("выход: %s", run_dir)
    return totals
