"""Сборка выходов (блок B11): таблица мутаций, ранжированные кандидаты, JSON, HTML.

Главный результат проекта — «таблица мутаций» и ранжированный список кандидатов
(потенциально важных замен). Здесь — их запись + манифест прогона и компактный
самодостаточный HTML-дашборд.
"""
from __future__ import annotations

import html
from pathlib import Path

from . import io
from .model import Mutation
from .ranking import Candidate

MUT_COLS = ["group", "branch", "position", "region", "ref", "alt", "kind",
            "branch_internal", "support"]
CAND_COLS = ["group", "position", "region", "ref", "alt", "kind",
             "n_branches", "n_internal", "max_support", "score",
             "contributions", "branches"]


def write_mutations(muts: list[Mutation], path: str | Path) -> Path:
    return io.write_tsv([m.as_row() for m in muts], path, columns=MUT_COLS)


def write_candidates(cands: list[Candidate], path: str | Path) -> Path:
    return io.write_tsv([c.as_row() for c in cands], path, columns=CAND_COLS)


def write_group_report(path: str | Path, *, group_key: str, tree, n_records: int,
                       n_mutations: int, candidates: list[Candidate],
                       clades: list[dict], top: int = 20) -> Path:
    obj = {
        "group": group_key,
        "n_records": n_records,
        "tree": {"method": tree.method, "model": tree.model,
                 "outgroup": tree.outgroup, "runtime_s": round(tree.runtime_s, 2),
                 "tool_version": tree.tool_version},
        "n_mutations": n_mutations,
        "n_candidates": len(candidates),
        "confident_clades": clades,
        "top_candidates": [c.as_row() for c in candidates[:top]],
    }
    return io.write_json(obj, path)


def write_manifest(path: str | Path, *, config: dict, tool_versions: dict,
                   group_statuses: list[dict], totals: dict) -> Path:
    return io.write_json({
        "config": config,
        "tool_versions": tool_versions,
        "totals": totals,
        "groups": group_statuses,
    }, path)


# ─────────────────────────────── HTML ─────────────────────────────────────

def _esc(x) -> str:
    return html.escape(str(x))


def write_summary_html(path: str | Path, *, totals: dict, group_statuses: list[dict],
                       top_candidates: list[Candidate], top: int = 40) -> Path:
    rows_groups = "".join(
        f"<tr><td>{_esc(g['key'])}</td><td>{_esc(g['status'])}</td>"
        f"<td class=n>{_esc(g.get('size',''))}</td>"
        f"<td class=n>{_esc(g.get('n_mutations',''))}</td>"
        f"<td class=n>{_esc(g.get('n_confident_clades',''))}</td>"
        f"<td>{_esc(g.get('reason',''))}</td></tr>"
        for g in group_statuses)
    rows_cand = "".join(
        f"<tr><td>{_esc(c.group)}</td><td class=n>{_esc(c.position)}</td>"
        f"<td><span class='reg {('cdr' if c.region.startswith('CDR') else 'fr')}'>{_esc(c.region)}</span></td>"
        f"<td>{_esc(c.ref)}→{_esc(c.alt)}</td><td>{_esc(c.kind)}</td>"
        f"<td class=n>{_esc(c.n_branches)}</td><td class=n>{c.max_support:.0f}</td>"
        f"<td class=n><b>{c.score:.2f}</b></td></tr>"
        for c in top_candidates[:top])
    doc = f"""<!doctype html><meta charset="utf-8"><title>BIOCODE — сводка прогона</title>
<style>
 body{{font-family:Manrope,system-ui,sans-serif;background:#0f1420;color:#e8edf3;margin:0;padding:22px}}
 h1{{font-size:20px;margin:0 0 2px}} .sub{{color:#9fb0c8;font-size:13px;margin-bottom:16px}}
 .cards{{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:18px}}
 .card{{background:#182034;border-radius:14px;padding:12px 16px;min-width:120px}}
 .card .k{{color:#9fb0c8;font-size:12px}} .card .v{{font-size:24px;font-weight:800}}
 h2{{font-size:15px;margin:18px 0 8px}}
 table{{border-collapse:collapse;width:100%;font-size:13px}}
 th,td{{text-align:left;padding:5px 10px;border-bottom:1px solid #223}} th{{color:#9fb0c8;position:sticky;top:0;background:#0f1420}}
 td.n{{text-align:right;font-variant-numeric:tabular-nums}}
 .reg{{padding:1px 7px;border-radius:9px;color:#10151f;font-weight:700;font-size:12px}}
 .reg.cdr{{background:#d94f3d}} .reg.fr{{background:#cfd9e6}}
 .wrap{{overflow-x:auto}}
</style>
<h1>BIOCODE — сводка прогона</h1>
<div class="sub">Филогенетический анализ эволюции антител · ML-путь (IQ-TREE)</div>
<div class="cards">
 <div class="card"><div class="k">групп проанализировано</div><div class="v">{_esc(totals.get('groups_ok',0))}</div></div>
 <div class="card"><div class="k">пропущено/ошибок</div><div class="v">{_esc(totals.get('groups_skipped',0)+totals.get('groups_failed',0))}</div></div>
 <div class="card"><div class="k">мутаций всего</div><div class="v">{_esc(totals.get('mutations',0))}</div></div>
 <div class="card"><div class="k">кандидатов</div><div class="v">{_esc(totals.get('candidates',0))}</div></div>
 <div class="card"><div class="k">уверенных клад</div><div class="v">{_esc(totals.get('confident_clades',0))}</div></div>
</div>
<h2>Топ кандидатов-мутаций (по скору значимости)</h2>
<div class="wrap"><table><tr><th>группа</th><th>поз.</th><th>регион</th><th>замена</th><th>тип</th><th>ветвей</th><th>UFBoot</th><th>скор</th></tr>{rows_cand}</table></div>
<h2>Группы</h2>
<div class="wrap"><table><tr><th>группа</th><th>статус</th><th>размер</th><th>мутаций</th><th>уверенных клад</th><th>причина</th></tr>{rows_groups}</table></div>
"""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(doc, encoding="utf-8")
    return p
