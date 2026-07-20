#!/usr/bin/env python3
"""Юнит-тесты analyze_mutations БЕЗ igblast/mb — проверяют СОДЕРЖИМОЕ таблиц.

Регрессия, которую они ловят: is_silent должен быть булевым на всём пути
(compute_mutations -> write_mutations/write_summary). Если его записать строкой
'yes'/'no', то `'yes' if m['is_silent'] else 'no'` и `if not m['is_silent']`
ломаются молча: колонка is_silent всегда 'yes', а summary весь в нулях.
Существующий tests/test_analyze_mutations проверяет только НАЛИЧИЕ файлов и такой
баг не увидит.
"""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "scripts"))
sys.path.insert(0, str(BASE / "scripts" / "06_analyze_mutations"))

import analyze_mutations as am  # noqa: E402

# AAG(K) TTT(F) GGA(G)  vs  AAA(K) TAT(Y) GGG(G)
# -> K->K silent, Y->F non-silent, G->G silent
QSEQ = "AAGTTTGGA"
SSEQ = "AAATATGGG"
DB = {"FR1-IMGT": {"qstart": 1, "qend": 9}}


def _result():
    return [{
        "query_id": "s1", "chain_type": "VH",
        "v_gene": "IGHV1", "d_gene": "-", "j_gene": "IGHJ1",
        "domain_boundaries": DB,
        "v_hit": {"qseq": QSEQ, "sseq": SSEQ, "qstart": 1},
        "j_hit": None,
    }]


def test_is_silent_is_bool():
    muts = am.compute_mutations(QSEQ, SSEQ, 1, DB)
    assert muts, "ожидались мутации"
    for m in muts:
        assert isinstance(m["is_silent"], bool), (
            f"is_silent должен быть bool, а не {type(m['is_silent']).__name__}")


def test_mutations_tsv_silent_column(tmp_path):
    out = tmp_path / "mutations.tsv"
    am.write_mutations(_result(), str(out))
    rows = [l.rstrip("\n").split("\t") for l in out.read_text().splitlines()]
    header, data = rows[0], rows[1:]
    idx = header.index("is_silent")
    silent_by_change = {(r[header.index("ref_aa")], r[header.index("mut_aa")]): r[idx]
                        for r in data}
    assert silent_by_change[("K", "K")] == "yes"     # синонимная
    assert silent_by_change[("Y", "F")] == "no"      # несинонимная — НЕ должна быть 'yes'


def test_summary_counts_nonsynonymous(tmp_path):
    out = tmp_path / "mutations_summary.tsv"
    am.write_summary(_result(), str(out))
    rows = [l.rstrip("\n").split("\t") for l in out.read_text().splitlines()]
    header, row = rows[0], rows[1]
    d = dict(zip(header, row))
    assert d["fr1_muts"] == "1", f"ожидалась 1 несинонимная в FR1, получено {d['fr1_muts']}"
    assert d["total_muts"] == "1", f"ожидался total=1, получено {d['total_muts']}"


if __name__ == "__main__":
    import tempfile
    test_is_silent_is_bool()
    with tempfile.TemporaryDirectory() as d:
        test_mutations_tsv_silent_column(Path(d))
        test_summary_counts_nonsynonymous(Path(d))
    print("OK: все юнит-тесты analyze_mutations прошли")
