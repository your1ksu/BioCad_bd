"""Извлечение мутаций на ветвях дерева (эволюционных событий).

Дерево укоренено в germline-outgroup (IQ-TREE ставит его в базальную развилку
``.treefile``), поэтому направление parent→child = предок→потомок = ход
соматического гипермутагенеза. Для каждой ветви сравниваем предковую
последовательность (ASR внутреннего узла) с потомком → список замен с контекстом
FR/CDR и типом (synonymous/replacement).

Метки внутренних узлов IQ-TREE имеют вид ``NodeName/aLRT/UFboot`` — ``NodeName``
совпадает с именами узлов в ``.state`` (см. :mod:`biocode.ancestral`), что и даёт
привязку предковых последовательностей к топологии.
"""
from __future__ import annotations

from io import StringIO

from .logging_ import get_logger
from .model import Mutation

log = get_logger("mutations")

_GAP = set("-.")


def _node_name(clade) -> str | None:
    """Имя узла IQ-TREE из метки клады ('Node3/0/29' → 'Node3'); лист → его id."""
    if clade.is_terminal():
        return clade.name
    if clade.name:
        return clade.name.split("/")[0]
    return None


def _support(clade) -> dict[str, float]:
    if clade.is_terminal() or not clade.name:
        return {}
    nums = []
    for tok in clade.name.split("/"):
        try:
            nums.append(float(tok))
        except ValueError:
            pass
    if len(nums) >= 2:
        return {"alrt": nums[-2], "ufboot": nums[-1]}
    if len(nums) == 1:
        return {"ufboot": nums[0]}
    return {}


def _classify(child_aligned: str, col: int, ref_char: str) -> str:
    """synonymous | replacement | unknown — по кодону в рамке потомка.

    Рамка: ungapped sequence_alignment начинается с границы кодона (FR1 = IMGT 1),
    поэтому индекс остатка // 3 = индекс кодона. Предковый кодон реконструируем,
    заменив производный нуклеотид на предковый (ref_char) в позиции внутри кодона.
    """
    r = sum(1 for ch in child_aligned[:col] if ch not in _GAP)   # индекс остатка в потомке
    ungapped = "".join(ch for ch in child_aligned if ch not in _GAP)
    cs = r - r % 3
    child_codon = ungapped[cs:cs + 3]
    if len(child_codon) != 3 or "N" in child_codon:
        return "unknown"
    pos = r % 3
    ref_codon = child_codon[:pos] + ref_char + child_codon[pos + 1:]
    if "N" in ref_codon:
        return "unknown"
    try:
        from Bio.Seq import Seq
        aa_child = str(Seq(child_codon).translate())
        aa_ref = str(Seq(ref_codon).translate())
    except Exception:
        return "unknown"
    return "synonymous" if aa_child == aa_ref else "replacement"


def extract_mutations(group_key: str, newick: str, ancestral: dict[str, str],
                      msa: dict[str, str], region_track: list[dict]) -> list[Mutation]:
    """Все замены на ветвях дерева (parent→child)."""
    from Bio import Phylo

    if not newick or not ancestral:
        return []
    tree = Phylo.read(StringIO(newick), "newick")

    # какие имена внутренних узлов присутствуют в метках дерева
    named = {nk for cl in tree.get_nonterminals() if (nk := _node_name(cl))}
    # корень базальной развилки в .treefile обычно без метки → берём оставшийся ASR-узел
    leftover = [n for n in ancestral if n not in named]

    def seq_of(clade) -> str | None:
        if clade.is_terminal():
            return msa.get(clade.name)
        nk = _node_name(clade)
        if nk is None and leftover:
            nk = leftover[0]
        return ancestral.get(nk) if nk else None

    def region_at(c: int) -> str:
        return region_track[c]["region"] if 0 <= c < len(region_track) else "other"

    muts: list[Mutation] = []

    def walk(parent) -> None:
        pseq = seq_of(parent)
        p_nk = _node_name(parent) or "ROOT"
        for child in parent.clades:
            cseq = seq_of(child)
            c_nk = _node_name(child) or (child.name or "?")
            sup = _support(child)
            if pseq and cseq:
                width = min(len(pseq), len(cseq))
                for c in range(width):
                    a, b = pseq[c], cseq[c]
                    if a != b and a not in _GAP and b not in _GAP:
                        muts.append(Mutation(
                            group=group_key, branch=f"{p_nk}→{c_nk}",
                            position=c, ref=a, alt=b,
                            kind=_classify(cseq, c, a), region=region_at(c),
                            support=sup, branch_internal=not child.is_terminal()))
            walk(child)

    walk(tree.root)
    log.info("группа %s: извлечено %d мутаций на ветвях", group_key, len(muts))
    return muts
