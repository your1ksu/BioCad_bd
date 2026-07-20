#!/usr/bin/env python3
"""
Конструктор синтетических BCR-данных из РЕАЛЬНЫХ зародышевых антител.

Идея: берём настоящие IMGT-зародышевые гены V(-D)-J из data/references,
собираем наивного предка клона (без «мусора» на концах — строго от начала
FR1 до конца FR4), затем накладываем контролируемые соматические мутации.
Для каждой мутации точно известно, в какой регион (FR1/CDR1/FR2/CDR2/FR3) она
попала — это независимая «истина» по IMGT-нумерации, с которой потом
сравнивается разметка пайплайна (IgBLAST + 06_analyze_mutations).

Выходы:
  synthetic_BCR.fasta   — чистые последовательности (вход для IgBLAST/пайплайна)
  ground_truth.tsv      — по одной строке на кодон-мутацию: регион, позиция,
                          зародышевая/наблюдаемая а.к., синонимичность
  groups.tsv            — какая группа какой germline получила
  summary.txt           — сводка + эхо конфига

Запуск:
  python make_synthetic.py -c synthetic_config.txt -r ../../data/references -o OUTDIR
"""

import argparse
import os
import random
from pathlib import Path


# --- Стандартный генетический код (как в 06_analyze_mutations) ---------------
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), 'X')


# --- Границы регионов по IMGT-нумерации (кодоны, 1-based) --------------------
# Совпадает с IMGT_DOMAIN_STARTS в 06_analyze_mutations/analyze_mutations.py.
def region_of_codon(k: int) -> str:
    if k <= 26:
        return 'FR1'
    if k <= 38:
        return 'CDR1'
    if k <= 55:
        return 'FR2'
    if k <= 65:
        return 'CDR2'
    if k <= 104:
        return 'FR3'
    if k <= 117:
        return 'CDR3'
    return 'FR4'


CDR_REGIONS = {'CDR1', 'CDR2'}
FR_REGIONS = {'FR1', 'FR2', 'FR3'}


# --- Чтение FASTA ------------------------------------------------------------
def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name = line[1:]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs


def parse_imgt_header(header):
    """'M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|...' -> (allele, functionality)."""
    parts = header.split('|')
    allele = parts[1] if len(parts) > 1 else header
    func = parts[3].strip() if len(parts) > 3 else ''
    return allele, func


def load_germlines(ref_dir, locus_prefix, seg):
    """Читает {IGHV,IGKV,...}.fasta -> {allele: gapped_seq} только функциональные (F)."""
    path = Path(ref_dir) / f"{locus_prefix}{seg}.fasta"
    out = {}
    if not path.exists():
        return out
    for header, seq in read_fasta(path).items():
        allele, func = parse_imgt_header(header)
        if func != 'F':                       # только функциональные аллели
            continue
        out[allele] = seq.lower()
    return out


# --- Карта «позиция в ungapped V -> регион» по IMGT-гэпам --------------------
def v_region_map(gapped_v):
    """
    Возвращает (ungapped_seq, region_at, codon_at):
      region_at[i] — регион нуклеотида i в ungapped-последовательности,
      codon_at[i]  — IMGT-кодон (1-based) этого нуклеотида.
    Каждые 3 колонки IMGT-выравнивания = один кодон (нумерация фиксирована,
    гэпы '.' занимают колонки, поэтому индекс кодона считается по колонке).
    """
    ungapped = []
    region_at = []
    codon_at = []
    for col, ch in enumerate(gapped_v):
        codon = col // 3 + 1
        if ch != '.':
            ungapped.append(ch)
            region_at.append(region_of_codon(codon))
            codon_at.append(codon)
    return ''.join(ungapped), region_at, codon_at


def is_complete_v(gapped_v):
    """V должен покрывать FR1..FR3 (до кодона 104 → колонка 312) и быть в рамке."""
    ungapped, region_at, _ = v_region_map(gapped_v)
    if len(gapped_v) < 306:                    # не доходит до конца FR3
        return False
    if len(ungapped) < 270 or len(ungapped) % 3 != 0:
        return False
    prot = ''.join(translate(ungapped[i:i + 3]) for i in range(0, len(ungapped), 3))
    return '*' not in prot                     # чистая рамка без стопов


# --- Сборка наивного предка --------------------------------------------------
NUC = 'ACGT'


def clean_codon(rng, avoid_stop=True):
    while True:
        c = ''.join(rng.choice(NUC) for _ in range(3))
        if not avoid_stop or translate(c) != '*':
            return c


def stops_in_frame(seq):
    seq = seq[:len(seq) - (len(seq) % 3)]
    return sum(1 for i in range(0, len(seq), 3) if translate(seq[i:i + 3]) == '*')


def j_coding_frame(j):
    """
    J-REGION кодирует конец CDR3 + FR4, но его рамка сдвинута (мотив [WF]G-x-G
    не с позиции 0). Выбираем рамку f∈{0,1,2} без стоп-кодонов, предпочитая ту,
    где виден мотив FR4. Возвращает j, обрезанный к этой рамке (кратно 3).
    """
    best = None
    for f in range(3):
        sub = j[f:]
        sub = sub[:len(sub) - (len(sub) % 3)]
        if not sub:
            continue
        prot = ''.join(translate(sub[i:i + 3]) for i in range(0, len(sub), 3))
        stops = prot.count('*')
        # мотив начала FR4: [WF]-G-x-G  (WGQG/WGRG у тяжёлых, FGxG у лёгких)
        has_motif = any(prot[k] in 'WF' and prot[k + 1] == 'G' and prot[k + 3] == 'G'
                        for k in range(len(prot) - 3))
        score = (stops, 0 if has_motif else 1, -len(sub))
        if best is None or score < best[0]:
            best = (score, sub)
    return best[1] if best else j


def build_ancestor(rng, v_ungapped, j_ungapped, d_seqs, is_heavy, mutate_fr4):
    """
    Наивный предок = V-REGION + CDR3-стык + J-REGION, строго в рамке (frame 0),
    без стоп-кодонов, без случайного «мусора» на концах.
    Возвращает (seq, v_len, j_start, j_len).
    CDR3-стык: для тяжёлых — кусок реального D + N-нуклеотиды; для лёгких — только N.
    Каждый сегмент кратен 3 и стыкуется по границе кодона → сквозная рамка 0.
    """
    v = v_ungapped                              # кратно 3, без стопов (is_complete_v)
    j = j_coding_frame(j_ungapped)              # в правильной рамке, кратно 3, без стопов
    for _ in range(200):
        n_codons = rng.randint(2, 4)            # 2-4 кодона N-области стыка
        junction = ''.join(clean_codon(rng) for _ in range(n_codons))
        if is_heavy and d_seqs:
            d = rng.choice(list(d_seqs.values()))
            if len(d) >= 6:
                dl = rng.choice([3, 6, 9])      # рамочный кусок реального D
                start = rng.randint(0, len(d) - dl)
                dseg = d[start:start + dl]
                junction = clean_codon(rng) + dseg + clean_codon(rng)
        if stops_in_frame(junction):            # D-кусок мог внести стоп
            continue
        seq = v + junction + j
        if '*' in ''.join(translate(seq[i:i + 3]) for i in range(0, len(seq), 3)):
            continue
        v_len = len(v)
        j_start = len(v) + len(junction)
        return seq, v_len, j_start, len(j)
    raise RuntimeError("Не удалось собрать предка без стоп-кодонов")


# --- Наложение мутаций с истиной по регионам ---------------------------------
def choose_positions(rng, region_at, v_len, n_mut, region_bias, mutate_fr4, total_len):
    """Выбирает n_mut различных позиций (0-based) с учётом смещения по регионам."""
    v_positions = list(range(v_len))
    if mutate_fr4:
        # FR4 = хвост J; позиции добавим, но они не в region_at (V-карта) — метим отдельно
        pass
    if region_bias == 'cdr':
        pool_pref = [p for p in v_positions if region_at[p] in CDR_REGIONS]
        pool_rest = [p for p in v_positions if region_at[p] in FR_REGIONS]
    elif region_bias == 'fr':
        pool_pref = [p for p in v_positions if region_at[p] in FR_REGIONS]
        pool_rest = [p for p in v_positions if region_at[p] in CDR_REGIONS]
    else:
        pool_pref, pool_rest = v_positions, []

    chosen = set()
    if pool_pref and (region_bias in ('cdr', 'fr')):
        n_pref = min(len(pool_pref), int(round(n_mut * 0.75)))
        chosen.update(rng.sample(pool_pref, n_pref))
    # добить остаток равномерно из всего V
    remaining = [p for p in v_positions if p not in chosen]
    need = n_mut - len(chosen)
    if need > 0 and remaining:
        chosen.update(rng.sample(remaining, min(need, len(remaining))))
    return sorted(chosen)


def random_other_base(base, rng):
    alt = [b for b in NUC if b != base.upper()]
    return rng.choice(alt) if alt else base


def mutate_at(rng, chars, p):
    """
    Меняет нуклеотид в позиции p на другой, НЕ создавая стоп-кодон: нонсенс-замены
    летальны для B-клетки, продуктивные антитела их не несут. Правит chars на месте.
    """
    cstart = (p // 3) * 3
    if cstart + 3 > len(chars):
        chars[p] = random_other_base(chars[p], rng)
        return
    codon = [chars[cstart], chars[cstart + 1], chars[cstart + 2]]
    off = p - cstart
    alts = [b for b in NUC if b != chars[p].upper()]
    rng.shuffle(alts)
    for b in alts:
        codon[off] = b
        if translate(''.join(codon)) != '*':
            chars[p] = b
            return
    # все три варианта дают стоп — крайне редко; позицию пропускаем


def mutate_seq(rng, seq, region_at, v_len, n, region_bias, mutate_fr4):
    """Новая строка: n различных не-нонсенс замен в V-регионе с учётом смещения."""
    positions = choose_positions(rng, region_at, v_len, n, region_bias,
                                 mutate_fr4, len(seq))
    chars = list(seq)
    for p in positions:
        mutate_at(rng, chars, p)
    return ''.join(chars)


def codon_mutations(germline, observed, v_len, region_at, codon_at):
    """
    Сравнивает germline и observed по КОДОНАМ V-региона (как это делает пайплайн:
    на уровне а.к.). Возвращает список dict истинных мутаций региона.
    """
    muts = []
    for cstart in range(0, v_len - 2, 3):
        g_cod = germline[cstart:cstart + 3].upper()
        o_cod = observed[cstart:cstart + 3].upper()
        if g_cod == o_cod:
            continue
        g_aa = translate(g_cod)
        o_aa = translate(o_cod)
        region = region_at[cstart]              # регион по первому нуклеотиду кодона
        imgt = codon_at[cstart]                 # IMGT-кодон (== aa-позиция)
        muts.append({
            'region': region,
            'imgt_codon': imgt,
            'nt_pos': cstart,
            'ref_codon': g_cod, 'obs_codon': o_cod,
            'ref_aa': g_aa, 'obs_aa': o_aa,
            'silent': 'yes' if g_aa == o_aa else 'no',
        })
    return muts


# --- Линия клона -------------------------------------------------------------
def gen_lineage(rng, ancestor, region_at, v_len, n_members, muts_per_leaf,
                region_bias, mutate_fr4):
    """
    Ветвящаяся линия: рекурсивно делим лист на поддеревья, на каждом ребре
    добавляя порцию мутаций. Бюджет = muts_per_leaf распределяется по пути
    так, что КАЖДЫЙ лист накапливает ровно ~muts_per_leaf замен (контроль
    нагрузки), при этом близкие листья делят общие «стволовые» мутации → клады.
    """
    def build(seq, budget, k):
        if k <= 1:
            return [mutate_seq(rng, seq, region_at, v_len, budget,
                               region_bias, mutate_fr4)]
        # часть бюджета — на общее ребро (делится всеми потомками)
        e = rng.randint(0, max(0, budget // 2)) if budget > 1 else 0
        shared = mutate_seq(rng, seq, region_at, v_len, e, region_bias, mutate_fr4)
        left = rng.randint(1, k - 1)
        return (build(shared, budget - e, left)
                + build(shared, budget - e, k - left))

    return build(ancestor, muts_per_leaf, n_members)[:n_members]


def gen_star(rng, ancestor, region_at, v_len, n_members, muts_per_leaf,
             region_bias, mutate_fr4):
    return [mutate_seq(rng, ancestor, region_at, v_len, muts_per_leaf,
                       region_bias, mutate_fr4)
            for _ in range(n_members)]


# --- Конфиг ------------------------------------------------------------------
def parse_config(path):
    cfg = {
        'groups': 5, 'sequences': 250, 'mutation_rate': 3.0,
        'germlines': 'auto', 'chain': 'mix', 'region_bias': 'none',
        'phylogeny': 'lineage', 'mutate_fr4': 'no', 'seed': 42,
    }
    for raw in open(path):
        line = raw.split('#', 1)[0].strip()
        if not line or '=' not in line:
            continue
        key, val = line.split('=', 1)
        cfg[key.strip()] = val.strip()
    cfg['groups'] = int(cfg['groups'])
    cfg['sequences'] = int(cfg['sequences'])
    cfg['mutation_rate'] = float(cfg['mutation_rate'])
    cfg['seed'] = int(cfg['seed'])
    cfg['mutate_fr4'] = str(cfg['mutate_fr4']).lower() in ('yes', 'true', '1')
    return cfg


LOCI = {
    'IGH': ('IGH', True),   # (префикс, тяжёлая?)
    'IGK': ('IGK', False),
    'IGL': ('IGL', False),
}


def pick_germlines(rng, ref_dir, cfg):
    """Возвращает список (group_name, locus, v_allele, v_gapped, j_allele, j_ungapped)."""
    chain = cfg['chain'].upper()
    wanted = cfg['germlines']
    # доступные локусы
    if chain == 'MIX':
        loci = ['IGH', 'IGK', 'IGL']
    elif chain in LOCI:
        loci = [chain]
    else:
        loci = ['IGH', 'IGK', 'IGL']

    # загрузим все нужные germline
    vbank, jbank = {}, {}
    for loc in loci:
        pref, _ = LOCI[loc]
        vbank[loc] = {a: s for a, s in load_germlines(ref_dir, pref, 'V').items()
                      if is_complete_v(s)}
        jbank[loc] = load_germlines(ref_dir, pref, 'J')

    picks = []
    if wanted != 'auto':
        for pair in wanted.split(','):
            pair = pair.strip()
            if not pair:
                continue
            vall, jall = [x.strip() for x in pair.split('/')]
            loc = vall[:3]
            vg = vbank.get(loc, {}).get(vall)
            jg = jbank.get(loc, {}).get(jall)
            if vg is None or jg is None:
                raise SystemExit(f"Не найдена пара {pair} среди функциональных germline")
            picks.append((loc, vall, vg, jall, jg))
        return picks

    # auto: раскидываем группы по локусам по кругу
    order = [loc for loc in ['IGK', 'IGH', 'IGL'] if vbank.get(loc) and jbank.get(loc)]
    if not order:
        raise SystemExit("Нет функциональных germline в data/references")
    for i in range(cfg['groups']):
        loc = order[i % len(order)]
        vall = rng.choice(sorted(vbank[loc]))
        jall = rng.choice(sorted(jbank[loc]))
        picks.append((loc, vall, vbank[loc][vall], jall, jbank[loc][jall]))
    return picks


def ungap_j(gapped_j):
    return gapped_j.replace('.', '')


def main():
    ap = argparse.ArgumentParser(description="Конструктор синтетических BCR из реальных germline.")
    ap.add_argument('-c', '--config', required=True)
    ap.add_argument('-r', '--references', default=None,
                    help="Каталог data/references (по умолчанию — относительно репозитория)")
    ap.add_argument('-o', '--output', required=True)
    args = ap.parse_args()

    cfg = parse_config(args.config)
    rng = random.Random(cfg['seed'])

    ref_dir = args.references
    if ref_dir is None:
        ref_dir = Path(__file__).resolve().parents[2] / 'data' / 'references'
    ref_dir = str(ref_dir)

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    d_seqs = {a: s for a, s in load_germlines(ref_dir, 'IGH', 'D').items()}
    if not d_seqs:
        d_seqs = {k: v for k, v in read_fasta(os.path.join(ref_dir, 'IGHD.fasta')).items()} \
            if os.path.exists(os.path.join(ref_dir, 'IGHD.fasta')) else {}
        d_seqs = {k: v.lower() for k, v in d_seqs.items()}

    picks = pick_germlines(rng, ref_dir, cfg)
    if len(picks) < cfg['groups']:
        cfg['groups'] = len(picks)

    # раскидать sequences по группам
    per = [cfg['sequences'] // cfg['groups']] * cfg['groups']
    for i in range(cfg['sequences'] % cfg['groups']):
        per[i] += 1

    fasta_lines = []
    gt_rows = ["sequence_id\tgroup\tv_gene\tj_gene\tregion\timgt_codon\tnt_pos\t"
               "ref_codon\tobs_codon\tref_aa\tobs_aa\tsilent"]
    groups_rows = ["group\tlocus\tv_allele\tj_allele\tn_sequences\tv_len\tancestor_len"]

    region_counter = {}
    seq_total = 0

    for gi, (loc, vall, vgap, jall, jgap) in enumerate(picks):
        v_ungapped, region_at, codon_at = v_region_map(vgap)
        j_ungapped = ungap_j(jgap)
        is_heavy = (loc == 'IGH')
        ancestor, v_len, j_start, j_len = build_ancestor(
            rng, v_ungapped, j_ungapped, d_seqs, is_heavy, cfg['mutate_fr4'])
        # приведём region_at/codon_at к фактической длине V в предке
        region_at = region_at[:v_len]
        codon_at = codon_at[:v_len]

        n_members = per[gi]
        muts_per_leaf = max(1, int(round(cfg['mutation_rate'] / 100.0 * v_len)))
        gname = f"G{gi+1:02d}_{vall.replace('*','_')}_{jall.replace('*','_')}"

        if cfg['phylogeny'] == 'star':
            observed_list = gen_star(rng, ancestor, region_at, v_len, n_members,
                                     muts_per_leaf, cfg['region_bias'], cfg['mutate_fr4'])
        else:
            observed_list = gen_lineage(rng, ancestor, region_at, v_len, n_members,
                                        muts_per_leaf, cfg['region_bias'], cfg['mutate_fr4'])

        groups_rows.append(f"{gname}\t{loc}\t{vall}\t{jall}\t{n_members}\t{v_len}\t{len(ancestor)}")

        for mi, observed in enumerate(observed_list):
            sid = f"{gname}_s{mi+1:04d}"          # без '|'/'*' — переживает IgBLAST/MSA
            fasta_lines.append(f">{sid}")
            fasta_lines.append(observed)
            seq_total += 1
            for m in codon_mutations(ancestor, observed, v_len, region_at, codon_at):
                gt_rows.append(
                    f"{sid}\t{gname}\t{vall}\t{jall}\t{m['region']}\t{m['imgt_codon']}\t"
                    f"{m['nt_pos']}\t{m['ref_codon']}\t{m['obs_codon']}\t"
                    f"{m['ref_aa']}\t{m['obs_aa']}\t{m['silent']}")
                region_counter[m['region']] = region_counter.get(m['region'], 0) + 1

    (out / 'synthetic_BCR.fasta').write_text('\n'.join(fasta_lines) + '\n')
    (out / 'ground_truth.tsv').write_text('\n'.join(gt_rows) + '\n')
    (out / 'groups.tsv').write_text('\n'.join(groups_rows) + '\n')

    # сводка + эхо конфига
    total_muts = len(gt_rows) - 1
    lines = []
    lines.append("=== СИНТЕТИЧЕСКИЙ ДАТАСЕТ BCR (из реальных germline) ===")
    lines.append(f"групп:                {cfg['groups']}")
    lines.append(f"последовательностей:  {seq_total}")
    lines.append(f"частота мутаций:      {cfg['mutation_rate']} % V-региона")
    lines.append(f"смещение по регионам: {cfg['region_bias']}")
    lines.append(f"форма клона:          {cfg['phylogeny']}")
    lines.append(f"локус:                {cfg['chain']}")
    lines.append(f"seed:                 {cfg['seed']}")
    lines.append("")
    lines.append(f"всего кодон-мутаций (истина): {total_muts}")
    lines.append("по регионам:")
    for reg in ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']:
        if reg in region_counter:
            lines.append(f"   {reg:5s} {region_counter[reg]}")
    lines.append("")
    lines.append("файлы: synthetic_BCR.fasta, ground_truth.tsv, groups.tsv")
    summary = '\n'.join(lines)
    (out / 'summary.txt').write_text(summary + '\n')
    print(summary)


if __name__ == '__main__':
    main()
