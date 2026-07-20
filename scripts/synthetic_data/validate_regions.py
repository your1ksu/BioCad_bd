#!/usr/bin/env python3
"""
Проверка: правильно ли пайплайн определил РЕГИОНЫ мутаций.

Сравнивает истину конструктора (ground_truth.tsv) с разметкой пайплайна
(06_analyze_mutations → mutation_tables/*/mutations.tsv). Даёт:
  * position→region: совпадает ли reported region с region_of_codon(imgt_position)
    (независимо от синтетики — чистая корректность IMGT-разметки);
  * матрицу ошибок «истинный регион × определённый регион»;
  * recall по регионам (сколько заложенных мутаций найдено в верном регионе);
  * разбор «Unknown» и стыка CDR3.

Запуск:
  python validate_regions.py --ground-truth OUT/ground_truth.tsv \
      --mutations RESULT/mutation_tables
"""

import argparse
import csv
import glob
import os
from collections import Counter, defaultdict


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


def norm_region(r: str) -> str:
    """'FR1-IMGT' / 'CDR3-IMGT (germline)' / 'FR4-IMGT' -> 'FR1'/'CDR3'/'FR4'."""
    r = (r or '').strip()
    if r.startswith('CDR3'):
        return 'CDR3'
    return r.split('-')[0] if r and r != 'Unknown' else r


REGION_ORDER = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4', 'Unknown', '']


def gen_allele(sid):
    """Аллель, из которой синтезировали: G05_IGHV4-61_05_IGHJ6_02_s0041 -> IGHV4-61*05."""
    import re
    m = re.match(r'G\d+_(IG[HKL]V[\w-]+?)_(\d+)_IG', sid)
    return m.group(1) + '*' + m.group(2) if m else None


def read_ground_truth(path):
    """(sequence_id, imgt_codon) -> dict; плюс counts по регионам."""
    gt = {}
    counts = Counter()
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = (row['sequence_id'], int(row['imgt_codon']))
            gt[key] = {'region': row['region'], 'silent': row['silent'],
                       'ref_aa': row['ref_aa'], 'mut_aa': row['obs_aa']}
            counts[row['region']] += 1
    return gt, counts


def read_pipeline_mutations(mut_dir):
    """Все строки mutations.tsv из mutation_tables/*/mutations.tsv."""
    rows = []
    for path in glob.glob(os.path.join(mut_dir, '*', 'mutations.tsv')):
        with open(path) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if not row.get('region') or not row.get('imgt_position'):
                    continue
                try:
                    imgt = int(row['imgt_position'])
                except ValueError:
                    continue
                rows.append({
                    'sequence_id': row.get('sequence_id', ''),
                    'imgt': imgt,
                    'region': norm_region(row['region']),
                    'silent': row.get('is_silent', ''),
                    'ref_aa': row.get('ref_aa', ''),
                    'mut_aa': row.get('mut_aa', ''),
                    'v_gene': row.get('v_gene', ''),
                })
    return rows


def bar(title):
    print('\n' + '=' * 62)
    print(title)
    print('=' * 62)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ground-truth', required=True)
    ap.add_argument('--mutations', required=True,
                    help='каталог mutation_tables (с подпапками */mutations.tsv)')
    args = ap.parse_args()

    gt, gt_counts = read_ground_truth(args.ground_truth)
    muts = read_pipeline_mutations(args.mutations)

    # одна и та же последовательность попадает в НЕСКОЛЬКО вложенных клад и
    # анализируется в каждой → строки дублируются. Схлопываем до различных
    # (sequence_id, imgt) для честного recall/precision.
    dedup = {}
    for m in muts:
        dedup[(m['sequence_id'], m['imgt'])] = m
    muts_unique = list(dedup.values())

    bar('ВХОД')
    print(f"истинных кодон-мутаций (ground truth):        {len(gt)}")
    print(f"мутаций в разметке пайплайна (все строки):    {len(muts)}")
    print(f"уникальных (sequence_id, imgt) в разметке:    {len(muts_unique)}")
    print("  (строк больше из-за вложенных клад: одна и та же последовательность")
    print("   анализируется в каждой кладе, куда входит)")

    # --- 1. position -> region consistency -----------------------------------
    consistent = 0
    checkable = 0
    for m in muts_unique:
        if not m['region'] or m['region'] == 'Unknown':
            continue
        checkable += 1
        if m['region'] == region_of_codon(m['imgt']):
            consistent += 1
    bar('1) СОГЛАСОВАННОСТЬ imgt_position → region (IMGT-корректность)')
    if checkable:
        print(f"согласовано: {consistent}/{checkable} = {100*consistent/checkable:.1f}%")
    else:
        print("нет мутаций с определённым регионом")

    # --- 2. матрица ошибок против истины конструктора ------------------------
    confusion = defaultdict(Counter)   # true_region -> Counter(detected_region)
    matched = 0
    unmatched = 0
    for m in muts_unique:
        key = (m['sequence_id'], m['imgt'])
        if key in gt:
            true_region = gt[key]['region']
            matched += 1
        else:
            # позиция не заложена нами (стык CDR3, FR4 или неточность позиции):
            # берём эталон по IMGT-номеру, чтобы всё равно оценить корректность
            true_region = region_of_codon(m['imgt'])
            unmatched += 1
        confusion[true_region][m['region'] or 'Unknown'] += 1

    detected_regions = REGION_ORDER
    bar('2) МАТРИЦА «истинный регион × определённый пайплайном»')
    print(f"совпало по (sequence_id, imgt): {matched}   прочих: {unmatched}")
    header = 'true \\ det'.ljust(11) + ''.join(d.rjust(8) for d in detected_regions if d)
    print(header)
    correct_total = 0
    grand_total = 0
    for tr in ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']:
        if tr not in confusion:
            continue
        line = tr.ljust(11)
        for d in detected_regions:
            if not d:
                continue
            c = confusion[tr].get(d, 0)
            line += str(c).rjust(8)
            grand_total += c
            if d == tr:
                correct_total += c
        print(line)
    if grand_total:
        print(f"\nобщая точность региона: {correct_total}/{grand_total} = "
              f"{100*correct_total/grand_total:.1f}%")

    # --- 3. recall/precision по замене а.к. (устойчиво к нумерации imgt) ------
    # imgt_position пайплайна дрейфует от истинного IMGT после scaffold-гэпов,
    # поэтому сопоставляем НЕ по номеру, а по факту замены: (регион, ref>mut) в
    # пределах одной последовательности. Ограничиваемся той же аллелью, чтобы
    # убрать шум переназначения аллели IgBLAST (см. диагностику).
    bar('3) RECALL / PRECISION ПО ЗАМЕНЕ А.К. (та же аллель, устойчиво к imgt)')
    detected_by_seq = defaultdict(Counter)   # sid -> Counter((region, ref>mut))
    for m in muts_unique:
        detected_by_seq[m['sequence_id']][(m['region'], m['ref_aa'] + '>' + m['mut_aa'])] += 1
    # какая аллель у последовательности в разметке пайплайна
    pipe_allele = {}
    for m in muts_unique:
        pipe_allele.setdefault(m['sequence_id'], m['v_gene'])

    planted_by_seq = defaultdict(Counter)
    for (sid, _imgt), g in gt.items():
        planted_by_seq[sid][(g['region'], g['ref_aa'] + '>' + g['mut_aa'])] += 1

    same = {sid for sid in detected_by_seq
            if pipe_allele.get(sid) == gen_allele(sid)}
    recall_hit = Counter(); planted_c = Counter(); detected_c = Counter()
    for sid in same:
        p = planted_by_seq[sid]; d = detected_by_seq[sid]
        for key, pc in p.items():
            reg = key[0]
            planted_c[reg] += pc
            recall_hit[reg] += min(pc, d.get(key, 0))
        for key, dc in d.items():
            detected_c[key[0]] += dc
    print(f"последовательностей с той же аллелью: {len(same)}")
    print('регион'.ljust(8) + 'заложено'.rjust(9) + 'найдено'.rjust(9) + 'recall'.rjust(9))
    for reg in ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3']:
        p = planted_c.get(reg, 0); h = recall_hit.get(reg, 0)
        rec = f"{100*h/p:.1f}%" if p else "—"
        print(reg.ljust(8) + str(p).rjust(9) + str(recall_hit.get(reg, 0)).rjust(9) + rec.rjust(9))
    tp = sum(planted_c.values()); th = sum(recall_hit.values())
    if tp:
        print(f"\nИТОГ recall (заложенные замены, найденные в верном регионе): "
              f"{th}/{tp} = {100*th/tp:.1f}%")


if __name__ == '__main__':
    main()
