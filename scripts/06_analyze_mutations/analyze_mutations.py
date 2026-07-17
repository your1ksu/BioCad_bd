#!/usr/bin/env python3
"""
Analyze antibody sequences: find mutations relative to germline,
annotate by FR/CDR regions (IMGT numbering), output mutation table.

Uses IgBLAST for germline alignment and domain annotation.

Usage:
  python analyze_mutations.py -i sequences.fasta -o mutations.tsv [--ref-dir references/]
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# Standard genetic code
CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}


def translate(codon: str) -> str:
    """Translate a DNA codon to amino acid (single letter)."""
    codon = codon.upper().replace('U', 'T')
    return CODON_TABLE.get(codon, 'X')


def find_conda_base() -> Path | None:
    """Find miniconda/anaconda base directory."""
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        p = Path(conda_prefix)
        if p.parent.name == 'envs':
            return p.parent.parent
        return p.parent

    conda_root = os.environ.get('CONDA_ROOT')
    if conda_root:
        return Path(conda_root)

    home = Path.home()
    for candidate in [home / 'miniconda3', home / 'anaconda3',
                      home / 'miniconda', home / 'anaconda',
                      Path('/opt') / 'miniconda3', Path('/opt') / 'anaconda3']:
        if candidate.exists():
            return candidate
    return None


def find_igblast(conda_base: Path | None = None) -> str | None:
    """Find igblastn binary."""
    import shutil
    which = shutil.which('igblastn')
    if which:
        return which

    if not conda_base:
        conda_base = find_conda_base()
    if not conda_base:
        return None

    candidates = [
        conda_base / 'bin' / 'igblastn',
        conda_base / 'share' / 'igblast' / 'bin' / 'igblastn',
    ]
    for c in candidates:
        if c.exists():
            return str(c)
    return None


def run_igblast(query_fasta: str, ref_dir: str, out_file: str) -> str:
    """Run igblastn on query sequences."""
    igblast_bin = find_igblast(find_conda_base())
    if not igblast_bin:
        raise RuntimeError("igblastn not found. Install with: conda install -c bioconda igblast")

    v_db = os.path.join(ref_dir, 'all_V')
    d_db = os.path.join(ref_dir, 'all_D')
    j_db = os.path.join(ref_dir, 'all_J')
    aux_file = os.path.join(ref_dir, 'human_gl.aux')

    for db in [v_db, d_db, j_db]:
        if not os.path.exists(db + '.nsq'):
            raise RuntimeError(f"BLAST database not found: {db}. Run makeblastdb first.")

    cmd = [
        igblast_bin,
        '-germline_db_V', v_db,
        '-germline_db_D', d_db,
        '-germline_db_J', j_db,
        '-query', query_fasta,
        '-out', out_file,
        '-outfmt', '7 std qseq sseq',
        '-num_alignments_V', '1',
        '-num_alignments_D', '1',
        '-num_alignments_J', '1',
        '-organism', 'human',
        '-domain_system', 'imgt',
    ]
    if os.path.exists(aux_file):
        cmd.extend(['-auxiliary_data', aux_file])

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"igblastn failed:\n{result.stderr}")

    return out_file


def parse_igblast_output(igblast_file: str) -> list[dict]:
    """Parse IgBLAST -outfmt '7 std qseq sseq' output."""
    with open(igblast_file) as f:
        content = f.read()

    results = []
    query_blocks = re.split(r'^# Query:\s+', content, flags=re.MULTILINE)

    for block in query_blocks[1:]:
        lines = block.strip().split('\n')
        query_id = lines[0].strip()
        if not query_id:
            continue

        result = {
            'query_id': query_id,
            'v_gene': None, 'd_gene': None, 'j_gene': None,
            'chain_type': None,
            'domain_boundaries': {},
            'v_hit': None, 'j_hit': None,
        }

        current_section = 'header'
        for line in lines[1:]:
            stripped = line.strip()

            if stripped.startswith('#') or not stripped:
                if 'Alignment summary between query and top germline V gene hit' in stripped:
                    current_section = 'domain_summary'
                elif 'Hit table' in stripped and 'the first field indicates' in stripped:
                    current_section = 'hit_table'
                elif 'V-(D)-J junction details' in stripped:
                    current_section = 'junction'
                continue

            if current_section == 'header':
                if not result['v_gene']:
                    fields = stripped.split('\t')
                    if len(fields) >= 4:
                        v_part = fields[0].split('|')
                        d_part = fields[1].split('|')
                        j_part = fields[2].split('|')
                        result['v_gene'] = '|'.join(v_part[1:2]) if len(v_part) > 1 else fields[0]
                        result['d_gene'] = '|'.join(d_part[1:2]) if len(d_part) > 1 else fields[1]
                        result['j_gene'] = '|'.join(j_part[1:2]) if len(j_part) > 1 else fields[2]
                        result['chain_type'] = fields[3]
                        current_section = 'junction'
                continue

            elif current_section == 'junction':
                result['junction_details'] = stripped.split('\t')
                current_section = 'domain_summary'
                continue

            elif current_section == 'domain_summary':
                if stripped.startswith('Total') or stripped.startswith('CDR3'):
                    continue
                fields = stripped.split('\t')
                if len(fields) >= 8:
                    domain_name = fields[0]
                    qstart = int(fields[1])
                    qend = int(fields[2])
                    length = int(fields[3])
                    matches = int(fields[4])
                    mismatches = int(fields[5])
                    gaps = int(fields[6])
                    pct_id = float(fields[7]) if fields[7] not in ('N/A', '') else 0.0

                    if domain_name in ('FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT',
                                       'FR3-IMGT', 'CDR3-IMGT (germline)'):
                        result['domain_boundaries'][domain_name] = {
                            'qstart': qstart, 'qend': qend,
                            'length': length, 'matches': matches,
                            'mismatches': mismatches, 'gaps': gaps,
                            'pct_id': pct_id,
                        }
                continue

            elif current_section == 'hit_table':
                if stripped.startswith('V\t') or stripped.startswith('J\t'):
                    fields = stripped.split('\t')
                    if len(fields) >= 16:
                        hit_type = fields[0]
                        hit_data = {
                            'subject_id': fields[2],
                            'pident': float(fields[3]),
                            'align_length': int(fields[4]),
                            'mismatches': int(fields[5]),
                            'gap_opens': int(fields[6]),
                            'gaps': int(fields[7]),
                            'qstart': int(fields[8]),
                            'qend': int(fields[9]),
                            'sstart': int(fields[10]),
                            'send': int(fields[11]),
                            'evalue': float(fields[12]),
                            'bitscore': float(fields[13]),
                            'qseq': fields[14] if len(fields) > 14 else '',
                            'sseq': fields[15] if len(fields) > 15 else '',
                        }
                        if hit_type == 'V':
                            result['v_hit'] = hit_data
                        elif hit_type == 'J':
                            result['j_hit'] = hit_data
                continue

        results.append(result)

    return results


IMGT_DOMAIN_STARTS = {
    'FR1-IMGT': 1,
    'CDR1-IMGT': 27,
    'FR2-IMGT': 39,
    'CDR2-IMGT': 56,
    'FR3-IMGT': 66,
    'CDR3-IMGT (germline)': 105,
    'FR4-IMGT': 118,
}


def compute_mutations(qseq: str, sseq: str, qstart: int, domain_boundaries: dict) -> list[dict]:
    """Find amino acid mutations from aligned nucleotide sequences."""
    mutations = []
    alen = len(qseq)

    domain_map = {}
    for dname, bounds in domain_boundaries.items():
        for pos in range(bounds['qstart'], bounds['qend'] + 1):
            domain_map[pos] = dname

    prev_domain = None
    local_aa = 0

    i = 0
    aa_pos = 1
    while i < alen:
        q_positions = []
        s_positions = []
        has_indel = False

        for _ in range(3):
            while i < alen and qseq[i] == '-' and sseq[i] == '-':
                i += 1

            if i >= alen:
                break

            non_gap_up_to_i = sum(1 for j in range(i) if qseq[j] != '-')
            query_pos = qstart + non_gap_up_to_i

            if qseq[i] == '-' or sseq[i] == '-':
                has_indel = True
                if qseq[i] != '-':
                    q_positions.append((query_pos, qseq[i]))
                if sseq[i] != '-':
                    s_positions.append((query_pos, sseq[i]))
            else:
                q_positions.append((query_pos, qseq[i]))
                s_positions.append((query_pos, sseq[i]))

            i += 1

        if len(q_positions) < 3 or len(s_positions) < 3:
            break

        if has_indel:
            aa_pos += 1
            continue

        q_codon = ''.join(p[1] for p in q_positions)
        s_codon = ''.join(p[1] for p in s_positions)

        q_aa = translate(q_codon)
        s_aa = translate(s_codon)

        query_nt_pos = q_positions[0][0]
        domain = domain_map.get(query_nt_pos, 'Unknown')

        if domain != prev_domain:
            prev_domain = domain
            local_aa = 1
        else:
            local_aa += 1

        imgt_position = IMGT_DOMAIN_STARTS.get(domain, aa_pos) + local_aa - 1

        if q_aa != s_aa:
            is_silent = (q_aa == '*' and s_aa == '*') or (q_aa == s_aa)
            mutations.append({
                'aa_pos': aa_pos,
                'imgt_position': imgt_position,
                'domain': domain,
                'ref_aa': s_aa,
                'mut_aa': q_aa,
                'ref_codon': s_codon.lower(),
                'mut_codon': q_codon.lower(),
                'query_nt_pos': query_nt_pos,
                'is_silent': is_silent or (q_aa == s_aa),
            })

        aa_pos += 1

    return mutations


def write_mutations(results: list[dict], out_file: str) -> None:
    """Write mutation table as TSV."""
    fields = [
        'sequence_id', 'chain_type', 'v_gene', 'd_gene', 'j_gene',
        'region', 'aa_position', 'imgt_position',
        'ref_aa', 'mut_aa', 'ref_codon', 'mut_codon',
        'query_nt_pos', 'is_silent',
    ]

    with open(out_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(fields)

        for r in results:
            qid = r['query_id']
            chain = r['chain_type'] or ''
            vg = r['v_gene'] or ''
            dg = r['d_gene'] or ''
            jg = r['j_gene'] or ''

            mut_count = 0

            if r['v_hit']:
                v_domains = {k: v for k, v in r['domain_boundaries'].items()
                             if k in ('FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT',
                                      'FR3-IMGT', 'CDR3-IMGT (germline)')}
                v_muts = compute_mutations(r['v_hit']['qseq'], r['v_hit']['sseq'],
                                           r['v_hit']['qstart'], v_domains)
                for m in v_muts:
                    writer.writerow([
                        qid, chain, vg, dg, jg,
                        m['domain'], m['aa_pos'], m['imgt_position'],
                        m['ref_aa'], m['mut_aa'],
                        m['ref_codon'], m['mut_codon'],
                        m['query_nt_pos'], 'yes' if m['is_silent'] else 'no',
                    ])
                    mut_count += 1

            if r['j_hit']:
                j_domains = {'FR4-IMGT': {'qstart': r['j_hit']['qstart'], 'qend': r['j_hit']['qend']}}
                j_muts = compute_mutations(r['j_hit']['qseq'], r['j_hit']['sseq'],
                                           r['j_hit']['qstart'], j_domains)
                for m in j_muts:
                    writer.writerow([
                        qid, chain, vg, dg, jg,
                        'FR4-IMGT', m['aa_pos'], m['imgt_position'],
                        m['ref_aa'], m['mut_aa'],
                        m['ref_codon'], m['mut_codon'],
                        m['query_nt_pos'], 'yes' if m['is_silent'] else 'no',
                    ])
                    mut_count += 1

            if mut_count == 0:
                writer.writerow([qid, chain, vg, dg, jg, '', '', '', '', '', '', '', '', ''])


def write_summary(results: list[dict], out_file: str) -> None:
    """Write per-sequence summary TSV."""
    fields = [
        'sequence_id', 'chain_type', 'v_gene', 'd_gene', 'j_gene',
        'fr1_muts', 'cdr1_muts', 'fr2_muts', 'cdr2_muts', 'fr3_muts',
        'cdr3_muts', 'fr4_muts', 'total_muts',
    ]

    with open(out_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(fields)

        for r in results:
            qid = r['query_id']
            chain = r['chain_type'] or ''
            vg = r['v_gene'] or ''
            dg = r['d_gene'] or ''
            jg = r['j_gene'] or ''

            region_counts = {
                'FR1-IMGT': 0, 'CDR1-IMGT': 0, 'FR2-IMGT': 0,
                'CDR2-IMGT': 0, 'FR3-IMGT': 0, 'CDR3-IMGT': 0,
                'FR4-IMGT': 0,
            }
            total = 0

            if r['v_hit']:
                v_domains = {k: v for k, v in r['domain_boundaries'].items()
                             if k in region_counts}
                v_muts = compute_mutations(r['v_hit']['qseq'], r['v_hit']['sseq'],
                                           r['v_hit']['qstart'], v_domains)
                for m in v_muts:
                    if not m['is_silent'] and m['domain'] in region_counts:
                        region_counts[m['domain']] += 1
                        total += 1

            if r['j_hit']:
                j_domains = {'FR4-IMGT': {'qstart': r['j_hit']['qstart'], 'qend': r['j_hit']['qend']}}
                j_muts = compute_mutations(r['j_hit']['qseq'], r['j_hit']['sseq'],
                                           r['j_hit']['qstart'], j_domains)
                for m in j_muts:
                    if not m['is_silent']:
                        region_counts['FR4-IMGT'] += 1
                        total += 1

            writer.writerow([
                qid, chain, vg, dg, jg,
                region_counts['FR1-IMGT'], region_counts['CDR1-IMGT'],
                region_counts['FR2-IMGT'], region_counts['CDR2-IMGT'],
                region_counts['FR3-IMGT'], region_counts['CDR3-IMGT'],
                region_counts['FR4-IMGT'], total,
            ])


def main() -> None:
    parser = argparse.ArgumentParser(description='Analyze antibody mutations using IgBLAST')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('--ref-dir', default='references', help='Reference germline database directory')
    parser.add_argument('--format', choices=['mutations', 'summary'], default='mutations',
                        help='Output format (default: mutations)')
    args = parser.parse_args()

    ref_dir = args.ref_dir
    if not os.path.isdir(ref_dir):
        print(f"Error: reference directory not found: {ref_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Running IgBLAST on {args.input}...", file=sys.stderr)
    with tempfile.NamedTemporaryFile(suffix='_igblast.tsv', delete=False) as tmp:
        igblast_out = tmp.name

    try:
        run_igblast(args.input, ref_dir, igblast_out)

        print("Parsing IgBLAST output...", file=sys.stderr)
        results = parse_igblast_output(igblast_out)

        if not results:
            print("Warning: no results from IgBLAST", file=sys.stderr)
            open(args.output, 'w').close()
            sys.exit(0)

        print(f"Processed {len(results)} sequences", file=sys.stderr)

        if args.format == 'summary':
            write_summary(results, args.output)
            print(f"Summary written to {args.output}", file=sys.stderr)
        else:
            write_mutations(results, args.output)
            print(f"Mutation table written to {args.output}", file=sys.stderr)

        for r in results:
            mut_count = 0
            if r['v_hit']:
                v_domains = {k: v for k, v in r['domain_boundaries'].items()
                             if k in ('FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT',
                                      'FR3-IMGT', 'CDR3-IMGT (germline)')}
                v_muts = compute_mutations(r['v_hit']['qseq'], r['v_hit']['sseq'],
                                           r['v_hit']['qstart'], v_domains)
                mut_count = len([m for m in v_muts if not m['is_silent']])
            print(f"  {r['query_id']}: V={r['v_gene'] or '-'}, "
                  f"D={r['d_gene'] or '-'}, J={r['j_gene'] or '-'}, "
                  f"mutations={mut_count}", file=sys.stderr)

    finally:
        if os.path.exists(igblast_out):
            os.unlink(igblast_out)


if __name__ == '__main__':
    main()