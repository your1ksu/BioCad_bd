#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
from pathlib import Path

FASTA_EXTS = (".fa", ".fasta", ".fas", ".aln")
_SPLIT_RE = re.compile(r"Average standard deviation of split frequencies:\s*([0-9.]+)")


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].strip().split()[0]
            seqs[cur] = ""
        elif cur is not None:
            seqs[cur] += line.strip()
    return seqs


def group_key_of(fasta_path: Path) -> str:
    name = fasta_path.stem
    if name.endswith("_aligned"):
        name = name[: -len("_aligned")]
    return name


def safe_names(ids: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    fwd = {sid: f"T{i:04d}" for i, sid in enumerate(ids)}
    rev = {v: k for k, v in fwd.items()}
    return fwd, rev


def build_nexus(key: str, msa: dict[str, str], fwd: dict[str, str],
                outgroup_safe: str | None, mb_ngen: int,
                mb_burnin_frac: float, seed: int) -> str:
    ids = list(msa)
    width = len(next(iter(msa.values())))
    samplefreq = max(1, mb_ngen // 1000)
    diagnfreq = max(1, mb_ngen // 10)
    printfreq = max(1, mb_ngen // 10)

    lines = ["#NEXUS", "", "begin data;",
             f"  dimensions ntax={len(ids)} nchar={width};",
             "  format datatype=DNA gap=- missing=? interleave=no;",
             "  matrix"]
    for sid in ids:
        lines.append(f"  {fwd[sid]}  {msa[sid]}")
    lines += ["  ;", "end;", "", "begin mrbayes;",
              f"  set autoclose=yes nowarn=yes seed={seed} swapseed={seed};",
              "  lset nst=6 rates=invgamma;                 [ GTR+I+G ]"]
    if outgroup_safe:
        lines.append(f"  outgroup {outgroup_safe};")
    lines += [
        f"  mcmc ngen={mb_ngen} samplefreq={samplefreq} nchains=4 nruns=2 "
        f"printfreq={printfreq} diagnfreq={diagnfreq} "
        f"starttree=random savebrlens=yes;",
        f"  sumt burninfrac={mb_burnin_frac} conformat=simple;",
        f"  sump burninfrac={mb_burnin_frac};",
        "end;", ""]
    return "\n".join(lines)


def parse_con_tre(path: str | Path, rev: dict[str, str]) -> tuple[str, dict[str, dict]]:
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
            pp = pp / 100.0 if pp > 1.0 else pp
            supports["|".join(sorted(leaves))] = {"posterior": round(pp, 4)}
    out = StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue().strip(), supports


def convergence(stdout: str) -> tuple[float | None, bool]:
    vals = _SPLIT_RE.findall(stdout or "")
    if not vals:
        return None, False
    last = float(vals[-1])
    return last, last < 0.01


def find_mb(explicit: str = "") -> Path | None:
    if explicit:
        p = Path(explicit)
        if p.is_file() and p.stat().st_mode & 0o111:
            return p
        return None
    for cand in ("mb", "mrbayes"):
        which = shutil.which(cand)
        if which:
            return Path(which)
    for envbin in [
        Path("/opt/miniconda3/envs/bv2026_msa/bin"),
        Path("/opt/miniconda3/envs/bv2026/bin"),
    ]:
        for cand in ("mb", "mrbayes"):
            p = envbin / cand
            if p.is_file() and p.stat().st_mode & 0o111:
                return p
    return None


def detect_nproc() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        pass
    n = os.cpu_count()
    return n if n else 1


def write_nexus(key: str, msa: dict[str, str], outgroup_id: str | None,
                mb_ngen: int, mb_burnin_frac: float, seed: int,
                out_dir: Path) -> tuple[Path, dict[str, str]]:
    fwd, rev = safe_names(list(msa))
    names_tsv = out_dir / f"{key}.names.tsv"
    names_tsv.write_text(
        "\n".join(f"{safe}\t{orig}" for orig, safe in fwd.items()) + "\n",
        encoding="utf-8")
    outg = fwd.get(outgroup_id) if outgroup_id else None
    nexus_text = build_nexus(key, msa, fwd, outg, mb_ngen, mb_burnin_frac, seed)
    nex_path = out_dir / f"{key}.nex"
    nex_path.write_text(nexus_text, encoding="utf-8")
    return nex_path, rev


def run_mb(key: str, nex_path: Path, rev: dict[str, str],
           mb_bin: str, timeout_s: int, out_dir: Path) -> str:
    binpath = find_mb(mb_bin)
    if binpath is None:
        return "nexus_only"

    t0 = time.time()
    try:
        proc = subprocess.run([str(binpath), nex_path.name],
                              cwd=str(out_dir), capture_output=True, text=True,
                              timeout=timeout_s)
    except subprocess.TimeoutExpired:
        return "failed"
    if proc.returncode != 0:
        return "failed"

    (out_dir / f"{key}.mb.log").write_text(proc.stdout or "", encoding="utf-8")
    runtime = time.time() - t0

    newick, supports = parse_con_tre(out_dir / f"{key}.nex.con.tre", rev)
    avg_std, converged = convergence(proc.stdout or "")
    return f"ok|{runtime}|{len(supports)}|{converged}|{avg_std}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="MrBayes — байесовские деревья по группам")
    ap.add_argument("input_dir", nargs="?", default="aligned_sequences",
                    help="папка с *.fasta (одно выравнивание = одна группа)")
    ap.add_argument("--out", default="mrbayes", help="выходная папка")
    ap.add_argument("--outgroup", default=None,
                    help="id таксона-outgroup (если есть в fasta)")
    ap.add_argument("--mb-ngen", type=int, default=200_000, help="длина цепи MCMC")
    ap.add_argument("--mb-burnin-frac", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--mb-bin", default="", help="явный путь к бинарю mb")
    ap.add_argument("--timeout-s", type=int, default=3600)
    ap.add_argument("--nexus-only", action="store_true",
                    help="только фаза 1 (сгенерировать .nex)")
    ap.add_argument("--workers", type=int, default=0,
                    help="сколько групп запускать параллельно (по умолчанию: число ядер)")
    ap.add_argument("--max-workers", type=int, default=0,
                    help="максимум параллельных mb (чтобы не убить диск/RAM)")
    args = ap.parse_args(argv)

    in_dir = Path(args.input_dir)
    out_dir = Path(args.out)
    if not in_dir.is_dir():
        print(f"Ошибка: папка не найдена: {in_dir}", file=sys.stderr)
        return 1

    fasta_files = sorted(p for p in in_dir.iterdir() if p.suffix.lower() in FASTA_EXTS)
    if not fasta_files:
        print(f"В {in_dir} не найдено *.fa/*.fasta/*.fas/*.aln", file=sys.stderr)
        return 0

    nproc = detect_nproc()
    workers = args.workers if args.workers > 0 else nproc
    if args.max_workers > 0:
        workers = min(workers, args.max_workers)
    else:
        workers = min(workers, 4)

    out_dir.mkdir(parents=True, exist_ok=True)
    counts = {"ok": 0, "nexus_only": 0, "failed": 0, "skipped": 0}

    # Фаза 1: генерация nexus (последовательно, быстро)
    tasks: list[tuple[str, Path, dict[str, str]]] = []
    for fasta in fasta_files:
        key = group_key_of(fasta)
        msa = read_fasta(fasta)
        if len(msa) < 3:
            print(f"[{key}] пропущено: меньше 3 последовательностей ({len(msa)})")
            counts["skipped"] += 1
            continue
        nex_path, rev = write_nexus(key, msa, args.outgroup,
                                     args.mb_ngen, args.mb_burnin_frac,
                                     args.seed, out_dir)
        print(f"[{key}] nexus готов: {nex_path}")
        tasks.append((key, nex_path, rev))

    if args.nexus_only or not tasks:
        print(f"\nИтого: nexus_only={len(tasks)} пропущено={counts['skipped']}")
        return 0

    # Фаза 2: запуск mb (параллельно)
    workers = min(workers, len(tasks))
    print(f"\n--- MrBayes: {len(tasks)} групп, {workers} воркеров ---")

    if workers <= 1:
        for key, nex_path, rev in tasks:
            status = run_mb(key, nex_path, rev, args.mb_bin, args.timeout_s, out_dir)
            if status.startswith("ok"):
                _, rt, nclades, ok_str, avg = status.split("|")
                converged_str = "да" if ok_str == "True" else "нет"
                print(f"[{key}] готово: {out_dir / (key + '.nex.con.tre')} "
                      f"(runtime={float(rt):.1f}s, клад={nclades}, "
                      f"сошлось={converged_str} avg_std={avg})")
                counts["ok"] += 1
            elif status == "failed":
                print(f"[{key}] ОШИБКА mb", file=sys.stderr)
                counts["failed"] += 1
            else:
                print(f"[{key}] mb не найден — nexus сохранён", file=sys.stderr)
                counts["nexus_only"] += 1
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            fut_to_key = {
                pool.submit(run_mb, key, nex_path, rev,
                            args.mb_bin, args.timeout_s, out_dir): key
                for key, nex_path, rev in tasks
            }
            for future in as_completed(fut_to_key):
                key = fut_to_key[future]
                try:
                    status = future.result()
                except Exception as e:
                    print(f"[{key}] ОШИБКА воркера: {e}", file=sys.stderr)
                    counts["failed"] += 1
                    continue

                if status.startswith("ok"):
                    _, rt, nclades, ok_str, avg = status.split("|")
                    converged_str = "да" if ok_str == "True" else "нет"
                    print(f"[{key}] готово: {out_dir / (key + '.nex.con.tre')} "
                          f"(runtime={float(rt):.1f}s, клад={nclades}, "
                          f"сошлось={converged_str} avg_std={avg})")
                    counts["ok"] += 1
                elif status == "failed":
                    print(f"[{key}] ОШИБКА mb", file=sys.stderr)
                    counts["failed"] += 1
                else:
                    print(f"[{key}] mb не найден — nexus сохранён", file=sys.stderr)
                    counts["nexus_only"] += 1

    print(f"\nИтого: ok={counts['ok']} только-nexus={counts['nexus_only']} "
          f"ошибок={counts['failed']} пропущено={counts['skipped']} из {len(fasta_files)} групп")
    return 0 if counts["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
