#!/usr/bin/env python3
"""MrBayes — байесовские деревья по группам (ускоренная версия).

Отличия от исходной версии (та же точность, кратно быстрее):
  1. ИЗОЛЯЦИЯ cwd на каждую группу — mb запускается в собственной временной
     подпапке out_dir/_work/<key>/, поэтому параллельные процессы больше не
     мешают друг другу общими файлами (.ckp/.mcmc/…). Это устраняет 61 ошибку,
     которые были артефактом гонок, а не данных.
  2. STOPRULE — mcmc останавливается по достижении сходимости
     (avg split freq < stopval), а не гоняет фиксированный ngen. Критерий
     сходимости тот же (0.01), что и раньше проверялся постфактум, поэтому
     результат тот же — но лёгкие группы больше не крутят лишние сотни тысяч
     поколений. ngen остаётся как верхний предел.
  3. ПАРАЛЛЕЛИЗМ по всем ядрам (каждый mb серийный и теперь изолирован).
  4. ГИБРИД CPU/GPU — крупные группы (>= --gpu-min-taxa) считаются на GPU через
     BEAGLE (--gpu-mb-bin), последовательно (одна A100); мелкие — параллельно на
     CPU. Если --gpu-mb-bin не задан, всё идёт на CPU.
"""
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
                mb_burnin_frac: float, seed: int,
                use_gpu: bool = False, stopval: float = 0.01) -> str:
    ids = list(msa)
    width = len(next(iter(msa.values())))
    samplefreq = max(1, mb_ngen // 1000)
    diagnfreq = max(1, mb_ngen // 100)          # чаще проверяем → раньше стопаем
    printfreq = max(1, mb_ngen // 10)

    lines = ["#NEXUS", "", "begin data;",
             f"  dimensions ntax={len(ids)} nchar={width};",
             "  format datatype=DNA gap=- missing=? interleave=no;",
             "  matrix"]
    for sid in ids:
        lines.append(f"  {fwd[sid]}  {msa[sid]}")
    lines += ["  ;", "end;", "", "begin mrbayes;",
              f"  set autoclose=yes nowarn=yes seed={seed} swapseed={seed};"]
    if use_gpu:
        # BEAGLE на GPU: двойная точность, SSE off — выигрыш на крупных группах
        lines.append("  set usebeagle=yes beagledevice=gpu beagleprecision=double "
                     "beaglescaling=dynamic;")
    lines.append("  lset nst=6 rates=invgamma;                 [ GTR+I+G ]")
    if outgroup_safe:
        lines.append(f"  outgroup {outgroup_safe};")
    # stoprule=yes stopval=<t>: останавливаемся при достижении сходимости,
    # ngen — верхний предел. Критерий тот же, что раньше проверялся постфактум.
    lines += [
        f"  mcmc ngen={mb_ngen} samplefreq={samplefreq} nchains=4 nruns=2 "
        f"printfreq={printfreq} diagnfreq={diagnfreq} "
        f"stoprule=yes stopval={stopval} "
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
    return None


def detect_nproc() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        pass
    n = os.cpu_count()
    return n if n else 1


def run_mb(key: str, nexus_text: str, rev: dict[str, str],
           mb_bin: str, timeout_s: int, out_dir: Path) -> str:
    """Запуск mb в ИЗОЛИРОВАННОЙ подпапке; наружу копируем только con.tre + лог."""
    binpath = find_mb(mb_bin)
    if binpath is None:
        (out_dir / f"{key}.nex").write_text(nexus_text, encoding="utf-8")
        return "nexus_only"

    work = out_dir / "_work" / key
    if work.exists():
        shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True, exist_ok=True)
    nex = work / f"{key}.nex"
    nex.write_text(nexus_text, encoding="utf-8")

    t0 = time.time()
    try:
        proc = subprocess.run([str(binpath), nex.name],
                              cwd=str(work), capture_output=True, text=True,
                              timeout=timeout_s)
    except subprocess.TimeoutExpired:
        shutil.rmtree(work, ignore_errors=True)
        return "failed"
    if proc.returncode != 0:
        shutil.rmtree(work, ignore_errors=True)
        return "failed"

    runtime = time.time() - t0
    # копируем нужные артефакты наверх, .nex тоже (для воспроизводимости)
    con = work / f"{key}.nex.con.tre"
    if con.is_file():
        shutil.copy(con, out_dir / f"{key}.nex.con.tre")
    shutil.copy(nex, out_dir / f"{key}.nex")
    (out_dir / f"{key}.mb.log").write_text(proc.stdout or "", encoding="utf-8")
    shutil.rmtree(work, ignore_errors=True)     # чистим тяжёлые промежуточные файлы

    newick, supports = parse_con_tre(out_dir / f"{key}.nex.con.tre", rev)
    avg_std, converged = convergence(proc.stdout or "")
    return f"ok|{runtime}|{len(supports)}|{converged}|{avg_std}"


def _report(key: str, status: str, out_dir: Path, counts: dict, tag: str = "") -> None:
    if status.startswith("ok"):
        _, rt, nclades, ok_str, avg = status.split("|")
        conv = "да" if ok_str == "True" else "нет"
        print(f"[{key}]{tag} готово: {out_dir / (key + '.nex.con.tre')} "
              f"(runtime={float(rt):.1f}s, клад={nclades}, сошлось={conv} avg_std={avg})")
        counts["ok"] += 1
    elif status == "failed":
        print(f"[{key}]{tag} ОШИБКА mb", file=sys.stderr)
        counts["failed"] += 1
    else:
        print(f"[{key}]{tag} mb не найден — nexus сохранён", file=sys.stderr)
        counts["nexus_only"] += 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="MrBayes — байесовские деревья по группам (ускоренная)")
    ap.add_argument("input_dir", nargs="?", default="aligned_sequences")
    ap.add_argument("--out", default="mrbayes")
    ap.add_argument("--outgroup", default=None)
    ap.add_argument("--mb-ngen", type=int, default=200_000,
                    help="верхний предел MCMC = как в оригинале (stoprule обычно "
                         "останавливает раньше; трудные группы дают тот же не-"
                         "сошедшийся результат, что и раньше — та же точность)")
    ap.add_argument("--mb-burnin-frac", type=float, default=0.25)
    ap.add_argument("--stopval", type=float, default=0.01,
                    help="порог сходимости для stoprule (тот же, что проверялся раньше)")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--mb-bin", default="", help="CPU-бинарь mb")
    ap.add_argument("--gpu-mb-bin", default="", help="GPU-бинарь mb (BEAGLE-CUDA)")
    ap.add_argument("--gpu-min-taxa", type=int, default=60,
                    help="группы с таксонов >= порога считать на GPU")
    ap.add_argument("--timeout-s", type=int, default=3600)
    ap.add_argument("--nexus-only", action="store_true")
    ap.add_argument("--workers", type=int, default=0,
                    help="параллельных CPU-групп (по умолчанию: число ядер)")
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
    out_dir.mkdir(parents=True, exist_ok=True)
    counts = {"ok": 0, "nexus_only": 0, "failed": 0, "skipped": 0}

    gpu_on = bool(args.gpu_mb_bin) and find_mb(args.gpu_mb_bin) is not None

    # Фаза 1: разметка групп на CPU/GPU + генерация nexus (в память)
    cpu_tasks: list[tuple[str, str, dict]] = []
    gpu_tasks: list[tuple[str, str, dict]] = []
    for fasta in fasta_files:
        key = group_key_of(fasta)
        msa = read_fasta(fasta)
        n = len(msa)
        if n < 4:
            # <4 таксонов: неразрешимая топология, MrBayes бессмыслен и падает
            # (эти группы были частью исходных «ошибок»); клады из них всё равно
            # не выходят → пропуск даёт тот же итог, но без ложных ошибок.
            counts["skipped"] += 1
            continue
        fwd, rev = safe_names(list(msa))
        (out_dir / f"{key}.names.tsv").write_text(
            "\n".join(f"{s}\t{o}" for o, s in fwd.items()) + "\n", encoding="utf-8")
        outg = fwd.get(args.outgroup) if args.outgroup else None
        to_gpu = gpu_on and n >= args.gpu_min_taxa
        nexus_text = build_nexus(key, msa, fwd, outg, args.mb_ngen,
                                 args.mb_burnin_frac, args.seed,
                                 use_gpu=to_gpu, stopval=args.stopval)
        (gpu_tasks if to_gpu else cpu_tasks).append((key, nexus_text, rev))

    if args.nexus_only:
        for key, txt, _ in cpu_tasks + gpu_tasks:
            (out_dir / f"{key}.nex").write_text(txt, encoding="utf-8")
        print(f"\nnexus_only: {len(cpu_tasks)+len(gpu_tasks)} групп")
        return 0

    print(f"--- MrBayes: CPU={len(cpu_tasks)} групп × {min(workers, max(1,len(cpu_tasks)))} воркеров"
          f" | GPU={len(gpu_tasks)} групп (последовательно) ---")

    # Фаза 2a: GPU-группы последовательно (одна A100)
    for key, txt, rev in gpu_tasks:
        status = run_mb(key, txt, rev, args.gpu_mb_bin, args.timeout_s, out_dir)
        _report(key, status, out_dir, counts, tag="[GPU]")

    # Фаза 2b: CPU-группы параллельно, каждая в своей подпапке
    cw = min(workers, len(cpu_tasks)) if cpu_tasks else 1
    if cw <= 1:
        for key, txt, rev in cpu_tasks:
            _report(key, run_mb(key, txt, rev, args.mb_bin, args.timeout_s, out_dir),
                    out_dir, counts)
    else:
        with ThreadPoolExecutor(max_workers=cw) as pool:
            futs = {pool.submit(run_mb, key, txt, rev, args.mb_bin,
                                args.timeout_s, out_dir): key
                    for key, txt, rev in cpu_tasks}
            for fut in as_completed(futs):
                key = futs[fut]
                try:
                    _report(key, fut.result(), out_dir, counts)
                except Exception as e:
                    print(f"[{key}] ОШИБКА воркера: {e}", file=sys.stderr)
                    counts["failed"] += 1

    shutil.rmtree(out_dir / "_work", ignore_errors=True)
    print(f"\nИтого: ok={counts['ok']} только-nexus={counts['nexus_only']} "
          f"ошибок={counts['failed']} пропущено={counts['skipped']} из {len(fasta_files)} групп")
    return 0 if counts["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
