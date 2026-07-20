#!/usr/bin/env python3
"""
Master pipeline runner for BCR analysis.

Usage:
    python scripts/run_pipeline.py -k batch1
    python scripts/run_pipeline.py -i data/batch1/BCR_data.tsv -o results/report_17072026_143000
"""

import argparse
import json
import os
import shutil
import sys
import subprocess
import time
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "results"
REFERENCES_DIR = DATA_DIR / "references"
ENV_NAME = "biocad_bcr_pipeline_environment"
ENV_YML = REPO_ROOT / "environment.yml"

PYTHON_BIN: str | None = None
BASH_BIN: str | None = None
CONDA_PREFIX: str | None = None

sys.path.insert(0, str(SCRIPTS_DIR))
from shared.config import PipelineConfig, DEFAULT_CONFIG


def ensure_env():
    global CONDA_PREFIX

    try:
        result = subprocess.run(
            ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
        )
        envs = json.loads(result.stdout).get("envs", [])
    except Exception:
        print("Error: conda not found. Install miniconda first.")
        sys.exit(1)

    env_paths = [p for p in envs if p.endswith("/" + ENV_NAME)]
    if not env_paths:
        print(f"Creating conda environment '{ENV_NAME}' from {ENV_YML} ...")
        subprocess.run(["conda", "env", "create", "-f", str(ENV_YML)], check=True)
        env_paths = [p for p in envs if p.endswith("/" + ENV_NAME)]
        print("Environment created.")

    CONDA_PREFIX = env_paths[0]

    global PYTHON_BIN, BASH_BIN
    PYTHON_BIN = subprocess.run(
        ["conda", "run", "-n", ENV_NAME, "which", "python"],
        capture_output=True, text=True, check=True
    ).stdout.strip()
    BASH_BIN = subprocess.run(
        ["conda", "run", "-n", ENV_NAME, "which", "bash"],
        capture_output=True, text=True, check=True
    ).stdout.strip()


def conda_env():
    env = os.environ.copy()
    env["PATH"] = f"{CONDA_PREFIX}/bin:{env.get('PATH', '')}"
    return env


ALL_STEPS = ["filter", "group", "filter_groups", "msa", "iqtree", "mrbayes", "viz", "clades", "mutations"]


def run_cmd(cmd, cwd=None, log_file=None, description=None):
    if description:
        print(f"\n{'='*60}")
        print(f"STEP: {description}")
        print(f"CMD:  {' '.join(cmd)}")
        print(f"{'='*60}")

    start = time.time()
    env = conda_env()
    try:
        if log_file:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            with open(log_file, "w") as lf:
                lf.write(f"CMD: {' '.join(cmd)}\n")
                lf.write(f"CWD: {cwd or REPO_ROOT}\n")
                lf.write(f"START: {datetime.now()}\n\n")
                result = subprocess.run(
                    cmd, cwd=cwd or REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
                    text=True, env=env,
                )
        else:
            result = subprocess.run(
                cmd, cwd=cwd or REPO_ROOT, capture_output=True, text=True, env=env,
            )
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(result.stderr, file=sys.stderr)

        elapsed = time.time() - start
        if result.returncode == 0:
            print(f"\u2713 DONE in {elapsed:.1f}s")
            return True
        else:
            print(f"\u2717 FAILED (exit {result.returncode}) in {elapsed:.1f}s")
            if log_file:
                print(f"  See log: {log_file}")
            return False
    except subprocess.TimeoutExpired:
        print(f"\u2717 TIMEOUT after {time.time() - start:.1f}s")
        return False
    except Exception as e:
        print(f"\u2717 ERROR: {e}")
        return False


def step1_filter(input_tsv, output_tsv, log_dir):
    return run_cmd([
        PYTHON_BIN, str(SCRIPTS_DIR / "01_filter_sequences" / "filter_sequences.py"),
        "-i", str(input_tsv), "-o", str(output_tsv), "-r", str(REFERENCES_DIR),
    ], log_file=log_dir / "step1_filter.log", description="Step 1: Filter BCR data")


def step2_group(filtered_tsv, output_dir, log_dir, grouping_strategy):
    return run_cmd([
        PYTHON_BIN, str(SCRIPTS_DIR / "02_group_by_germlines" / "group_by_germlines.py"),
        "-i", str(filtered_tsv), "-o", str(output_dir), "-r", str(REFERENCES_DIR),
        "--grouping-strategy", grouping_strategy,
    ], log_file=log_dir / "step2_group.log", description="Step 2: Group by germlines")


def step2b_filter_groups(vj_dir, filtered_dir, log_dir, min_size, max_size):
    return run_cmd([
        PYTHON_BIN, str(SCRIPTS_DIR / "filter_by_symbol_count" / "filter_by_symbol_count.py"),
        "-i", str(vj_dir), "-o", str(filtered_dir),
        "--min", str(min_size), "--max", str(max_size),
    ], log_file=log_dir / "step2b_filter_groups.log", description="Step 2b: Filter groups by size")


def step3_msa(grouped_vj_dir, output_dir, log_dir):
    return run_cmd([
        PYTHON_BIN, str(SCRIPTS_DIR / "03_multiple_alignment" / "multiple_alignment.py"),
        "-i", str(grouped_vj_dir), "-o", str(output_dir),
    ], log_file=log_dir / "step3_msa.log", description="Step 3: MSA (MAFFT)")


def step4a_iqtree(aligned_dir, trees_dir, log_dir, iqtree_model):
    return run_cmd([
        PYTHON_BIN, str(SCRIPTS_DIR / "04a_build_trees_iqtree" / "build_trees_iqtree.py"),
        "-i", str(aligned_dir), "-o", str(trees_dir),
        "--model", iqtree_model,
    ], log_file=log_dir / "step4a_iqtree.log", description="Step 4a: IQ-TREE (ML)")


def step4b_mrbayes(aligned_dir, mrbayes_dir, log_dir, config):
    return run_cmd([
        PYTHON_BIN, str(SCRIPTS_DIR / "04b_build_trees_mrbayes" / "build_trees_mrbayes.py"),
        str(aligned_dir), "--out", str(mrbayes_dir),
        "--mb-ngen", str(config.mb_ngen_default),
    ], log_file=log_dir / "step4b_mrbayes.log", description="Step 4b: MrBayes (Bayesian)")


def step5_viz(trees_dir, viz_dir, log_dir):
    return run_cmd([
        PYTHON_BIN, str(SCRIPTS_DIR / "visualize_trees" / "visualize_trees.py"),
        "-i", str(trees_dir), "-o", str(viz_dir),
    ], log_file=log_dir / "step5_viz.log", description="Step 5: Tree visualization")


def step6_clades(trees_dir, mrbayes_dir, groups_dir, aligned_dir, log_dir):
    clades_dir = groups_dir / "clades"
    return run_cmd([
        PYTHON_BIN, str(SCRIPTS_DIR / "05_clade_search" / "clade_search.py"),
        "--iqtree-dir", str(trees_dir), "--mrbayes-dir", str(mrbayes_dir),
        "--out", str(groups_dir / "report.json"),
        "--clades-fasta-dir", str(clades_dir), "--aligned-dir", str(aligned_dir),
    ], log_file=log_dir / "step6_clades.log", description="Step 6: Confident clades report")


def step7_mutations(clades_dir, mutations_dir, ref_dir, log_dir):
    return run_cmd([
        PYTHON_BIN, str(SCRIPTS_DIR / "06_analyze_mutations" / "run_mutations.py"),
        "-i", str(clades_dir), "-o", str(mutations_dir), "-r", str(ref_dir),
    ], log_file=log_dir / "step7_mutations.log", description="Step 7: Analyze mutations")


def cleanup_iqtree(trees_dir: Path):
    if not trees_dir.exists():
        return
    removed = 0
    for group_dir in trees_dir.iterdir():
        if not group_dir.is_dir():
            continue
        for f in group_dir.iterdir():
            if f.suffix in (".log", ".bionj", ".mldist", ".iqtree", ".ckp.gz",
                            ".model.gz", ".splits.nex", ".vcf", ".contree"):
                f.unlink(missing_ok=True)
                removed += 1
    if removed:
        print(f"  IQ-TREE: удалено {removed} промежуточных файлов")


def cleanup_mrbayes(mrbayes_dir: Path):
    if not mrbayes_dir.exists():
        return
    removed = 0
    for f in mrbayes_dir.iterdir():
        if not f.is_file():
            continue
        if f.suffix in (".parts", ".mcmc", ".trprobs", ".tstat", ".vstat", ".lstat"):
            f.unlink(missing_ok=True)
            removed += 1
    if removed:
        print(f"  MrBayes: удалено {removed} промежуточных файлов")


def cleanup_grouped(grouped_dir: Path):
    for subdir in ("v", "d", "j"):
        p = grouped_dir / subdir
        if p.exists():
            shutil.rmtree(p)
            print(f"  Удалена папка grouped_by_germlines/{subdir}/ (не используется)")


def main():
    ensure_env()

    parser = argparse.ArgumentParser(description="Run full BCR analysis pipeline.")
    parser.add_argument("-k", "--key", help="Data key (folder under data/)")
    parser.add_argument("-i", "--input", help="Input BCR_data.tsv path (alternative to -k)")
    parser.add_argument("-o", "--output", help="Output report directory (auto-generated if not provided)")
    parser.add_argument("--config", type=Path, default=None,
                        help="JSON config file (see shared/config.py)")
    parser.add_argument("--skip", nargs="*", default=[],
                        choices=ALL_STEPS, help="Steps to skip")
    parser.add_argument("--only", nargs="*", default=[],
                        choices=ALL_STEPS, help="Run only these steps")
    parser.add_argument("--parallel-trees", action="store_true",
                        help="Run IQ-TREE and MrBayes in parallel")
    parser.add_argument("--grouping-strategy", choices=["allele", "gene", "v_only"], default=None,
                        help="VJ grouping: allele (IGHV1-2*01_IGHJ4*01), "
                             "gene (IGHV1-2_IGHJ4, default), v_only (IGHV1-2)")
    parser.add_argument("--min-group-size", type=int, default=None,
                        help="Minimum sequences per group (filter_groups)")
    parser.add_argument("--max-group-size", type=int, default=None,
                        help="Maximum sequences per group (filter_groups)")
    parser.add_argument("--iqtree-model", default=None,
                        help="IQ-TREE substitution model (default: GTR+F+I+G4)")
    args = parser.parse_args()

    # Load config
    if args.config and args.config.exists():
        config = PipelineConfig.load(args.config)
        print(f"Loaded config from {args.config}")
    else:
        config = DEFAULT_CONFIG

    # CLI overrides config
    if args.grouping_strategy is not None:
        config.grouping_strategy = args.grouping_strategy
    if args.min_group_size is not None:
        config.min_group_size = args.min_group_size
    if args.max_group_size is not None:
        config.max_group_size = args.max_group_size
    if args.iqtree_model is not None:
        config.iqtree_model = args.iqtree_model
    if args.parallel_trees:
        config.parallel_trees = True

    # Resolve input
    if args.key:
        input_tsv = DATA_DIR / args.key / "BCR_data.tsv"
        base_name = args.key
    elif args.input:
        input_tsv = Path(args.input)
        base_name = input_tsv.parent.name
    else:
        parser.error("Either -k/--key or -i/--input is required")

    if not input_tsv.exists():
        parser.error(f"Input file not found: {input_tsv}")

    # Resolve output directory
    if args.output:
        report_dir = Path(args.output)
    else:
        timestamp = datetime.now().strftime("%d%m%Y_%H%M%S")
        report_dir = RESULTS_DIR / f"report_{timestamp}"

    report_dir.mkdir(parents=True, exist_ok=True)
    log_dir = report_dir / "logs"
    log_dir.mkdir(exist_ok=True)

    # Save config snapshot
    config.report_dir = str(report_dir)
    config.save(report_dir / "config.json")

    # Pipeline paths
    filtered_tsv = report_dir / "BCR_data_filtered.tsv"
    grouped_dir = report_dir / "grouped_by_germlines"
    grouped_vj_dir = grouped_dir / "vj"
    vj_filtered_dir = grouped_dir / "vj_filtered"
    aligned_dir = report_dir / "aligned_sequences"
    trees_dir = report_dir / "trees"
    mrbayes_dir = report_dir / "mrbayes"
    viz_dir = report_dir / "trees_visualization"
    groups_dir = report_dir / "groups"
    clades_dir = groups_dir / "clades"
    mutations_dir = report_dir / "mutation_tables"

    # Determine steps to run
    if args.only:
        steps_to_run = set(args.only)
    else:
        steps_to_run = set(ALL_STEPS) - set(args.skip)

    # Dependency checks
    if "mutations" in steps_to_run and "clades" not in steps_to_run:
        print("WARNING: mutations requires clades. Adding clades step.")
        steps_to_run.add("clades")
    if "clades" in steps_to_run and "iqtree" not in steps_to_run and "mrbayes" not in steps_to_run:
        print("WARNING: clades requires iqtree and/or mrbayes output.")
    if "viz" in steps_to_run and "iqtree" not in steps_to_run and "mrbayes" not in steps_to_run:
        print("WARNING: viz requires trees from iqtree or mrbayes.")
    if "filter_groups" in steps_to_run and "group" not in steps_to_run:
        print("WARNING: filter_groups requires group. Adding group step.")
        steps_to_run.add("group")
    if "msa" in steps_to_run and "group" not in steps_to_run and "filter_groups" not in steps_to_run:
        print("WARNING: msa requires group or filter_groups output.")

    print(f"\n{'#'*60}")
    print(f"BIOCAD BCR PIPELINE")
    print(f"Input:  {input_tsv}")
    print(f"Output: {report_dir}")
    print(f"Config: grouping={config.grouping_strategy}, "
          f"iqtree_model={config.iqtree_model}, "
          f"min_group={config.min_group_size}")
    print(f"Steps:  {', '.join(s for s in ALL_STEPS if s in steps_to_run)}")
    print(f"{'#'*60}")

    pipeline_start = time.time()
    results = {}

    # Step 1: Filter
    if "filter" in steps_to_run:
        results["filter"] = step1_filter(input_tsv, filtered_tsv, log_dir)
        if not results["filter"]:
            print("Pipeline stopped at step 1")
            return 1

    # Step 2: Group
    if "group" in steps_to_run:
        results["group"] = step2_group(filtered_tsv, grouped_dir, log_dir,
                                       config.grouping_strategy)
        if not results["group"]:
            print("Pipeline stopped at step 2")
            return 1
        cleanup_grouped(grouped_dir)

    # Step 2b: Filter groups by size
    msa_input_dir = grouped_vj_dir
    if "filter_groups" in steps_to_run:
        results["filter_groups"] = step2b_filter_groups(
            grouped_vj_dir, vj_filtered_dir, log_dir,
            config.min_group_size, config.max_group_size)
        if results["filter_groups"]:
            msa_input_dir = vj_filtered_dir
        else:
            print("Pipeline stopped at filter_groups")
            return 1
    else:
        results["filter_groups"] = True

    # Step 3: MSA
    if "msa" in steps_to_run:
        results["msa"] = step3_msa(msa_input_dir, aligned_dir, log_dir)
        if not results["msa"]:
            print("Pipeline stopped at step 3")
            return 1

    # Steps 4a + 4b: Trees
    trees_ok = True
    if "iqtree" in steps_to_run or "mrbayes" in steps_to_run:
        if config.parallel_trees and "iqtree" in steps_to_run and "mrbayes" in steps_to_run:
            import threading
            iqtree_result = [False]
            mrbayes_result = [False]

            def run_iqtree():
                iqtree_result[0] = step4a_iqtree(aligned_dir, trees_dir, log_dir,
                                                 config.iqtree_model)

            def run_mrbayes():
                mrbayes_result[0] = step4b_mrbayes(aligned_dir, mrbayes_dir, log_dir,
                                                   config)

            t1 = threading.Thread(target=run_iqtree)
            t2 = threading.Thread(target=run_mrbayes)
            t1.start()
            t2.start()
            t1.join()
            t2.join()
            results["iqtree"] = iqtree_result[0]
            results["mrbayes"] = mrbayes_result[0]
            trees_ok = iqtree_result[0] or mrbayes_result[0]
        else:
            if "iqtree" in steps_to_run:
                results["iqtree"] = step4a_iqtree(aligned_dir, trees_dir, log_dir,
                                                  config.iqtree_model)
                trees_ok = trees_ok and results["iqtree"]
                cleanup_iqtree(trees_dir)
            if "mrbayes" in steps_to_run:
                results["mrbayes"] = step4b_mrbayes(aligned_dir, mrbayes_dir, log_dir,
                                                    config)
                trees_ok = trees_ok and results["mrbayes"]
                cleanup_mrbayes(mrbayes_dir)

    # Step 5: Viz
    if "viz" in steps_to_run and trees_ok:
        results["viz"] = step5_viz(trees_dir, viz_dir, log_dir)
    elif "viz" in steps_to_run:
        print("Skipping viz: no tree output available")
        results["viz"] = False

    # Step 6: Clades
    if "clades" in steps_to_run:
        has_iqtree = "iqtree" in steps_to_run and results.get("iqtree")
        has_mrbayes = "mrbayes" in steps_to_run and results.get("mrbayes")
        if has_iqtree or has_mrbayes:
            results["clades"] = step6_clades(trees_dir, mrbayes_dir, groups_dir,
                                             aligned_dir, log_dir)
        else:
            print("Skipping clades: no tree output available")
            results["clades"] = False

    # Step 7: Mutations
    if "mutations" in steps_to_run:
        if results.get("clades") and clades_dir.exists() and any(clades_dir.iterdir()):
            results["mutations"] = step7_mutations(clades_dir, mutations_dir,
                                                   DATA_DIR / "references", log_dir)
        else:
            print("Skipping mutations: no clade FASTAs available")
            results["mutations"] = False

    # Summary
    elapsed = time.time() - pipeline_start
    print(f"\n{'#'*60}")
    print(f"PIPELINE COMPLETE in {elapsed:.1f}s")
    print(f"Report: {report_dir}")
    for step, ok in results.items():
        status = "\u2713" if ok else "\u2717"
        print(f"  {status} {step}")
    print(f"{'#'*60}")

    return 0 if all(results.values()) else 1


if __name__ == "__main__":
    sys.exit(main())