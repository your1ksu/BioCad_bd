#!/usr/bin/env python3
"""
Master pipeline runner for BCR analysis.

Usage:
    python scripts/run_pipeline.py -k batch1
    python scripts/run_pipeline.py -i data/batch1/BCR_data.tsv -o results/report_17072026_143000
"""

import argparse
import os
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


def run_cmd(cmd, cwd=None, log_file=None, description=None):
    """Run command and log output."""
    if description:
        print(f"\n{'='*60}")
        print(f"STEP: {description}")
        print(f"CMD:  {' '.join(cmd)}")
        print(f"{'='*60}")

    start = time.time()
    elapsed = 0.0
    try:
        if log_file:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            with open(log_file, "w") as lf:
                lf.write(f"CMD: {' '.join(cmd)}\n")
                lf.write(f"CWD: {cwd or REPO_ROOT}\n")
                lf.write(f"START: {datetime.now()}\n\n")
                result = subprocess.run(
                    cmd, cwd=cwd or REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT, text=True
                )
        else:
            result = subprocess.run(
                cmd, cwd=cwd or REPO_ROOT, capture_output=True, text=True
            )
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(result.stderr, file=sys.stderr)

        elapsed = time.time() - start
        if result.returncode == 0:
            print(f"✓ DONE in {elapsed:.1f}s")
            return True
        else:
            print(f"✗ FAILED (exit {result.returncode}) in {elapsed:.1f}s")
            if log_file:
                print(f"  See log: {log_file}")
            return False
    except subprocess.TimeoutExpired:
        print(f"✗ TIMEOUT after {elapsed:.1f}s")
        return False
    except Exception as e:
        print(f"✗ ERROR: {e}")
        return False


def step1_filter(input_tsv, output_tsv, log_dir):
    """Step 1: Filter BCR data."""
    return run_cmd([
        sys.executable, str(SCRIPTS_DIR / "01_filter_sequences" / "filter_sequences.py"),
        "-i", str(input_tsv),
        "-o", str(output_tsv),
        "-r", str(REFERENCES_DIR),
    ], log_file=log_dir / "step1_filter.log", description="Step 1: Filter BCR data")


def step2_group(filtered_tsv, output_dir, log_dir):
    """Step 2: Group by germlines."""
    return run_cmd([
        sys.executable, str(SCRIPTS_DIR / "02_group_by_germlines" / "group_by_germlines.py"),
        "-i", str(filtered_tsv),
        "-o", str(output_dir),
        "-r", str(REFERENCES_DIR),
    ], log_file=log_dir / "step2_group.log", description="Step 2: Group by germlines")


def step3_msa(grouped_vj_dir, output_dir, log_dir):
    """Step 3: Multiple sequence alignment."""
    return run_cmd([
        sys.executable, str(SCRIPTS_DIR / "03_multiple_alignment" / "multiple_alignment.py"),
        "-i", str(grouped_vj_dir),
        "-o", str(output_dir),
    ], log_file=log_dir / "step3_msa.log", description="Step 3: MSA (MAFFT)")


def step4a_iqtree(aligned_dir, trees_dir, log_dir):
    """Step 4a: IQ-TREE (ML trees)."""
    return run_cmd([
        "bash", str(SCRIPTS_DIR / "04a_build_trees_iqtree" / "build_trees_iqtree.sh"),
        str(aligned_dir), str(trees_dir),
    ], log_file=log_dir / "step4a_iqtree.log", description="Step 4a: IQ-TREE (ML)")


def step4b_mrbayes(aligned_dir, mrbayes_dir, log_dir):
    """Step 4b: MrBayes (Bayesian trees)."""
    return run_cmd([
        "conda", "run", "-n", "biocad_pipeline", "python",
        str(SCRIPTS_DIR / "04b_build_trees_mrbayes" / "build_trees_mrbayes.py"),
        str(aligned_dir),
        "--out", str(mrbayes_dir),
    ], log_file=log_dir / "step4b_mrbayes.log", description="Step 4b: MrBayes (Bayesian)")


def step5_viz(trees_dir, viz_dir, log_dir):
    """Step 5: Visualize trees."""
    return run_cmd([
        "conda", "run", "-n", "biocad_pipeline", "bash",
        str(SCRIPTS_DIR / "visualize_trees" / "visualize_trees.sh"),
        str(trees_dir), str(viz_dir),
    ], log_file=log_dir / "step5_viz.log", description="Step 5: Tree visualization")


def step6_clades(trees_dir, mrbayes_dir, groups_dir, aligned_dir, log_dir):
    """Step 6: Confident clades report + extract clade FASTAs."""
    clades_dir = groups_dir / "clades"
    return run_cmd([
        "conda", "run", "-n", "biocad_pipeline", "python",
        str(SCRIPTS_DIR / "05_clade_search" / "clade_search.py"),
        "--iqtree-dir", str(trees_dir),
        "--mrbayes-dir", str(mrbayes_dir),
        "--out", str(groups_dir / "report.json"),
        "--clades-fasta-dir", str(clades_dir),
        "--aligned-dir", str(aligned_dir),
    ], log_file=log_dir / "step6_clades.log", description="Step 6: Confident clades report")


def step7_mutations(clades_dir, mutations_dir, ref_dir, log_dir):
    """Step 7: Analyze mutations."""
    return run_cmd([
        "bash", str(SCRIPTS_DIR / "06_analyze_mutations" / "analyze_mutations.sh"),
        str(clades_dir), str(mutations_dir), str(ref_dir),
    ], log_file=log_dir / "step7_mutations.log", description="Step 7: Analyze mutations")


def main():
    parser = argparse.ArgumentParser(description="Run full BCR analysis pipeline.")
    parser.add_argument("-k", "--key", help="Data key (folder under data/)")
    parser.add_argument("-i", "--input", help="Input BCR_data.tsv path (alternative to -k)")
    parser.add_argument("-o", "--output", help="Output report directory (auto-generated if not provided)")
    parser.add_argument("--skip", nargs="*", default=[],
                        choices=["filter", "group", "msa", "iqtree", "mrbayes", "viz", "clades", "mutations"],
                        help="Steps to skip")
    parser.add_argument("--only", nargs="*", default=[],
                        choices=["filter", "group", "msa", "iqtree", "mrbayes", "viz", "clades", "mutations"],
                        help="Run only these steps")
    parser.add_argument("--parallel-trees", action="store_true",
                        help="Run IQ-TREE and MrBayes in parallel (experimental)")
    args = parser.parse_args()

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

    # Pipeline paths
    filtered_tsv = report_dir / "BCR_data_filtered.tsv"
    grouped_dir = report_dir / "grouped_by_germlines"
    grouped_vj_dir = grouped_dir / "vj"
    aligned_dir = report_dir / "aligned_sequences"
    trees_dir = report_dir / "trees"
    mrbayes_dir = report_dir / "mrbayes"
    viz_dir = report_dir / "trees_visualization"
    groups_dir = report_dir / "groups"
    clades_dir = groups_dir / "clades"
    mutations_dir = report_dir / "mutation_tables"

    # Determine steps to run
    all_steps = ["filter", "group", "msa", "iqtree", "mrbayes", "viz", "clades", "mutations"]
    if args.only:
        steps_to_run = set(args.only)
    else:
        steps_to_run = set(all_steps) - set(args.skip)

    # Dependency check
    if "mutations" in steps_to_run and "clades" not in steps_to_run:
        print("WARNING: mutations step requires clades output. Adding clades step.")
        steps_to_run.add("clades")
    if "clades" in steps_to_run and "iqtree" not in steps_to_run and "mrbayes" not in steps_to_run:
        print("WARNING: clades step requires iqtree and/or mrbayes output.")
    if "viz" in steps_to_run and "iqtree" not in steps_to_run and "mrbayes" not in steps_to_run:
        print("WARNING: viz step requires trees from iqtree or mrbayes.")

    print(f"\n{'#'*60}")
    print(f"BIOCAD BCR PIPELINE")
    print(f"Input:  {input_tsv}")
    print(f"Output: {report_dir}")
    print(f"Steps:  {', '.join(s for s in all_steps if s in steps_to_run)}")
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
        results["group"] = step2_group(filtered_tsv, grouped_dir, log_dir)
        if not results["group"]:
            print("Pipeline stopped at step 2")
            return 1

    # Step 3: MSA
    if "msa" in steps_to_run:
        results["msa"] = step3_msa(grouped_vj_dir, aligned_dir, log_dir)
        if not results["msa"]:
            print("Pipeline stopped at step 3")
            return 1

    # Steps 4a + 4b: Trees (can run in parallel)
    trees_ok = True
    if "iqtree" in steps_to_run or "mrbayes" in steps_to_run:
        if args.parallel_trees and "iqtree" in steps_to_run and "mrbayes" in steps_to_run:
            # Parallel execution
            import threading
            iqtree_result = [False]
            mrbayes_result = [False]

            def run_iqtree():
                iqtree_result[0] = step4a_iqtree(aligned_dir, trees_dir, log_dir)

            def run_mrbayes():
                mrbayes_result[0] = step4b_mrbayes(aligned_dir, mrbayes_dir, log_dir)

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
            # Sequential
            if "iqtree" in steps_to_run:
                results["iqtree"] = step4a_iqtree(aligned_dir, trees_dir, log_dir)
                trees_ok = trees_ok and results["iqtree"]
            if "mrbayes" in steps_to_run:
                results["mrbayes"] = step4b_mrbayes(aligned_dir, mrbayes_dir, log_dir)
                trees_ok = trees_ok and results["mrbayes"]

    # Step 5: Viz (needs trees)
    if "viz" in steps_to_run and trees_ok:
        results["viz"] = step5_viz(trees_dir, viz_dir, log_dir)
    elif "viz" in steps_to_run:
        print("Skipping viz: no tree output available")
        results["viz"] = False

    # Step 6: Clades (needs trees + mrbayes)
    if "clades" in steps_to_run:
        has_iqtree = "iqtree" in steps_to_run and results.get("iqtree")
        has_mrbayes = "mrbayes" in steps_to_run and results.get("mrbayes")
        if has_iqtree or has_mrbayes:
            results["clades"] = step6_clades(trees_dir, mrbayes_dir, groups_dir, aligned_dir, log_dir)
        else:
            print("Skipping clades: no tree output available")
            results["clades"] = False

    # Step 7: Mutations (needs clades)
    if "mutations" in steps_to_run:
        if results.get("clades") and clades_dir.exists() and any(clades_dir.iterdir()):
            results["mutations"] = step7_mutations(clades_dir, mutations_dir, DATA_DIR / "references", log_dir)
        else:
            print("Skipping mutations: no clade FASTAs available")
            results["mutations"] = False

    # Summary
    elapsed = time.time() - pipeline_start
    print(f"\n{'#'*60}")
    print(f"PIPELINE COMPLETE in {elapsed:.1f}s")
    print(f"Report: {report_dir}")
    for step, ok in results.items():
        status = "✓" if ok else "✗"
        print(f"  {status} {step}")
    print(f"{'#'*60}")

    return 0 if all(results.values()) else 1


if __name__ == "__main__":
    sys.exit(main())