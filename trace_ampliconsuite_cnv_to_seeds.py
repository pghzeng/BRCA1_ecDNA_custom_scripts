#!/usr/bin/env python3
"""Trace how AmpliconSuite/AA selects CNV seeds from a CNVkit .cns (variant v3).

Compared to v2, this variant makes the two tracing relaxations configurable:
- Stage1: optionally disable the arm-median-based threshold: Thresh = Med + ccg - 2
- Stage2: optionally disable the BAM coverage-ratio check ("BAM Coverage Check")

Behavior summary:
- Default: preserve the original Stage1/Stage2 checks from ref_code.
- To reproduce v2 behavior: pass both
    --disable_stage1_med_ccg_thresh and --disable_stage2_bam_coverage_check.
- If Stage2 BAM coverage checking remains enabled, --bam is required.

Parameters:
- --cns: input CNVkit .cns file.
- --bam: BAM/CRAM used by the Stage2 BAM coverage-ratio check. Optional only when
    --disable_stage2_bam_coverage_check is set; otherwise required.
- --ref: reference name passed to AA/ref_code, e.g. hg19, GRCh37, GRCh38, GRCh38_viral, mm10, GRCm38.
- --aa_data_repo: AmpliconArchitect AA_DATA_REPO root containing per-reference annotation files.
- --gain: CN gain cutoff used by Stage1/Stage2 filtering. Default: 4.5.
- --cnsize_min: minimum merged seed-cluster size retained at final Stage2 output. Default: 50000.
- --sample / -s: sample name used in output file names and final AA_CNV_SEEDS bed name field.
- --outdir: output directory for final reports and seed bed.
- --tmpdir: optional temp working directory for intermediate files.
- --keep_intermediates: keep the temp working directory under outdir for debugging.
- --disable_stage1_med_ccg_thresh: disable the Stage1 chromosome-arm median threshold term
    and use only ccg-based thresholding.
- --disable_stage2_bam_coverage_check: disable the Stage2 BAM coverage-ratio filter;
    when enabled, Stage2 does not require --bam.

Typical usage:
- Original behavior (keep both checks):
    trace_ampliconsuite_cnv_to_seeds_v3.py --cns sample.cns --bam sample.bam --ref GRCh38 \
            --aa_data_repo /path/to/AA_DATA_REPO -s sample --outdir out
- v2-like behavior (disable both checks, no bam required):
    trace_ampliconsuite_cnv_to_seeds_v3.py --cns sample.cns --ref GRCh38 \
            --aa_data_repo /path/to/AA_DATA_REPO -s sample --outdir out \
            --disable_stage1_med_ccg_thresh --disable_stage2_bam_coverage_check
- Mixed behavior (keep Stage1 original threshold, disable only BAM coverage check):
    trace_ampliconsuite_cnv_to_seeds_v3.py --cns sample.cns --ref GRCh38 \
            --aa_data_repo /path/to/AA_DATA_REPO -s sample --outdir out \
            --disable_stage2_bam_coverage_check

Outputs (in --outdir):
- <sample>.trace_SEED_stage1.tsv
- <sample>.trace_SEED_stage2.tsv
- <sample>.AA_CNV_SEEDS.bed

Implementation note:
This script executes the workspace-copied, instrumented AA/paalib code under ref_code/.
Reasons are emitted by those instrumented scripts into AA_TRACE_LOG_FILE and then summarized here.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from typing import Dict, List, Tuple


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_CODE_DIR = os.path.join(BASE_DIR, "ref_code")
sys.path.append(REF_CODE_DIR)

import cnv_prefilter  # noqa: E402


def cns_to_cnv_calls_bed(cns_path: str, bed_path: str) -> None:
    """Match the pipeline's awk conversion: NR>1 {print chr start end CNVkit 2*(2^log2)}"""
    with open(cns_path) as f, open(bed_path, "w") as out:
        _ = next(f)  # header
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            log2 = float(parts[4])
            cn = 2 * (2 ** log2)
            out.write(f"{chrom}\t{start}\t{end}\tCNVkit\t{cn}\n")


def read_tsv_events(trace_path: str) -> List[Dict[str, str]]:
    with open(trace_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def index_events_by_coords(events: List[Dict[str, str]]) -> Dict[Tuple[str, int, int], List[Dict[str, str]]]:
    d: Dict[Tuple[str, int, int], List[Dict[str, str]]] = {}
    for e in events:
        key = (e["Chrom"], int(e["Start"]), int(e["End"]))
        d.setdefault(key, []).append(e)
    return d


def load_unfiltered_gains_by_chrom(unfiltered_bed: str) -> Dict[str, List[Tuple[int, int, float]]]:
    by_chrom: Dict[str, List[Tuple[int, int, float]]] = defaultdict(list)
    with open(unfiltered_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(float(parts[1]))
            end = int(float(parts[2]))
            cn = float(parts[3])
            by_chrom[chrom].append((start, end, cn))
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: (x[0], x[1]))
    return by_chrom


def overlaps_bed(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def annotate_penalty_mult(details: str, seg_len: int) -> str:
    """If FailDetails contains PenaltyMult=xx, append a parenthesized reason.

    PenaltyMult is emitted by ref_code/cnv_prefilter.py to inflate ccg (and thus the Stage1 threshold)
    for long segments and/or segments overlapping long continuous high-CN runs.
    """
    if not details or details == "-":
        return details
    if "PenaltyMult=" not in details:
        return details
    if re.search(r"PenaltyMult=[0-9.]+\(", details):
        return details

    mult_m = re.search(r"PenaltyMult=([0-9.]+)", details)
    if not mult_m:
        return details

    try:
        mult = float(mult_m.group(1))
    except ValueError:
        return details

    max_run_len = None
    mrl_m = re.search(r"MaxHighCNRunLen=([0-9]+)", details)
    if mrl_m:
        try:
            max_run_len = int(mrl_m.group(1))
        except ValueError:
            max_run_len = None

    reasons: List[str] = []
    if seg_len > 20_000_000 and mult >= 2.0 - 1e-9:
        reasons.append(f"Len>20Mb (Len={seg_len})")
    if (max_run_len is not None) and (max_run_len > 10_000_000 or abs(max_run_len - 10_000_000) <= 1):
        if mult >= 1.5 - 1e-9:
            reasons.append(f"HighCNRun>=10Mb (MaxHighCNRunLen={max_run_len})")

    if not reasons:
        reasons.append("Stage1 penalty applied")

    reason_str = "; ".join(reasons)
    return re.sub(r"PenaltyMult=([0-9.]+)", lambda m: f"PenaltyMult={m.group(1)}({reason_str})", details, count=1)


def chrom_sort_key(chrom: str) -> Tuple[int, int, str]:
    c = chrom
    if c.startswith("chr"):
        c = c[3:]
    if c.isdigit():
        return (0, int(c), "")
    upper = c.upper()
    if upper == "X":
        return (1, 23, "")
    if upper == "Y":
        return (1, 24, "")
    if upper in ("M", "MT"):
        return (1, 25, "")
    return (2, 10**9, chrom)


def sort_tsv_by_bed_coords(rows: List[List[str]], chrom_i: int, start_i: int, end_i: int) -> List[List[str]]:
    def key(r: List[str]):
        return (chrom_sort_key(r[chrom_i]), int(float(r[start_i])), int(float(r[end_i])))

    return sorted(rows, key=key)


def rewrite_seeds_bed(in_bed: str, out_bed: str, sample: str) -> None:
    """Rewrite ref_code amplified_intervals output to: chrom start end sample cngain.

    ref_code/amplified_intervals.py writes: chrom start end <cn_float> <source_path>
    We keep coordinates, move CN float to last column, and set name=sample.
    """
    with open(in_bed) as src, open(out_bed, "w") as dst:
        for line in src:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            cn = parts[3]
            dst.write(f"{chrom}\t{start}\t{end}\t{sample}\t{cn}\n")


def get_ref_centromeres(ref_name: str, aa_data_repo: str) -> Dict[str, Tuple[str, str]]:
    centromere_dict: Dict[str, Tuple[str, str]] = {}
    fnameD = {
        "GRCh38": "GRCh38_centromere.bed",
        "GRCh37": "human_g1k_v37_centromere.bed",
        "hg19": "hg19_centromere.bed",
        "mm10": "mm10_centromere.bed",
        "GRCm38": "GRCm38_centromere.bed",
        "GRCh38_viral": "GRCh38_centromere.bed",
    }

    path = os.path.join(aa_data_repo, ref_name, fnameD.get(ref_name, "GRCh38_centromere.bed"))
    if not os.path.exists(path):
        return {}

    with open(path) as infile:
        for line in infile:
            fields = line.strip().split()
            if fields:
                centromere_dict[fields[0]] = (fields[1], fields[2])
    return centromere_dict


def get_ref_sizes(ref_name: str, aa_data_repo: str) -> Dict[str, int]:
    file_list_path = os.path.join(aa_data_repo, ref_name, "file_list.txt")
    chr_len_filename = ""

    if os.path.exists(file_list_path):
        with open(file_list_path) as f:
            for line in f:
                parts = line.strip().split()
                if parts and parts[0] == "chrLen_file":
                    chr_len_filename = parts[1]
                    break

    if not chr_len_filename:
        chr_len_filename = f"{ref_name}full.fa.fai"

    path = os.path.join(aa_data_repo, ref_name, chr_len_filename)
    sizes: Dict[str, int] = {}
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    sizes[parts[0]] = int(parts[1])
    return sizes


def configure_trace_mode(args: argparse.Namespace) -> Dict[str, str]:
    os.environ["AA_DISABLE_STAGE1_MED_CCGB_THRESH"] = "1" if args.disable_stage1_med_ccg_thresh else "0"
    os.environ["AA_DISABLE_STAGE2_BAM_COVERAGE_CHECK"] = "1" if args.disable_stage2_bam_coverage_check else "0"
    env = os.environ.copy()
    return env


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cns", required=True)
    p.add_argument(
        "--bam",
        default="",
        help="Optional BAM/CRAM. Required when Stage2 BAM coverage-ratio checking is enabled.",
    )
    p.add_argument("--ref", required=True)
    p.add_argument("--aa_data_repo", required=True)
    p.add_argument("--gain", type=float, default=4.5)
    p.add_argument("--cnsize_min", type=int, default=50000)
    p.add_argument("-s", "--sample", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--tmpdir", default="")
    p.add_argument("--keep_intermediates", action="store_true")
    p.add_argument(
        "--disable_stage1_med_ccg_thresh",
        action="store_true",
        help="Disable the Stage1 arm-median-based threshold (Thresh = Med + ccg - 2).",
    )
    p.add_argument(
        "--disable_stage2_bam_coverage_check",
        action="store_true",
        help="Disable the Stage2 BAM coverage-ratio check. When omitted, --bam is required.",
    )
    args = p.parse_args()

    if not args.disable_stage2_bam_coverage_check and not args.bam:
        p.error("--bam is required unless --disable_stage2_bam_coverage_check is set")

    return args


def main() -> None:
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    tmp_parent = args.tmpdir if args.tmpdir else None
    tmp_cm = tempfile.TemporaryDirectory(dir=tmp_parent)
    workdir = tmp_cm.name

    base = os.path.splitext(os.path.basename(args.cns))[0]
    cnv_calls_bed = os.path.join(workdir, f"{base}_CNV_CALLS.bed")
    trace_path = os.path.join(workdir, "trace_report_v3.tsv")

    os.environ["AA_DATA_REPO"] = args.aa_data_repo
    os.environ["AA_TRACE_LOG_FILE"] = trace_path
    env = configure_trace_mode(args)

    # cnv_prefilter caches TRACE_LOG_FILE at import time; update explicitly.
    cnv_prefilter.TRACE_LOG_FILE = trace_path

    if os.path.exists(trace_path):
        os.remove(trace_path)
    with open(trace_path, "w") as f:
        f.write("Chrom\tStart\tEnd\tGene\tLog2\tCN\tStatus\tStage\tReason\tDetails\n")

    # 1) Convert CNS -> CNV_CALLS.bed
    cns_to_cnv_calls_bed(args.cns, cnv_calls_bed)

    # 2) Stage1: cnv_prefilter -> unfiltered_gains
    centromere_dict = get_ref_centromeres(args.ref, args.aa_data_repo)
    chr_sizes = get_ref_sizes(args.ref, args.aa_data_repo)
    unfiltered_gains_bed = cnv_prefilter.prefilter_bed(
        cnv_calls_bed,
        args.ref,
        centromere_dict,
        chr_sizes,
        args.gain,
        workdir,
    )

    # 3) Stage2: amplified_intervals -> seeds
    seeds_prefix = os.path.join(workdir, f"{base}_AA_CNV_SEEDS")
    cmd = [
        sys.executable,
        os.path.join(REF_CODE_DIR, "amplified_intervals.py"),
        "--bed",
        unfiltered_gains_bed,
        "--ref",
        args.ref,
        "--gain",
        str(args.gain),
        "--cnsize_min",
        str(args.cnsize_min),
        "--out",
        seeds_prefix,
    ]
    if args.bam:
        cmd.extend(["--bam", args.bam])
    subprocess.check_call(cmd, env=env, cwd=workdir)

    seeds_bed = seeds_prefix + ".bed"
    out_seeds = os.path.join(args.outdir, f"{args.sample}_AA_CNV_SEEDS.bed")
    if os.path.exists(seeds_bed):
        rewrite_seeds_bed(seeds_bed, out_seeds, args.sample)

    # 4) Summaries from trace log
    events = read_tsv_events(trace_path)
    events_by_coord = index_events_by_coords(events)
    unfiltered_by_chrom = load_unfiltered_gains_by_chrom(unfiltered_gains_bed)

    # Stage2 report: one line per unfiltered_gains interval.
    stage2_rows: List[List[str]] = []
    with open(unfiltered_gains_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(float(parts[1]))
            end = int(float(parts[2]))
            cn = parts[3]
            key = (chrom, start, end)
            evs = events_by_coord.get(key, [])

            seed_evs = [e for e in evs if e.get("Status") == "SEED"]
            if seed_evs:
                stage2_rows.append([chrom, str(start), str(end + 1), cn, "IN_SEED_CLUSTER", "Final", "-", seed_evs[0].get("Details", "-")])
                continue

            cand = [e for e in evs if e.get("Status") == "FILTERED" and e.get("Stage") in ("AmplifiedInit", "Amplified", "Final")]
            if cand:
                e0 = cand[0]
                stage2_rows.append([
                    chrom,
                    str(start),
                    str(end + 1),
                    cn,
                    "FILTERED",
                    e0.get("Stage", "-"),
                    e0.get("Reason", "-"),
                    e0.get("Details", "-"),
                ])
            else:
                stage2_rows.append([chrom, str(start), str(end + 1), cn, "NO_STAGE2_DECISION", "-", "-", "-"])

    stage2_rows = sort_tsv_by_bed_coords(stage2_rows, 0, 1, 2)
    out_stage2 = os.path.join(args.outdir, f"{args.sample}.trace_SEED_stage2.tsv")
    with open(out_stage2, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["Chrom", "Start", "End", "CN", "Stage2Status", "FailStage", "FailReason", "FailDetails"])
        w.writerows(stage2_rows)

    # Stage1 report: one line per .cns record.
    stage1_rows: List[List[str]] = []
    with open(args.cns) as f:
        _ = next(f)
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            log2 = float(parts[4])
            cn = 2 * (2 ** log2)

            key = (chrom, start, end)
            evs = events_by_coord.get(key, [])
            pre = [e for e in evs if e.get("Stage") == "Prefilter"]

            stage1_status = "NO_PREFILTER_DECISION"
            kept_subtype = "-"
            fail_reason = "-"
            fail_details = "-"
            match_count = 0
            match_str = "-"

            filtered = [e for e in pre if e.get("Status") == "FILTERED"]
            if filtered:
                stage1_status = "FILTERED"
                fail_reason = filtered[0].get("Reason", "-")
                fail_details = annotate_penalty_mult(filtered[0].get("Details", "-"), end - start)
            else:
                kept = [e for e in pre if e.get("Status") in ("KEPT_S1", "KEPT_S1_MOD")]
                if kept:
                    last_status = kept[-1].get("Status", "KEPT_S1")
                    kept_subtype = "INTACT" if last_status == "KEPT_S1" else "SPLIT_TRIM"

                    hits = []
                    for us, ue, ucn in unfiltered_by_chrom.get(chrom, []):
                        if overlaps_bed(start, end, us, ue):
                            hits.append((us, ue, ucn))
                    match_count = len(hits)
                    if hits:
                        stage1_status = "IN_UNFILTERED_GAINS"
                        show = hits[:5]
                        parts2 = [f"{chrom}:{s}-{e}({cnv:.3f})" for s, e, cnv in show]
                        suffix = "" if len(hits) <= 5 else f"; +{len(hits)-5} more"
                        match_str = ";".join(parts2) + suffix
                    else:
                        stage1_status = "FILTERED"
                        if cn <= args.gain:
                            fail_reason = "CN <= Gain (final Stage1 cutoff)"
                            fail_details = f"CN={cn:.3f} <= gain={args.gain}"
                        else:
                            fail_reason = "Not present in unfiltered_gains"
                            fail_details = "Passed intermediate Prefilter but did not overlap any final unfiltered_gains interval"

            stage1_rows.append([
                chrom,
                str(start),
                str(end),
                f"{cn}",
                stage1_status,
                kept_subtype,
                fail_reason,
                fail_details,
                str(match_count),
                match_str,
            ])

    stage1_rows = sort_tsv_by_bed_coords(stage1_rows, 0, 1, 2)
    out_stage1 = os.path.join(args.outdir, f"{args.sample}.trace_SEED_stage1.tsv")
    with open(out_stage1, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow([
            "Chrom",
            "Start",
            "End",
            "CN",
            "Stage1Status",
            "Stage1KeptSubtype",
            "FailReason",
            "FailDetails",
            "UnfilteredGainsMatchCount",
            "UnfilteredGainsMatches",
        ])
        w.writerows(stage1_rows)

    if args.keep_intermediates:
        keep_path = os.path.join(args.outdir, f"{args.sample}.trace_intermediates")
        if not os.path.exists(keep_path):
            os.rename(workdir, keep_path)
            tmp_cm.cleanup = lambda: None  # type: ignore

    tmp_cm.cleanup()


if __name__ == "__main__":
    main()