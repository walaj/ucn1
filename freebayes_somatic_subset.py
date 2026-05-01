#!/usr/bin/env python3
"""
freebayes_somatic_subset.py
===========================

Filter FreeBayes tumor-vs-normal VCFs down to a high-confidence somatic-only
subset using a uniform set of cutoffs. Every cutoff is exposed as a CLI flag so
you can sweep up and down a ROC curve by re-running with different values.

Inputs
------
One or more FreeBayes VCFs (raw .vcf or bgzipped .vcf.gz). The script does not
require the index. Tumor and normal samples are detected by regex match against
the sample names in the #CHROM line (defaults: "tumor" and "normal", case
insensitive). Override with --tumor-regex / --normal-regex if your samples are
named otherwise.

Outputs
-------
For each input file, a filtered VCF is written to --output-dir (or alongside
the input, if --in-place is given) with suffix ``.somatic_subset.vcf`` (or
``.vcf.gz`` if --gzip-output). A summary TSV (``somatic_subset_summary.tsv``)
is written to --output-dir capturing per-file counts and the cutoffs used. The
exact cutoffs are also recorded as ``##sarek_subset_*`` header lines in each
output VCF so the result is self-documenting.

Filtering logic (per record, per ALT allele)
--------------------------------------------
A record is KEPT only when ALL of the following are true:

  1. FILTER is "PASS" or "." (or any allowed value in --pass-values).
  2. QUAL >= --qual.
  3. (default) The record is biallelic SNV/indel, i.e. exactly one ALT.
     Disable with --keep-multiallelic.
  4. (optional) REF and ALT are both length 1 if --snvs-only.
  5. Tumor genotype is non-reference and not missing.
     Normal genotype is reference (0/0) or missing-but-low-VAF.
  6. Tumor depth (DP):    --min-tumor-dp <= DP <= --max-tumor-dp
     Normal depth (DP):   --min-normal-dp <= DP <= --max-normal-dp
  7. Tumor ALT depth >= --min-tumor-ad
     Normal ALT depth <= --max-normal-ad
  8. Tumor VAF >= --min-tumor-vaf
     Normal VAF <= --max-normal-vaf

Optional add-ons (off by default; enable with corresponding flag):
  9. --min-mqm   : INFO/MQM   (mean MQ of ALT reads)
 10. --min-mqmr  : INFO/MQMR  (mean MQ of REF reads)
 11. --max-sap   : INFO/SAP   (Phred strand-bias on ALT)
 12. --max-srp   : INFO/SRP   (Phred strand-bias on REF)
 13. --max-epp   : INFO/EPP   (Phred read-end placement bias on ALT)
 14. --min-paired: INFO/PAIRED (fraction of ALT reads in proper pairs)

Re-running for ROC sweeps
-------------------------
Loosen / tighten any of these knobs and re-run. Example bash sweep::

  for vaf in 0.02 0.05 0.10 0.15 0.20; do
    for q in 1 10 30 100; do
      ./freebayes_somatic_subset.py *.freebayes.vcf \\
          -o subsets/vaf${vaf}_q${q} \\
          --min-tumor-vaf ${vaf} --qual ${q}
    done
  done

The summary TSV in each output dir lets you plot kept-count vs. cutoff to
trace the curve.

Author: scripted for Jeremiah Wala. No external dependencies; pure stdlib.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import glob
import gzip
import io
import os
import re
import shlex
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def smart_open_read(path: str) -> io.TextIOBase:
    """Open .vcf or .vcf.gz transparently for streaming text reads."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def smart_open_write(path: str, gzip_output: bool) -> io.TextIOBase:
    """Open output file as plain text or gzip depending on flag/suffix."""
    if gzip_output or path.endswith(".gz"):
        if not path.endswith(".gz"):
            path = path + ".gz"
        return gzip.open(path, "wt")
    return open(path, "wt")


def out_path_for(input_path: str, out_dir: Optional[str], in_place: bool,
                 gzip_output: bool) -> str:
    base = os.path.basename(input_path)
    base = re.sub(r"\.vcf(\.gz)?$", "", base, flags=re.IGNORECASE)
    suffix = ".somatic_subset.vcf"
    if gzip_output:
        suffix += ".gz"
    if in_place:
        return os.path.join(os.path.dirname(input_path) or ".", base + suffix)
    assert out_dir is not None
    os.makedirs(out_dir, exist_ok=True)
    return os.path.join(out_dir, base + suffix)


# ---------------------------------------------------------------------------
# VCF parsing helpers
# ---------------------------------------------------------------------------

INFO_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_.]*)=([^;]*)")


def parse_info(info_field: str) -> Dict[str, str]:
    if info_field == "." or not info_field:
        return {}
    out: Dict[str, str] = {}
    for kv in info_field.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = "1"  # flag-only INFO key
    return out


def parse_sample(fmt_keys: List[str], sample_field: str) -> Dict[str, str]:
    vals = sample_field.split(":")
    return dict(zip(fmt_keys, vals))


def to_int(s: Optional[str], default: Optional[int] = None) -> Optional[int]:
    if s is None or s == "." or s == "":
        return default
    try:
        return int(s)
    except ValueError:
        return default


def to_float(s: Optional[str], default: Optional[float] = None) -> Optional[float]:
    if s is None or s == "." or s == "":
        return default
    try:
        return float(s)
    except ValueError:
        return default


def first_int(s: Optional[str], default: Optional[int] = None) -> Optional[int]:
    """For FORMAT fields with comma-separated values, return the first int."""
    if s is None or s == "." or s == "":
        return default
    try:
        return int(s.split(",")[0])
    except ValueError:
        return default


def gt_is_ref(gt: Optional[str]) -> bool:
    """0/0 or 0|0; treats './.' as NOT ref (missing)."""
    if not gt:
        return False
    parts = re.split(r"[/|]", gt)
    return all(p == "0" for p in parts)


def gt_is_missing(gt: Optional[str]) -> bool:
    if not gt:
        return True
    parts = re.split(r"[/|]", gt)
    return all(p == "." for p in parts)


def gt_is_alt(gt: Optional[str]) -> bool:
    """At least one non-zero, non-missing allele."""
    if not gt:
        return False
    parts = re.split(r"[/|]", gt)
    if any(p == "." for p in parts):
        return False
    return any(p != "0" and p.isdigit() and int(p) > 0 for p in parts)


def detect_sample_indices(header_line: str,
                          tumor_re: re.Pattern,
                          normal_re: re.Pattern,
                          tumor_hint: Optional[str],
                          normal_hint: Optional[str],
                          ) -> Tuple[int, int, str, str]:
    cols = header_line.rstrip("\n").split("\t")
    if len(cols) < 11:
        raise ValueError(f"Header has only {len(cols)} columns; need >= 11 for T/N.")
    samples = cols[9:]

    def pick(rx: re.Pattern, hint: Optional[str], label: str) -> int:
        if hint is not None:
            if hint not in samples:
                raise ValueError(f"--{label}-sample {hint!r} not in #CHROM samples {samples}")
            return samples.index(hint)
        hits = [i for i, s in enumerate(samples) if rx.search(s)]
        if len(hits) == 1:
            return hits[0]
        if len(hits) == 0:
            raise ValueError(
                f"No sample in {samples} matches --{label}-regex {rx.pattern!r}. "
                f"Use --{label}-sample to set explicitly."
            )
        raise ValueError(
            f"Multiple samples in {samples} match --{label}-regex {rx.pattern!r}: "
            f"{[samples[i] for i in hits]}. Use --{label}-sample to disambiguate."
        )

    t = pick(tumor_re, tumor_hint, "tumor")
    n = pick(normal_re, normal_hint, "normal")
    if t == n:
        raise ValueError(f"Tumor and normal regex resolved to the same sample: {samples[t]!r}")
    return 9 + t, 9 + n, samples[t], samples[n]


# ---------------------------------------------------------------------------
# Per-record evaluation
# ---------------------------------------------------------------------------

REJECT_REASONS = [
    "filter_not_pass",
    "low_qual",
    "multiallelic",
    "not_snv",
    "tumor_gt_not_alt",
    "normal_gt_not_ref",
    "low_tumor_dp",
    "high_tumor_dp",
    "low_normal_dp",
    "high_normal_dp",
    "low_tumor_ad",
    "high_normal_ad",
    "low_tumor_vaf",
    "high_normal_vaf",
    "low_mqm",
    "low_mqmr",
    "high_sap",
    "high_srp",
    "high_epp",
    "low_paired",
    "parse_error",
]


def evaluate_record(fields: List[str], t_idx: int, n_idx: int,
                    args: argparse.Namespace,
                    counters: Dict[str, int]) -> bool:
    """Return True if the record passes all filters."""
    try:
        chrom, pos, _id, ref, alt, qual, flt, info_str, fmt = fields[:9]
    except ValueError:
        counters["parse_error"] += 1
        return False

    # 1. FILTER
    if flt not in args.pass_values_set:
        counters["filter_not_pass"] += 1
        return False

    # 2. QUAL
    qv = to_float(qual, default=None)
    if qv is None or qv < args.qual:
        counters["low_qual"] += 1
        return False

    # 3. multiallelic
    alts = alt.split(",")
    if len(alts) > 1 and not args.keep_multiallelic:
        counters["multiallelic"] += 1
        return False

    # 4. SNV-only
    if args.snvs_only:
        if len(ref) != 1 or any(len(a) != 1 for a in alts):
            counters["not_snv"] += 1
            return False

    # FORMAT/sample parsing
    fmt_keys = fmt.split(":")
    try:
        t_field = fields[t_idx]
        n_field = fields[n_idx]
    except IndexError:
        counters["parse_error"] += 1
        return False
    t_s = parse_sample(fmt_keys, t_field)
    n_s = parse_sample(fmt_keys, n_field)

    # 5. Genotype constraints (single-ALT semantics; for multi-ALT we relax)
    t_gt = t_s.get("GT")
    n_gt = n_s.get("GT")
    if not gt_is_alt(t_gt):
        counters["tumor_gt_not_alt"] += 1
        return False
    if not (gt_is_ref(n_gt) or gt_is_missing(n_gt)):
        counters["normal_gt_not_ref"] += 1
        return False

    # Pull AD / AO / RO / DP from FORMAT.
    def alt_count(s: Dict[str, str]) -> Optional[int]:
        # Prefer AD (Number=R, "ref,alt1,alt2,...") then AO (alt counts).
        ad = s.get("AD")
        if ad and ad != ".":
            parts = ad.split(",")
            if len(parts) >= 2:
                try:
                    # Sum across all ALT alleles when multi-allelic kept; otherwise just the first ALT.
                    return sum(int(p) for p in parts[1:] if p not in (".", ""))
                except ValueError:
                    pass
        ao = s.get("AO")
        if ao and ao != ".":
            try:
                return sum(int(p) for p in ao.split(",") if p not in (".", ""))
            except ValueError:
                pass
        return None

    def ref_count(s: Dict[str, str]) -> Optional[int]:
        ad = s.get("AD")
        if ad and ad != ".":
            parts = ad.split(",")
            if parts and parts[0] not in (".", ""):
                try:
                    return int(parts[0])
                except ValueError:
                    pass
        ro = s.get("RO")
        return first_int(ro)

    def depth(s: Dict[str, str], ref_n: Optional[int], alt_n: Optional[int]) -> Optional[int]:
        dp = first_int(s.get("DP"))
        if dp is not None:
            return dp
        if ref_n is not None and alt_n is not None:
            return ref_n + alt_n
        return None

    t_alt = alt_count(t_s)
    t_ref = ref_count(t_s)
    n_alt = alt_count(n_s)
    n_ref = ref_count(n_s)
    t_dp = depth(t_s, t_ref, t_alt)
    n_dp = depth(n_s, n_ref, n_alt)

    if t_dp is None or n_dp is None or t_alt is None or n_alt is None:
        counters["parse_error"] += 1
        return False

    # 6. Depth bounds
    if t_dp < args.min_tumor_dp:
        counters["low_tumor_dp"] += 1
        return False
    if t_dp > args.max_tumor_dp:
        counters["high_tumor_dp"] += 1
        return False
    if n_dp < args.min_normal_dp:
        counters["low_normal_dp"] += 1
        return False
    if n_dp > args.max_normal_dp:
        counters["high_normal_dp"] += 1
        return False

    # 7. ALT-supporting depth
    if t_alt < args.min_tumor_ad:
        counters["low_tumor_ad"] += 1
        return False
    if n_alt > args.max_normal_ad:
        counters["high_normal_ad"] += 1
        return False

    # 8. VAF
    t_vaf = (t_alt / t_dp) if t_dp > 0 else 0.0
    n_vaf = (n_alt / n_dp) if n_dp > 0 else 0.0
    if t_vaf < args.min_tumor_vaf:
        counters["low_tumor_vaf"] += 1
        return False
    if n_vaf > args.max_normal_vaf:
        counters["high_normal_vaf"] += 1
        return False

    # 9-14. INFO-based optional filters
    needs_info = any(v is not None for v in (
        args.min_mqm, args.min_mqmr, args.max_sap,
        args.max_srp, args.max_epp, args.min_paired
    ))
    if needs_info:
        info = parse_info(info_str)

        def ifloat(k: str) -> Optional[float]:
            v = info.get(k)
            if v is None:
                return None
            return to_float(v.split(",")[0])

        if args.min_mqm is not None:
            mqm = ifloat("MQM")
            if mqm is None or mqm < args.min_mqm:
                counters["low_mqm"] += 1
                return False
        if args.min_mqmr is not None:
            mqmr = ifloat("MQMR")
            if mqmr is None or mqmr < args.min_mqmr:
                counters["low_mqmr"] += 1
                return False
        if args.max_sap is not None:
            sap = ifloat("SAP")
            if sap is not None and sap > args.max_sap:
                counters["high_sap"] += 1
                return False
        if args.max_srp is not None:
            srp = ifloat("SRP")
            if srp is not None and srp > args.max_srp:
                counters["high_srp"] += 1
                return False
        if args.max_epp is not None:
            epp = ifloat("EPP")
            if epp is not None and epp > args.max_epp:
                counters["high_epp"] += 1
                return False
        if args.min_paired is not None:
            paired = ifloat("PAIRED")
            if paired is None or paired < args.min_paired:
                counters["low_paired"] += 1
                return False

    return True


# ---------------------------------------------------------------------------
# Main per-file processing
# ---------------------------------------------------------------------------

def cutoff_lines(args: argparse.Namespace, tumor_name: str, normal_name: str) -> List[str]:
    """Generate ##sarek_subset_* header lines documenting the run."""
    stamp = _dt.datetime.now().isoformat(timespec="seconds")
    cmd = " ".join(shlex.quote(a) for a in sys.argv)
    items = [
        ("timestamp", stamp),
        ("command", cmd),
        ("tumor_sample", tumor_name),
        ("normal_sample", normal_name),
        ("qual", args.qual),
        ("min_tumor_vaf", args.min_tumor_vaf),
        ("min_tumor_ad", args.min_tumor_ad),
        ("min_tumor_dp", args.min_tumor_dp),
        ("max_tumor_dp", args.max_tumor_dp),
        ("max_normal_vaf", args.max_normal_vaf),
        ("max_normal_ad", args.max_normal_ad),
        ("min_normal_dp", args.min_normal_dp),
        ("max_normal_dp", args.max_normal_dp),
        ("pass_values", ",".join(sorted(args.pass_values_set))),
        ("keep_multiallelic", args.keep_multiallelic),
        ("snvs_only", args.snvs_only),
        ("min_mqm", args.min_mqm),
        ("min_mqmr", args.min_mqmr),
        ("max_sap", args.max_sap),
        ("max_srp", args.max_srp),
        ("max_epp", args.max_epp),
        ("min_paired", args.min_paired),
    ]
    return [f"##sarek_subset_{k}={v}\n" for k, v in items]


def process_file(input_path: str, out_path: str, args: argparse.Namespace,
                 log_every: int = 1_000_000) -> Dict[str, object]:
    counters = {r: 0 for r in REJECT_REASONS}
    n_in = 0
    n_out = 0
    t_start = time.time()
    tumor_name = normal_name = None
    t_idx = n_idx = -1

    tumor_re = re.compile(args.tumor_regex)
    normal_re = re.compile(args.normal_regex)

    with smart_open_read(input_path) as fh, \
         smart_open_write(out_path, args.gzip_output) as out:
        for line in fh:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                t_idx, n_idx, tumor_name, normal_name = detect_sample_indices(
                    line, tumor_re, normal_re, args.tumor_sample, args.normal_sample
                )
                # inject our annotations just before #CHROM
                out.write(f"##sarek_subset_input={os.path.basename(input_path)}\n")
                for hl in cutoff_lines(args, tumor_name, normal_name):
                    out.write(hl)
                out.write(line)
                continue

            n_in += 1
            fields = line.rstrip("\n").split("\t")
            if evaluate_record(fields, t_idx, n_idx, args, counters):
                out.write(line)
                n_out += 1

            if args.verbose and n_in % log_every == 0:
                elapsed = time.time() - t_start
                rate = n_in / elapsed if elapsed > 0 else 0.0
                sys.stderr.write(
                    f"  [{os.path.basename(input_path)}] read={n_in:,} kept={n_out:,} "
                    f"({rate:,.0f} rec/s)\n"
                )

    elapsed = time.time() - t_start
    return {
        "input": input_path,
        "output": out_path,
        "tumor_sample": tumor_name,
        "normal_sample": normal_name,
        "records_in": n_in,
        "records_out": n_out,
        "kept_fraction": (n_out / n_in) if n_in else 0.0,
        "elapsed_sec": round(elapsed, 1),
        **counters,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="freebayes_somatic_subset.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument(
        "inputs", nargs="*",
        help="One or more FreeBayes VCFs (.vcf or .vcf.gz). "
             "If empty, --input-dir is required.",
    )
    p.add_argument(
        "--input-dir",
        help="Directory of VCFs to process. Combined with --glob.",
    )
    p.add_argument(
        "--glob", default="*.vcf*",
        help="Glob (relative to --input-dir) for picking inputs. Default: %(default)s",
    )
    p.add_argument(
        "-o", "--output-dir",
        help="Where to write filtered VCFs. Created if needed. "
             "Mutually exclusive with --in-place.",
    )
    p.add_argument(
        "--in-place", action="store_true",
        help="Write output VCF next to each input instead of into --output-dir.",
    )
    p.add_argument(
        "--gzip-output", action="store_true",
        help="Write outputs as .vcf.gz (gzip; not bgzip-indexed).",
    )

    # Sample identification
    p.add_argument(
        "--tumor-regex", default=r"(?i)tumor|_T$|-T$|\.T$",
        help="Regex to match tumor sample column in #CHROM line. Default: %(default)s",
    )
    p.add_argument(
        "--normal-regex", default=r"(?i)normal|_N$|-N$|\.N$",
        help="Regex to match normal sample column. Default: %(default)s",
    )
    p.add_argument(
        "--tumor-sample",
        help="Force exact tumor sample name (overrides --tumor-regex).",
    )
    p.add_argument(
        "--normal-sample",
        help="Force exact normal sample name (overrides --normal-regex).",
    )

    # Core cutoffs
    p.add_argument("--qual", type=float, default=20.0,
                   help="Minimum site QUAL. Default: %(default)s")
    p.add_argument("--min-tumor-vaf", type=float, default=0.05,
                   help="Min tumor ALT VAF. Default: %(default)s")
    p.add_argument("--min-tumor-ad", type=int, default=4,
                   help="Min tumor ALT-supporting reads. Default: %(default)s")
    p.add_argument("--min-tumor-dp", type=int, default=10,
                   help="Min tumor total depth. Default: %(default)s")
    p.add_argument("--max-tumor-dp", type=int, default=1000,
                   help="Max tumor total depth (catch pileups). Default: %(default)s")
    p.add_argument("--max-normal-vaf", type=float, default=0.02,
                   help="Max normal ALT VAF. Default: %(default)s")
    p.add_argument("--max-normal-ad", type=int, default=1,
                   help="Max normal ALT-supporting reads. Default: %(default)s")
    p.add_argument("--min-normal-dp", type=int, default=10,
                   help="Min normal total depth. Default: %(default)s")
    p.add_argument("--max-normal-dp", type=int, default=1000,
                   help="Max normal total depth. Default: %(default)s")

    p.add_argument(
        "--pass-values", default="PASS,.",
        help="Comma-separated FILTER values to accept. Default: %(default)s",
    )
    p.add_argument(
        "--keep-multiallelic", action="store_true",
        help="Keep multi-ALT sites (default: drop them).",
    )
    p.add_argument(
        "--snvs-only", action="store_true",
        help="Restrict to SNVs (REF and every ALT have length 1).",
    )

    # Optional INFO-based filters (off unless set)
    p.add_argument("--min-mqm", type=float, default=None,
                   help="If set, require INFO/MQM (mean ALT MQ) >= this.")
    p.add_argument("--min-mqmr", type=float, default=None,
                   help="If set, require INFO/MQMR (mean REF MQ) >= this.")
    p.add_argument("--max-sap", type=float, default=None,
                   help="If set, require INFO/SAP (Phred ALT strand bias) <= this.")
    p.add_argument("--max-srp", type=float, default=None,
                   help="If set, require INFO/SRP (Phred REF strand bias) <= this.")
    p.add_argument("--max-epp", type=float, default=None,
                   help="If set, require INFO/EPP (Phred ALT end-placement bias) <= this.")
    p.add_argument("--min-paired", type=float, default=None,
                   help="If set, require INFO/PAIRED (proper-pair fraction) >= this.")

    p.add_argument("--summary-name", default="somatic_subset_summary.tsv",
                   help="Summary TSV filename written into --output-dir. Default: %(default)s")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Print progress to stderr.")
    p.add_argument("--dry-run", action="store_true",
                   help="Resolve inputs/outputs and print the plan; do not filter.")
    return p


def collect_inputs(args: argparse.Namespace) -> List[str]:
    inputs: List[str] = list(args.inputs)
    if args.input_dir:
        pattern = os.path.join(args.input_dir, args.glob)
        inputs.extend(sorted(glob.glob(pattern)))
    # Filter to vcf-ish only and dedupe
    seen = set()
    cleaned = []
    for p in inputs:
        if p in seen:
            continue
        seen.add(p)
        if not (p.endswith(".vcf") or p.endswith(".vcf.gz")):
            continue
        if not os.path.isfile(p):
            continue
        cleaned.append(p)
    return cleaned


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.in_place and not args.output_dir:
        sys.exit("ERROR: provide --output-dir DIR or pass --in-place.")
    if args.in_place and args.output_dir:
        sys.exit("ERROR: --in-place and --output-dir are mutually exclusive.")

    args.pass_values_set = {v.strip() for v in args.pass_values.split(",") if v.strip()}

    inputs = collect_inputs(args)
    if not inputs:
        sys.exit("ERROR: no input VCFs found. Pass paths or --input-dir.")

    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)

    print(f"# {len(inputs)} input VCF(s):", file=sys.stderr)
    for p in inputs:
        print(f"#   {p}", file=sys.stderr)

    if args.dry_run:
        for p in inputs:
            print(f"DRY-RUN -> {out_path_for(p, args.output_dir, args.in_place, args.gzip_output)}")
        return 0

    rows: List[Dict[str, object]] = []
    for ip in inputs:
        op = out_path_for(ip, args.output_dir, args.in_place, args.gzip_output)
        if args.verbose:
            sys.stderr.write(f"--> {ip}\n    {op}\n")
        try:
            row = process_file(ip, op, args)
        except Exception as e:
            sys.stderr.write(f"FAILED on {ip}: {e}\n")
            row = {
                "input": ip, "output": op, "tumor_sample": None, "normal_sample": None,
                "records_in": 0, "records_out": 0, "kept_fraction": 0.0,
                "elapsed_sec": 0.0, "error": str(e),
                **{r: 0 for r in REJECT_REASONS},
            }
        rows.append(row)
        sys.stderr.write(
            f"   kept {row['records_out']:,} / {row['records_in']:,} "
            f"({100*row['kept_fraction']:.2f}%) in {row['elapsed_sec']}s\n"
        )

    # Write summary TSV
    if args.output_dir:
        summary_path = os.path.join(args.output_dir, args.summary_name)
    else:
        summary_path = os.path.join(
            os.path.dirname(inputs[0]) or ".", args.summary_name
        )
    fieldnames = (
        ["input", "output", "tumor_sample", "normal_sample",
         "records_in", "records_out", "kept_fraction", "elapsed_sec"]
        + REJECT_REASONS
        + ["error"]
    )
    with open(summary_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    sys.stderr.write(f"\nSummary: {summary_path}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
