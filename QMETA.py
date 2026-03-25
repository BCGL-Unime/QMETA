# -*- coding: utf-8 -*-

"""

QIIME2 Interactive Metagenomics Pipeline - Fixed version

Fix applied: taxa collapse uses --o-collapsed-table instead of --output-dir
to avoid errors when the output directory already exists.

"""

import os
import re
import sys
import shutil
import subprocess
from datetime import datetime
from typing import Optional, List, Dict, Tuple

# -----------------------------
# Utilities and Validators
# -----------------------------

def log(msg: str, logfile: Optional[str] = None):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    if logfile:
        with open(logfile, "a", encoding="utf-8") as f:
            f.write(line + "\n")

def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

def run_cmd(cmd: List[str], logfile: Optional[str] = None) -> Tuple[int, str]:
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
        out = p.stdout or ""
        if logfile:
            log("Running: " + " ".join(cmd), logfile)
            if p.returncode != 0:
                # Log full output only on failure to avoid writing MBs for long-running QIIME2 commands
                log(out, logfile)
        return p.returncode, out
    except Exception as e:
        if logfile:
            log(f"ERROR executing command: {e}", logfile)
        return 1, str(e)

def ensure_int(prompt: str, default: Optional[int] = None, min_v: Optional[int] = None, max_v: Optional[int] = None, retries: int = 3) -> int:
    for _ in range(retries):
        raw = input(f"{prompt}{' ['+str(default)+']' if default is not None else ''}: ").strip()
        if raw.lower() in ("exit", "quit", "q"):
            print("Exit requested. Safe shutdown.")
            sys.exit(0)
        if raw == "" and default is not None:
            val = default
        else:
            if not re.match(r"^-?\d+$", raw):
                print("Invalid value: please enter an integer.")
                continue
            val = int(raw)
            if min_v is not None and val < min_v:
                print(f"Value too small (min {min_v}). Try again.")
                continue
            if max_v is not None and val > max_v:
                print(f"Value too large (max {max_v}). Try again.")
                continue
        return val
    print("Too many failed attempts. Safe shutdown.")
    sys.exit(1)

def ensure_yes_no(prompt: str, default: Optional[bool] = None) -> bool:
    suf = " [Y/n]" if default is True else (" [y/N]" if default is False else " [y/n]")
    while True:
        raw = input(f"{prompt}{suf}: ").strip().lower()
        if raw in ("exit", "quit", "q"):
            print("Exit requested. Safe shutdown.")
            sys.exit(0)
        if raw == "" and default is not None:
            return default
        if raw in ("y", "yes"):
            return True
        if raw in ("n", "no"):
            return False
        print("Please answer with 'y' or 'n'.")

def ensure_path(prompt: str, must_exist: bool = True, is_dir: Optional[bool] = None, default: Optional[str] = None) -> str:
    while True:
        raw = input(f"{prompt}{' ['+default+']' if default else ''}: ").strip()
        if raw.lower() in ("exit", "quit", "q"):
            print("Exit requested. Safe shutdown.")
            sys.exit(0)
        if raw == "" and default:
            raw = default
        path = os.path.abspath(raw)
        if must_exist:
            if not os.path.exists(path):
                print("Path does not exist. Try again.")
                continue
            if is_dir is True and not os.path.isdir(path):
                print("Expected a directory. Try again.")
                continue
            if is_dir is False and not os.path.isfile(path):
                print("Expected a file. Try again.")
                continue
        else:
            parent = os.path.dirname(path)
            if parent and not os.path.exists(parent):
                try:
                    os.makedirs(parent, exist_ok=True)
                except Exception as e:
                    print(f"Cannot create directory: {e}")
                    continue
        return path

def ensure_csv_list(prompt: str, allow_empty: bool = True) -> List[str]:
    raw = input(f"{prompt}").strip()
    if raw.lower() in ("exit", "quit", "q"):
        print("Exit requested. Safe shutdown.")
        sys.exit(0)
    if raw == "" and allow_empty:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]

TAXONOMY_LEVELS = [2, 5, 6]  # Taxonomic levels used for collapse, diversity, and differential analysis

def validate_reference_levels(refs_input: str, md_col: str) -> Optional[List[str]]:
    """Validate that reference levels follow the format column::value."""
    refs = [r.strip() for r in refs_input.split(",") if r.strip()]
    for r in refs:
        if "::" not in r:
            print(f"  Invalid format '{r}'. Expected '{md_col}::value'.")
            return None
        col, val = r.split("::", 1)
        if col != md_col:
            print(f"  Column mismatch: got '{col}', expected '{md_col}'.")
            return None
        if not val:
            print(f"  Empty value in '{r}'.")
            return None
    return refs

# -----------------------------
# Pipeline State
# -----------------------------

class PipelineState:
    def __init__(self):
        self.params: Dict[str, object] = {}
        self.status_msgs: List[str] = []
        self.errors: List[str] = []
        self.outdir: str = os.path.abspath("qiime2_out")
        self.logfile: str = os.path.join(self.outdir, "pipeline.log")
        os.makedirs(self.outdir, exist_ok=True)

    def add_status(self, msg: str):
        self.status_msgs.append(msg)
        log(msg, self.logfile)

    def add_error(self, msg: str):
        self.errors.append(msg)
        log("ERROR: " + msg, self.logfile)

# -----------------------------
# Pipeline
# -----------------------------

class QIIME2Pipeline:
    def __init__(self, state: PipelineState):
        self.st = state
        if not which("qiime"):
            self.st.add_error("QIIME 2 CLI not found in PATH.")
            print("Error: QIIME 2 not found. Make sure 'qiime' is in PATH.")
            sys.exit(1)

    def maybe_resume(self, targets: List[str], step_name: str) -> bool:
        """Check if output files exist first, then ask whether to reuse them."""
        all_exist = all(t and os.path.exists(t) for t in targets)

        if all_exist:
            reuse = ensure_yes_no(f"Existing outputs found for {step_name}. Reuse them without re-running?", default=True)
            if reuse:
                return True
            else:
                self.st.add_status(f"Recalculation requested for {step_name}.")
                return False
        else:
            self.st.add_status(f"Outputs for {step_name} not found: calculation needed.")
            return False

    # ---- STEP 1: Import ----
    def step_import(self):
        print("\n=== STEP 1/8: Data import ===")
        paired = ensure_yes_no("Is the data PAIRED-END?", default=True)
        self.st.params["paired_end"] = paired

        if paired:
            print("Paired-end import using MANIFEST V2 (PairedEndFastqManifestPhred33V2).")
            manifest = ensure_path("Path to existing TSV manifest (V2)", must_exist=True, is_dir=False, default="manifest.tsv")
            out_qza = ensure_path("Output path for imported reads .qza", must_exist=False, is_dir=False, default=os.path.join(self.st.outdir, "import.qza"))

            if self.maybe_resume([out_qza], "Import (paired)"):
                self.st.params["demux_qza"] = out_qza
                self.st.add_status(f"Paired-end import: reusing {out_qza}")
            else:
                cmd = [
                    "qiime", "tools", "import",
                    "--type", "SampleData[PairedEndSequencesWithQuality]",
                    "--input-path", manifest,
                    "--output-path", out_qza,
                    "--input-format", "PairedEndFastqManifestPhred33V2"
                ]
                code, _ = run_cmd(cmd, self.st.logfile)
                if code != 0:
                    self.st.add_error("Paired-end import failed. Check the V2 manifest.")
                    self.validate_and_continue(False)
                    return
                self.st.params["demux_qza"] = out_qza
                self.st.add_status(f"Paired-end import complete: {out_qza}")
        else:
            print("Single-end import using MANIFEST V2.")
            manifest = ensure_path("Path to existing TSV manifest", must_exist=True, is_dir=False, default="manifest.tsv")
            out_qza = ensure_path("Output path for imported reads .qza", must_exist=False, is_dir=False, default=os.path.join(self.st.outdir, "import.qza"))

            if self.maybe_resume([out_qza], "Import (single)"):
                self.st.params["demux_qza"] = out_qza
                self.st.add_status(f"Single-end import: reusing {out_qza}")
            else:
                cmd = [
                    "qiime", "tools", "import",
                    "--type", "SampleData[SequencesWithQuality]",
                    "--input-path", manifest,
                    "--output-path", out_qza,
                    "--input-format", "SingleEndFastqManifestPhred33V2"
                ]
                code, _ = run_cmd(cmd, self.st.logfile)
                if code != 0:
                    self.st.add_error("Single-end import failed.")
                    self.validate_and_continue(False)
                    return
                self.st.params["demux_qza"] = out_qza
                self.st.add_status(f"Single-end import complete: {out_qza}")

        demux_qzv_import = os.path.join(self.st.outdir, "summary_import.qzv")
        if not self.maybe_resume([demux_qzv_import], "demux summarize (import)"):
            code_sum, _ = run_cmd(["qiime", "demux", "summarize", "--i-data", self.st.params["demux_qza"], "--o-visualization", demux_qzv_import], self.st.logfile)
            if code_sum == 0:
                self.st.add_status(f"Import quality plot generated: {demux_qzv_import}")
            else:
                self.st.add_error("Failed to generate summary_import.qzv.")
        else:
            self.st.add_status(f"Reusing import quality plot: {demux_qzv_import}")

        self.validate_and_continue(True)

    # ---- STEP 2: Cutadapt (optional) ----
    def step_cutadapt(self):
        print("\n=== STEP 2/8: Trimming (optional) with cutadapt ===")
        do_trim = ensure_yes_no("Run trimming/primer removal with cutadapt?", default=False)
        self.st.params["trimming"] = "Yes" if do_trim else "No"
        input_qza = self.st.params.get("demux_qza")

        if not do_trim:
            self.st.add_status("Trimming skipped.")
            self.st.params["trimmed_qza"] = input_qza
            self.validate_and_continue(True)
            return

        out_qza = os.path.join(self.st.outdir, "import_cutadapt.qza")
        if self.maybe_resume([out_qza], "cutadapt"):
            self.st.params["trimmed_qza"] = out_qza
            self.st.add_status(f"Reusing trimming output: {out_qza}")
        else:
            if self.st.params.get("paired_end", True):
                primer_f = input("Enter forward primer (sequence): ").strip()
                primer_r = input("Enter reverse primer (sequence): ").strip()
                n_cores = ensure_int("Number of cores for cutadapt", default=30, min_v=1)
                self.st.params["primer_f"] = primer_f
                self.st.params["primer_r"] = primer_r

                cmd = [
                    "qiime", "cutadapt", "trim-paired",
                    "--i-demultiplexed-sequences", input_qza,
                    "--p-cores", str(n_cores),
                    "--p-front-f", primer_f,
                    "--p-front-r", primer_r,
                    "--p-match-adapter-wildcards",
                    "--p-discard-untrimmed",
                    "--o-trimmed-sequences", out_qza
                ]
            else:
                primer = input("Enter primer (single-end, sequence): ").strip()
                n_cores = ensure_int("Number of cores for cutadapt", default=30, min_v=1)
                self.st.params["primer_f"] = primer

                cmd = [
                    "qiime", "cutadapt", "trim-single",
                    "--i-demultiplexed-sequences", input_qza,
                    "--p-cores", str(n_cores),
                    "--p-front", primer,
                    "--p-match-adapter-wildcards",
                    "--p-discard-untrimmed",
                    "--o-trimmed-sequences", out_qza
                ]

            code, _ = run_cmd(cmd, self.st.logfile)
            if code != 0:
                self.st.add_error("Cutadapt trimming failed. Using untrimmed input as fallback.")
                self.st.params["trimmed_qza"] = input_qza
            else:
                self.st.params["trimmed_qza"] = out_qza
                self.st.add_status(f"Trimming complete. Output: {out_qza}")

        demux_qzv = os.path.join(self.st.outdir, "summary_import_cutadapt.qzv")
        if not self.maybe_resume([demux_qzv], "demux summarize (post-trim)"):
            code2, _ = run_cmd(["qiime", "demux", "summarize", "--i-data", self.st.params["trimmed_qza"], "--o-visualization", demux_qzv], self.st.logfile)
            if code2 == 0:
                self.st.add_status(f"Post-trim quality plot generated: {demux_qzv}")
            else:
                self.st.add_error("Failed to generate summary_import_cutadapt.qzv.")
        else:
            self.st.add_status(f"Reusing post-trim quality plot: {demux_qzv}")

        self.validate_and_continue(True)

    # ---- STEP 3: DADA2 ----
    def step_dada2(self):
        print("\n=== STEP 3/8: Denoising with DADA2 ===")
        demux_for_dada2 = self.st.params.get("trimmed_qza") or self.st.params.get("demux_qza")
        paired = self.st.params.get("paired_end", True)

        trim_left_f = ensure_int("Enter --p-trim-left-f (integer, >=0)", default=0, min_v=0)
        trim_left_r = ensure_int("Enter --p-trim-left-r (integer, >=0)", default=0, min_v=0) if paired else None
        trunc_len_f = ensure_int("Enter --p-trunc-len-f (integer, >=0; 0=no truncation)", default=248, min_v=0)
        trunc_len_r = ensure_int("Enter --p-trunc-len-r (integer, >=0; 0=no truncation)", default=248, min_v=0) if paired else None
        trunc_q = ensure_int("Enter --p-trunc-q (quality threshold)", default=2, min_v=0)
        min_fold = ensure_int("Enter --p-min-fold-parent-over-abundance", default=4, min_v=1)
        n_threads = ensure_int("Number of threads for DADA2 (0 = all cores)", default=0, min_v=0)

        table_qza = os.path.join(self.st.outdir, "table.qza")
        repseqs_qza = os.path.join(self.st.outdir, "rep-seqs.qza")
        stats_qza = os.path.join(self.st.outdir, "denoising-stats.qza")

        if not self.maybe_resume([table_qza, repseqs_qza, stats_qza], "DADA2"):
            if paired:
                cmd = [
                    "qiime", "dada2", "denoise-paired",
                    "--i-demultiplexed-seqs", demux_for_dada2,
                    "--p-trim-left-f", str(trim_left_f),
                    "--p-trim-left-r", str(trim_left_r),
                    "--p-trunc-len-f", str(trunc_len_f),
                    "--p-trunc-len-r", str(trunc_len_r),
                    "--p-trunc-q", str(trunc_q),
                    "--p-min-fold-parent-over-abundance", str(min_fold),
                    "--p-n-threads", str(n_threads),
                    "--o-table", table_qza,
                    "--o-representative-sequences", repseqs_qza,
                    "--o-denoising-stats", stats_qza
                ]
            else:
                cmd = [
                    "qiime", "dada2", "denoise-single",
                    "--i-demultiplexed-seqs", demux_for_dada2,
                    "--p-trim-left", str(trim_left_f),
                    "--p-trunc-len", str(trunc_len_f),
                    "--p-trunc-q", str(trunc_q),
                    "--p-min-fold-parent-over-abundance", str(min_fold),
                    "--p-n-threads", str(n_threads),
                    "--o-table", table_qza,
                    "--o-representative-sequences", repseqs_qza,
                    "--o-denoising-stats", stats_qza
                ]
            code, _ = run_cmd(cmd, self.st.logfile)
            if code != 0:
                self.st.add_error("DADA2 denoising failed. Review the parameters.")
                self.validate_and_continue(False)
                return
            self.st.add_status(f"DADA2 complete (n-threads={n_threads}, trunc-q={trunc_q}, min-fold={min_fold}). ASV table: {table_qza}")
        else:
            self.st.add_status(f"Reusing DADA2 results: {table_qza}, {repseqs_qza}")

        md_file = ensure_path("Path to QIIME2 metadata TSV/TXT file", must_exist=True, is_dir=False, default="metadata.txt")
        self.st.params["metadata_file"] = md_file

        table_qzv = os.path.join(self.st.outdir, "table.qzv")
        if not self.maybe_resume([table_qzv], "feature-table summarize"):
            cmd_tbl = ["qiime", "feature-table", "summarize", "--i-table", table_qza, "--o-visualization", table_qzv, "--m-sample-metadata-file", md_file]
            code_tbl, _ = run_cmd(cmd_tbl, self.st.logfile)
            if code_tbl == 0:
                self.st.add_status(f"Feature table summary generated: {table_qzv}")
            else:
                self.st.add_error("Failed to generate table.qzv.")
        else:
            self.st.add_status(f"Reusing table.qzv: {table_qzv}")

        stats_qzv = os.path.join(self.st.outdir, "denoising-stats.qzv")
        if not self.maybe_resume([stats_qzv], "metadata tabulate (denoising-stats)"):
            code_stats, _ = run_cmd(["qiime", "metadata", "tabulate", "--m-input-file", stats_qza, "--o-visualization", stats_qzv], self.st.logfile)
            if code_stats == 0:
                self.st.add_status(f"Denoising stats visualization generated: {stats_qzv}")
            else:
                self.st.add_error("Failed to generate denoising-stats.qzv.")
        else:
            self.st.add_status(f"Reusing denoising-stats.qzv: {stats_qzv}")

        self.st.params["table_qza"] = table_qza
        self.st.params["repseqs_qza"] = repseqs_qza

        print("\n" + "="*80)
        print("IMPORTANT: Check table.qzv to verify read counts per sample.")
        print("="*80)

        self.validate_and_continue(True)

    # ---- STEP 4: Taxonomy ----
    def step_taxonomy(self):
        print("\n=== STEP 4/8: Taxonomic classification ===")
        table_qza = self.st.params.get("table_qza")
        repseqs_qza = self.st.params.get("repseqs_qza")

        if not (table_qza and repseqs_qza):
            self.st.add_error("Missing table/rep-seqs for taxonomy classification.")
            self.validate_and_continue(False)
            return

        classif_dir = os.path.join(self.st.outdir, "classification")
        os.makedirs(classif_dir, exist_ok=True)
        taxonomy_qza = os.path.join(classif_dir, "taxonomy_denoised.qza")

        if self.maybe_resume([taxonomy_qza], "taxonomy classification"):
            self.st.add_status(f"Reusing taxonomy: {taxonomy_qza}")
            self.st.params["taxonomy_qza"] = taxonomy_qza
            self.validate_and_continue(True)
            return

        have_tax = ensure_yes_no("Do you already have a taxonomy.qza to use?", default=False)
        if have_tax:
            taxonomy_qza = ensure_path("Path to taxonomy.qza", must_exist=True, is_dir=False)
            self.st.params["taxonomy_qza"] = taxonomy_qza
            self.st.add_status(f"Using provided taxonomy: {taxonomy_qza}")
            self.validate_and_continue(True)
            return

        do_classify = ensure_yes_no("Classify now using a pre-trained classifier.qza?", default=True)
        if not do_classify:
            self.st.add_error("No taxonomy provided. Taxonomic filtering will not be applicable.")
            self.validate_and_continue(False)
            return

        classifier_qza = ensure_path("Path to pre-trained classifier.qza", must_exist=True, is_dir=False, default="/mnt/hdd2/database/qiime_pretrained_2024-5/silva-138-99-nb-classifier.qza")
        n_jobs = ensure_int("Number of jobs for classify-sklearn", default=30, min_v=1)

        cmd = [
            "qiime", "feature-classifier", "classify-sklearn",
            "--i-classifier", classifier_qza,
            "--i-reads", repseqs_qza,
            "--o-classification", taxonomy_qza,
            "--p-n-jobs", str(n_jobs)
        ]
        code, _ = run_cmd(cmd, self.st.logfile)
        if code != 0:
            self.st.add_error("Taxonomic classification failed.")
            self.validate_and_continue(False)
        else:
            self.st.params["taxonomy_qza"] = taxonomy_qza
            self.st.add_status(f"Taxonomy generated: {taxonomy_qza}")
            self.validate_and_continue(True)

    # ---- STEP 5: Taxonomic filter and visualizations ----
    def step_taxa_filter(self):
        print("\n=== STEP 5/8: Taxonomic filter and visualizations ===")
        table_qza = self.st.params.get("table_qza")
        repseqs_qza = self.st.params.get("repseqs_qza")
        taxonomy_qza = self.st.params.get("taxonomy_qza")
        md_file = self.st.params.get("metadata_file")

        if not (table_qza and repseqs_qza and taxonomy_qza and os.path.exists(taxonomy_qza)):
            self.st.add_error("Missing taxonomy: skipping taxonomic filter.")
            self.validate_and_continue(False)
            return

        # Exclude eukaryota and unassigned features (case-insensitive matching)
        exclude_terms = "eukaryota,unassigned"
        table_filt = os.path.join(self.st.outdir, "table-no-eukaryota-no-unassigned.qza")

        if not self.maybe_resume([table_filt], "taxa filter-table"):
            code1, _ = run_cmd([
                "qiime", "taxa", "filter-table",
                "--i-table", table_qza,
                "--i-taxonomy", taxonomy_qza,
                "--p-exclude", exclude_terms,
                "--o-filtered-table", table_filt
            ], self.st.logfile)
            if code1 != 0:
                self.st.add_error("Table filtering failed. Using unfiltered table.")
                table_filt = table_qza
            else:
                self.st.add_status(f"Filtered table generated: {table_filt}")
        else:
            self.st.add_status(f"Reusing filtered table: {table_filt}")

        table_filt_qzv = os.path.join(self.st.outdir, "table-no-eukaryota-no-unassigned.qzv")
        if not self.maybe_resume([table_filt_qzv], "feature-table summarize (filtered)"):
            cmd_tbl = ["qiime", "feature-table", "summarize", "--i-table", table_filt, "--o-visualization", table_filt_qzv]
            if md_file and os.path.exists(md_file):
                cmd_tbl += ["--m-sample-metadata-file", md_file]
            code_tbl, _ = run_cmd(cmd_tbl, self.st.logfile)
            if code_tbl == 0:
                self.st.add_status(f"Filtered feature table summary generated: {table_filt_qzv}")
            else:
                self.st.add_error("Failed to generate table-no-eukaryota-no-unassigned.qzv.")
        else:
            self.st.add_status(f"Reusing filtered table.qzv: {table_filt_qzv}")

        repseqs_qzv = os.path.join(self.st.outdir, "rep-seqs.qzv")
        if not self.maybe_resume([repseqs_qzv], "tabulate-seqs with taxonomy"):
            code_rep, _ = run_cmd([
                "qiime", "feature-table", "tabulate-seqs",
                "--i-data", repseqs_qza,
                "--i-taxonomy", taxonomy_qza,
                "--o-visualization", repseqs_qzv
            ], self.st.logfile)
            if code_rep == 0:
                self.st.add_status(f"Rep-seqs visualization (with taxonomy) generated: {repseqs_qzv}")
            else:
                self.st.add_error("Failed to generate rep-seqs.qzv.")
        else:
            self.st.add_status(f"Reusing rep-seqs.qzv: {repseqs_qzv}")

        classif_dir = os.path.join(self.st.outdir, "classification")
        taxonomy_metadata_qzv = os.path.join(classif_dir, "taxonomy_denoised_metadata.qzv")
        if not self.maybe_resume([taxonomy_metadata_qzv], "metadata tabulate (taxonomy)"):
            code_tax, _ = run_cmd([
                "qiime", "metadata", "tabulate",
                "--m-input-file", taxonomy_qza,
                "--o-visualization", taxonomy_metadata_qzv
            ], self.st.logfile)
            if code_tax == 0:
                self.st.add_status(f"Taxonomy metadata visualization generated: {taxonomy_metadata_qzv}")
            else:
                self.st.add_error("Failed to generate taxonomy_denoised_metadata.qzv.")
        else:
            self.st.add_status(f"Reusing taxonomy_denoised_metadata.qzv: {taxonomy_metadata_qzv}")

        taxa_barplot_qzv = os.path.join(classif_dir, "taxa-bar-plots.qzv")
        if not self.maybe_resume([taxa_barplot_qzv], "taxa barplot"):
            code_bar, _ = run_cmd([
                "qiime", "taxa", "barplot",
                "--i-table", table_filt,
                "--i-taxonomy", taxonomy_qza,
                "--m-metadata-file", md_file,
                "--o-visualization", taxa_barplot_qzv
            ], self.st.logfile)
            if code_bar == 0:
                self.st.add_status(f"Taxa barplot generated: {taxa_barplot_qzv}")
            else:
                self.st.add_error("Failed to generate taxa-bar-plots.qzv.")
        else:
            self.st.add_status(f"Reusing taxa barplot: {taxa_barplot_qzv}")

        self.st.params["table_for_downstream"] = table_filt

        print("\n" + "="*80)
        print("IMPORTANT: Check table-no-eukaryota-no-unassigned.qzv to determine sampling depth.")
        print("="*80)

        self.validate_and_continue(True)

    # ---- STEP 6: Taxa collapse and export ----
    def step_taxa_collapse(self):
        print("\n=== STEP 6/8: Taxa collapse and export (levels 2, 5, 6) ===")
        table_filt = self.st.params.get("table_for_downstream")
        taxonomy_qza = self.st.params.get("taxonomy_qza")
        classif_dir = os.path.join(self.st.outdir, "classification")

        if not (table_filt and taxonomy_qza):
            self.st.add_error("Missing filtered table or taxonomy for taxa collapse.")
            self.validate_and_continue(False)
            return

        levels = TAXONOMY_LEVELS

        for level in levels:
            level_dir = os.path.join(classif_dir, f"level{level}")
            os.makedirs(level_dir, exist_ok=True)

            collapsed_table = os.path.join(level_dir, "collapsed_table.qza")

            # FIX: use --o-collapsed-table instead of --output-dir
            if not self.maybe_resume([collapsed_table], f"taxa collapse level {level}"):
                code, _ = run_cmd([
                    "qiime", "taxa", "collapse",
                    "--i-table", table_filt,
                    "--i-taxonomy", taxonomy_qza,
                    "--p-level", str(level),
                    "--o-collapsed-table", collapsed_table
                ], self.st.logfile)
                if code != 0:
                    self.st.add_error(f"Taxa collapse failed for level {level}.")
                    continue
                self.st.add_status(f"Taxa collapse complete for level {level}: {collapsed_table}")
            else:
                self.st.add_status(f"Reusing taxa collapse level {level}: {collapsed_table}")

            relative_table = os.path.join(level_dir, f"level{level}_relative.qza")
            if not self.maybe_resume([relative_table], f"relative frequency level {level}"):
                code, _ = run_cmd([
                    "qiime", "feature-table", "relative-frequency",
                    "--i-table", collapsed_table,
                    "--o-relative-frequency-table", relative_table
                ], self.st.logfile)
                if code != 0:
                    self.st.add_error(f"Relative frequency failed for level {level}.")
                    continue
                self.st.add_status(f"Relative frequency generated for level {level}: {relative_table}")
            else:
                self.st.add_status(f"Reusing relative frequency level {level}: {relative_table}")

            relative_qzv = os.path.join(level_dir, f"level{level}_relative_tabulated.qzv")
            if not self.maybe_resume([relative_qzv], f"tabulate relative level {level}"):
                code, _ = run_cmd([
                    "qiime", "metadata", "tabulate",
                    "--m-input-file", relative_table,
                    "--o-visualization", relative_qzv
                ], self.st.logfile)
                if code != 0:
                    self.st.add_error(f"Tabulate relative frequency failed for level {level}.")
                else:
                    self.st.add_status(f"Relative frequency tabulated for level {level}: {relative_qzv}")
            else:
                self.st.add_status(f"Reusing relative frequency tabulated level {level}: {relative_qzv}")

            feature_table_biom = os.path.join(level_dir, "feature-table.biom")
            if not self.maybe_resume([feature_table_biom], f"export BIOM level {level}"):
                code, _ = run_cmd([
                    "qiime", "tools", "export",
                    "--input-path", collapsed_table,
                    "--output-path", level_dir
                ], self.st.logfile)
                if code != 0:
                    self.st.add_error(f"BIOM export failed for level {level}.")
                else:
                    self.st.add_status(f"BIOM export complete for level {level}: {feature_table_biom}")
            else:
                self.st.add_status(f"Reusing BIOM level {level}: {feature_table_biom}")

            tsv_file = os.path.join(level_dir, f"level{level}.tsv")
            if not self.maybe_resume([tsv_file], f"convert TSV level {level}"):
                code, _ = run_cmd([
                    "biom", "convert",
                    "-i", feature_table_biom,
                    "-o", tsv_file,
                    "--to-tsv"
                ], self.st.logfile)
                if code != 0:
                    self.st.add_error(f"TSV conversion failed for level {level}.")
                else:
                    self.st.add_status(f"TSV generated for level {level}: {tsv_file}")
            else:
                self.st.add_status(f"Reusing TSV level {level}: {tsv_file}")

        self.validate_and_continue(True)

    # ---- STEP 7: Diversity (core-metrics non-phylogenetic + alpha on collapsed tables) ----
    def step_diversity(self):
        print("\n=== STEP 7/8: Diversity analysis ===")
        table_filt = self.st.params.get("table_for_downstream")
        md_file = self.st.params.get("metadata_file")
        classif_dir = os.path.join(self.st.outdir, "classification")

        if not (table_filt and md_file):
            self.st.add_error("Missing filtered table or metadata for diversity analysis.")
            self.validate_and_continue(False)
            return

        sampling_depth = ensure_int("Enter --p-sampling-depth for core-metrics (integer > 0)", default=68037, min_v=1)
        self.st.params["p-sampling-depth"] = sampling_depth

        core_dir = os.path.join(self.st.outdir, "core-metrics-results")
        if self.maybe_resume([core_dir], "core-metrics") and os.path.isdir(core_dir):
            self.st.add_status(f"Reusing core-metrics: {core_dir}")
        else:
            cmd_core = [
                "qiime", "diversity", "core-metrics",
                "--i-table", table_filt,
                "--p-sampling-depth", str(sampling_depth),
                "--m-metadata-file", md_file,
                "--output-dir", core_dir
            ]
            code, _ = run_cmd(cmd_core, self.st.logfile)
            if code != 0:
                self.st.add_error("Core metrics failed.")
            else:
                self.st.add_status(f"Core metrics (non-phylogenetic) complete. Output: {core_dir}")

        levels = TAXONOMY_LEVELS
        metrics = ["chao1", "simpson", "shannon", "ace", "observed_features", "pielou_e", "simpson_e"]

        for level in levels:
            level_dir = os.path.join(classif_dir, f"level{level}")
            collapsed_table = os.path.join(level_dir, "collapsed_table.qza")

            if not os.path.exists(collapsed_table):
                self.st.add_error(f"Collapsed table not found for level {level}: {collapsed_table}")
                continue

            alpha_output_dir = os.path.join(level_dir, "alpha-diversity")
            os.makedirs(alpha_output_dir, exist_ok=True)

            all_alpha_files = []
            for metric in metrics:
                all_alpha_files.append(os.path.join(alpha_output_dir, f"{metric}_vector.qza"))
                all_alpha_files.append(os.path.join(alpha_output_dir, f"{metric}_group-significance.qzv"))

            if self.maybe_resume(all_alpha_files, f"Alpha diversity level {level} (all metrics)"):
                self.st.add_status(f"Reusing all alpha diversity indices for level {level}.")
                continue

            for metric in metrics:
                alpha_vector_file = os.path.join(alpha_output_dir, f"{metric}_vector.qza")
                significance_output = os.path.join(alpha_output_dir, f"{metric}_group-significance.qzv")

                code, _ = run_cmd([
                    "qiime", "diversity", "alpha",
                    "--i-table", collapsed_table,
                    "--p-metric", metric,
                    "--o-alpha-diversity", alpha_vector_file
                ], self.st.logfile)
                if code != 0:
                    self.st.add_error(f"Alpha diversity failed for level {level}, metric {metric}.")
                    continue
                self.st.add_status(f"Alpha diversity calculated for level {level}, metric {metric}: {alpha_vector_file}")

                code, _ = run_cmd([
                    "qiime", "diversity", "alpha-group-significance",
                    "--i-alpha-diversity", alpha_vector_file,
                    "--m-metadata-file", md_file,
                    "--o-visualization", significance_output
                ], self.st.logfile)
                if code != 0:
                    self.st.add_error(f"Alpha group significance failed for level {level}, metric {metric}.")
                else:
                    self.st.add_status(f"Alpha group significance generated for level {level}, metric {metric}: {significance_output}")

        max_depth = ensure_int("Enter --p-max-depth for alpha-rarefaction (default 52724)", default=52724, min_v=sampling_depth)
        self.st.params["p-max-depth"] = max_depth

        raref_html = os.path.join(core_dir, "alpha-rarefaction.qzv")
        if not self.maybe_resume([raref_html], "alpha-rarefaction"):
            code, _ = run_cmd([
                "qiime", "diversity", "alpha-rarefaction",
                "--i-table", table_filt,
                "--m-metadata-file", md_file,
                "--p-max-depth", str(max_depth),
                "--o-visualization", raref_html
            ], self.st.logfile)
            if code != 0:
                self.st.add_error("Alpha rarefaction failed.")
            else:
                self.st.add_status(f"Alpha rarefaction complete: {raref_html}")
        else:
            self.st.add_status(f"Reusing alpha-rarefaction: {raref_html}")

        self.validate_and_continue(True)

    # ---- STEP 8: Differential analysis (ANCOM-BC) on collapsed tables ----
    def step_differential(self):
        print("\n=== STEP 8/8: Differential analysis (ANCOM-BC) ===")
        md_file = self.st.params.get("metadata_file")
        classif_dir = os.path.join(self.st.outdir, "classification")

        md_col = input("Enter metadata column for ANCOM-BC formula: ").strip()
        self.st.params["metadata_column"] = md_col

        levels = TAXONOMY_LEVELS

        # Collect per-level reference levels interactively before running any commands
        for level in levels:
            print(f"\n{'='*80}")
            print(f"ANCOM-BC for LEVEL {level}")
            print(f"{'='*80}")

            refs = None
            while refs is None:
                refs_input = input(f"Enter reference level(s) for level {level} (format: {md_col}::value, comma-separated, leave empty for none): ").strip()
                if refs_input.lower() in ("exit", "quit", "q"):
                    print("Exit requested. Safe shutdown.")
                    sys.exit(0)

                if not refs_input:
                    refs = []
                    print("Warning: no reference level specified for this level.")
                    break

                refs = validate_reference_levels(refs_input, md_col)

                if refs is None:
                    print("\nWarning: please re-enter reference levels correctly.\n")

            self.st.params[f"references_level{level}"] = refs

        # Run ANCOM-BC for each level using the per-level reference levels collected above
        for level in levels:
            level_dir = os.path.join(classif_dir, f"level{level}")
            collapsed_table = os.path.join(level_dir, "collapsed_table.qza")

            if not os.path.exists(collapsed_table):
                self.st.add_error(f"Collapsed table not found for level {level}: {collapsed_table}")
                continue

            da_output_dir = os.path.join(level_dir, "DA")
            os.makedirs(da_output_dir, exist_ok=True)

            level_refs = self.st.params.get(f"references_level{level}", [])

            if not level_refs:
                ancombc_qza = os.path.join(da_output_dir, "ancombc_differentials.qza")
                if not self.maybe_resume([ancombc_qza], f"ANCOM-BC level {level} (no reference)"):
                    cmd = [
                        "qiime", "composition", "ancombc",
                        "--i-table", collapsed_table,
                        "--m-metadata-file", md_file,
                        "--p-formula", md_col,
                        "--o-differentials", ancombc_qza
                    ]
                    code, _ = run_cmd(cmd, self.st.logfile)
                    if code != 0:
                        self.st.add_error(f"ANCOM-BC failed for level {level}. Check formula.")
                        continue
                    self.st.add_status(f"ANCOM-BC complete for level {level} (formula '{md_col}', no reference). Output: {ancombc_qza}")
                else:
                    self.st.add_status(f"Reusing ANCOM-BC level {level}: {ancombc_qza}")

                barplot_qzv = os.path.join(da_output_dir, "ancombc_differentials.qzv")
                if not self.maybe_resume([barplot_qzv], f"da-barplot level {level}"):
                    code_plot, _ = run_cmd([
                        "qiime", "composition", "da-barplot",
                        "--i-data", ancombc_qza,
                        "--p-significance-threshold", "0.05",
                        "--p-level-delimiter", ";",
                        "--o-visualization", barplot_qzv
                    ], self.st.logfile)
                    if code_plot != 0:
                        self.st.add_error(f"da-barplot generation failed for level {level}.")
                    else:
                        self.st.add_status(f"da-barplot generated for level {level}: {barplot_qzv}")
                else:
                    self.st.add_status(f"Reusing da-barplot level {level}: {barplot_qzv}")
            else:
                for r in level_refs:
                    r = r.strip()
                    ref_value = r.split("::")[-1] if "::" in r else r
                    ancombc_qza = os.path.join(da_output_dir, f"{md_col}_{ref_value}_ref_differentials.qza")

                    if not self.maybe_resume([ancombc_qza], f"ANCOM-BC level {level} (ref={ref_value})"):
                        cmd = [
                            "qiime", "composition", "ancombc",
                            "--i-table", collapsed_table,
                            "--m-metadata-file", md_file,
                            "--p-formula", md_col,
                            "--p-reference-levels", r,
                            "--o-differentials", ancombc_qza
                        ]
                        code, _ = run_cmd(cmd, self.st.logfile)
                        if code != 0:
                            self.st.add_error(f"ANCOM-BC failed for level {level}, reference '{ref_value}'. Check formula and reference levels.")
                            continue
                        self.st.add_status(f"ANCOM-BC complete for level {level} (formula '{md_col}', ref='{ref_value}'). Output: {ancombc_qza}")
                    else:
                        self.st.add_status(f"Reusing ANCOM-BC level {level} (ref={ref_value}): {ancombc_qza}")

                    barplot_qzv = os.path.join(da_output_dir, f"{md_col}_{ref_value}_ref_differentials.qzv")
                    if not self.maybe_resume([barplot_qzv], f"da-barplot level {level} (ref={ref_value})"):
                        code_plot, _ = run_cmd([
                            "qiime", "composition", "da-barplot",
                            "--i-data", ancombc_qza,
                            "--p-significance-threshold", "0.05",
                            "--p-level-delimiter", ";",
                            "--o-visualization", barplot_qzv
                        ], self.st.logfile)
                        if code_plot != 0:
                            self.st.add_error(f"da-barplot generation failed for level {level}, reference '{ref_value}'.")
                        else:
                            self.st.add_status(f"da-barplot generated for level {level}, reference '{ref_value}': {barplot_qzv}")
                    else:
                        self.st.add_status(f"Reusing da-barplot level {level} (ref={ref_value}): {barplot_qzv}")

        self.st.params["differential_method"] = "ANCOM-BC"
        self.validate_and_continue(True)

    # ---- Final summary ----
    def step_summary(self):
        print("\n=== FINAL SUMMARY ===")
        p = self.st.params
        trimming = p.get("trimming", "No")
        samp_depth = p.get("p-sampling-depth", "-")
        max_depth = p.get("p-max-depth", "-")
        md_col = f'"{p.get("metadata_column")}"' if p.get("metadata_column") else "-"
        refs_by_level = {
            level: p.get(f"references_level{level}", [])
            for level in [2, 5, 6]
        }
        refs_str = "; ".join(
            f"L{level}: {', '.join(refs) if refs else 'none'}"
            for level, refs in refs_by_level.items()
        )

        print("\n---")
        print("CHOSEN PARAMETERS:")
        print(f"- Trimming: {trimming}")
        print(f"- p-sampling-depth: {samp_depth}")
        print(f"- p-max-depth: {max_depth}")
        print(f"- Metadata column: {md_col}")
        print(f"- References by level: {refs_str}")
        print("STATUS:")
        for s in self.st.status_msgs:
            print(f"- {s}")
        if self.st.errors:
            print("---")
            print("ERRORS:")
            for e in self.st.errors:
                print(f"- {e}")
        print("---")
        print(f"Full log: {self.st.logfile}")
        print("\n--- Pipeline complete! ---")

    # ---- Helper ----
    def validate_and_continue(self, ok: bool):
        if ok:
            print("Step validated. Proceeding to next step...")
        else:
            print("Step has errors. Proceeding with best-effort fallback.")

# -----------------------------
# Main
# -----------------------------

def main():
    print("== QIIME2 Interactive Metagenomics Pipeline ==")
    print("Fixed version based on meta.py, meta2_no_phylogeny.py, meta3_no_phylogeny.py")
    print("Note: type 'exit' at any prompt to safely interrupt the pipeline.\n")

    st = PipelineState()
    pipe = QIIME2Pipeline(st)

    pipe.step_import()
    pipe.step_cutadapt()
    pipe.step_dada2()
    pipe.step_taxonomy()
    pipe.step_taxa_filter()
    pipe.step_taxa_collapse()
    pipe.step_diversity()
    pipe.step_differential()
    pipe.step_summary()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nKeyboard interrupt. Safe exit.")
        sys.exit(1)
