#!/usr/bin/env bash
set -euo pipefail

BAM_DIR="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/bam"
OUT_DIR="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/pks"                             # where per-sample pks BAMs will be written
PKS_INDEX="pks_index"
GTF="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/clb_genes.gtf"
# parallelization:
JOBS=10
THREADS_PER_JOB=8
FEATURECOUNTS_THREADS=2 # threads for the final featureCounts run

# Strandedness setting for featureCounts:
# 0 = unstranded (default), 1 = stranded, 2 = reversely stranded
STRANDED=0

# Behavior
CLEANUP=true            # remove intermediate FASTQs / unsorted BAMs
TMPDIR="${OUT_DIR}/pks_temp"
FORCE=false             # if false and final aggregated file exists, script WILL SKIP
OUT_FC_BASENAME="clb_featureCounts_aggregated.txt"
OUT_FC="${OUT_DIR}/${OUT_FC_BASENAME}"

# === TESTING: Uncomment and set a few specific BAM paths to test the pipeline on a handful of samples.
# Leave commented to process all '*-sorted.bam' in $BAM_DIR.
TEST_BAMS_FILE=""
TEST_BAMS_FILE=("/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/bam/12b_D12_84_TR-sorted.bam")

# ========================================
# Tool lookups (must be on PATH)
# ========================================
SAMTOOLS="$(command -v samtools || true)"
BOWTIE2="$(command -v bowtie2 || true)"
FEATURECOUNTS="$(command -v featureCounts || true)"
PARALLEL="$(command -v parallel || true)"
XARGS="$(command -v xargs || true)"

if [[ -z "$SAMTOOLS" || -z "$BOWTIE2" || -z "$FEATURECOUNTS" ]]; then
  echo "ERROR: samtools, bowtie2 and featureCounts must be in PATH." >&2
  echo "Found: samtools='$SAMTOOLS' bowtie2='$BOWTIE2' featureCounts='$FEATURECOUNTS'" >&2
  exit 1
fi

# If aggregated output exists and we're not forcing, skip entire run
if [[ -f "$OUT_FC" && "$FORCE" != "true" ]]; then
  echo "Aggregated featureCounts file already exists at: $OUT_FC"
  echo "Set FORCE=true to regenerate. Exiting."
  exit 0
fi

# Choose parallel backend
USE_PARALLEL=false
if [[ -n "$PARALLEL" ]]; then
  USE_PARALLEL=true
  echo "Using GNU parallel ($PARALLEL) with $JOBS jobs."
elif [[ -n "$XARGS" ]]; then
  echo "GNU parallel not found; will use xargs -P ($XARGS) with $JOBS jobs."
else
  echo "Neither 'parallel' nor 'xargs' found; will run sequentially (JOBS=1)."
  JOBS=1
fi

mkdir -p "$TMPDIR"
mkdir -p "$OUT_DIR"

# helper function that processes one BAM (extract paired unmapped pairs, align, sort, index)
process_sample() {
  bam="$1"
  SAMTOOLS="$2"
  BOWTIE2="$3"
  PKS_INDEX="$4"
  OUT_DIR="$5"
  TMPDIR="$6"
  THREADS="$7"
  CLEANUP="$8"

  set -euo pipefail

  fname=$(basename "$bam")
  sample="${fname%-sorted.bam}"

  echo ">>> [$sample] START: $(date)"

  unsorted_bam="${TMPDIR}/${sample}-unsorted-unmapped.bam"
  fastq_r1="${TMPDIR}/${sample}_unmapped_R1.fastq"
  fastq_r2="${TMPDIR}/${sample}_unmapped_R2.fastq"
  pks_bam="${OUT_DIR}/${sample}-pks.bam"
  pks_sorted="${OUT_DIR}/${sample}-pks-sorted.bam"
  marker="${OUT_DIR}/${sample}.pks_done"

  # Skip if sample is already processed (marker exists)
  if [[ -f "$marker" ]]; then
    echo "[$sample] marker exists ($marker) — skipping sample."
    return 0
  fi

  # 1) extract reads where both mates are unmapped (flag 0xC = 12)
  echo "[$sample] extracting paired-unmapped reads -> $unsorted_bam"
  "$SAMTOOLS" view -b -f 12 "$bam" -o "$unsorted_bam" || true

  # If no reads were extracted, create an empty sorted BAM so featureCounts has a file to process
  if [[ ! -s "$unsorted_bam" ]]; then
    echo "[$sample] no paired unmapped reads found. Creating an empty sorted BAM as placeholder."
    # Create empty BAM by converting /dev/null (samtools view will error on /dev/null so use header trick)
    # Create a minimal empty bam by copying header from original bam and producing 0 alignments
    "$SAMTOOLS" view -H "$bam" > "${TMPDIR}/${sample}.header.sam"
    # write an empty BAM with header only (samtools view requires input, so pipe nothing)
    cat "${TMPDIR}/${sample}.header.sam" | "$SAMTOOLS" view -b -o "$pks_bam" -
    rm -f "${TMPDIR}/${sample}.header.sam"

    # sort & index
    "$SAMTOOLS" sort -@ "$THREADS" -o "$pks_sorted" "$pks_bam"
    "$SAMTOOLS" index "$pks_sorted"

    touch "$marker"
    rm -f "$unsorted_bam"
    echo "[$sample] placeholder pks-sorted BAM created: $pks_sorted"
    return 0
  fi

  # 2) convert paired BAM -> paired FASTQ (samtools fastq -1 -2)
  echo "[$sample] BAM -> paired FASTQ: $fastq_r1 , $fastq_r2"
  rm -f "$fastq_r1" "$fastq_r2"
  "$SAMTOOLS" fastq -1 "$fastq_r1" -2 "$fastq_r2" -0 /dev/null -s /dev/null -n "$unsorted_bam"

  if [[ ! -s "$fastq_r1" && ! -s "$fastq_r2" ]]; then
    echo "[$sample] produced empty FASTQs; creating placeholder BAM and skipping alignment."
    # create placeholder empty sorted BAM
    "$SAMTOOLS" view -H "$bam" > "${TMPDIR}/${sample}.header.sam"
    cat "${TMPDIR}/${sample}.header.sam" | "$SAMTOOLS" view -b -o "$pks_bam" -
    rm -f "${TMPDIR}/${sample}.header.sam"
    "$SAMTOOLS" sort -@ "$THREADS" -o "$pks_sorted" "$pks_bam"
    "$SAMTOOLS" index "$pks_sorted"
    touch "$marker"
    rm -f "$unsorted_bam"
    return 0
  fi

  # 3) align paired reads with bowtie2 -> BAM
  echo "[$sample] bowtie2 aligning -> $pks_bam"
  "$BOWTIE2" -x "$PKS_INDEX" -1 "$fastq_r1" -2 "$fastq_r2" -p "$THREADS" \
    | "$SAMTOOLS" view -b -o "$pks_bam" -

  # 4) sort & index
  echo "[$sample] sorting -> $pks_sorted"
  "$SAMTOOLS" sort -@ "$THREADS" -o "$pks_sorted" "$pks_bam"
  echo "[$sample] indexing -> ${pks_sorted}.bai"
  "$SAMTOOLS" index "$pks_sorted"

  # mark done
  touch "$marker"

  # cleanup intermediates if requested
  if [[ "$CLEANUP" == "true" ]]; then
    rm -f "$unsorted_bam" "$fastq_r1" "$fastq_r2" "$pks_bam"
  fi

  echo ">>> [$sample] DONE: $(date)"
}

export -f process_sample

# Gather list of BAMs to process
if [[ -n "$TEST_BAMS_FILE" && -f "$TEST_BAMS_FILE" ]]; then
  echo "Using TEST_BAMS_FILE: $TEST_BAMS_FILE"
  mapfile -t BAMS < "$TEST_BAMS_FILE"
else
  mapfile -t BAMS < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*-sorted.bam" | sort)
fi


if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "No '*-sorted.bam' files found in $BAM_DIR (and TEST_BAMS not set). Exiting."
  exit 1
fi

echo "Will process ${#BAMS[@]} BAM(s)."

# Run in parallel
if [[ "$USE_PARALLEL" == true ]]; then
  printf "%s\n" "${BAMS[@]}" \
    | parallel -j "$JOBS" --halt now,fail=1 \
       process_sample {} "$SAMTOOLS" "$BOWTIE2" "$PKS_INDEX" "$OUT_DIR" "$TMPDIR" "$THREADS_PER_JOB" "$CLEANUP"
else
  printf "%s\n" "${BAMS[@]}" \
    | xargs -n1 -P "$JOBS" -I{} bash -c 'process_sample "$@"' _ {} "$SAMTOOLS" "$BOWTIE2" "$PKS_INDEX" "$OUT_DIR" "$TMPDIR" "$THREADS_PER_JOB" "$CLEANUP"
fi

# wait a moment for jobs to finish
sleep 2

# collect pks-sorted bams
mapfile -t PKS_BAMS < <(find "$OUT_DIR" -maxdepth 1 -type f -name "*-pks-sorted.bam" | sort)

if [[ ${#PKS_BAMS[@]} -eq 0 ]]; then
  echo "No '*-pks-sorted.bam' files found after processing; nothing to run featureCounts on."
  exit 0
fi

echo "Running aggregated featureCounts on ${#PKS_BAMS[@]} BAMs with $FEATURECOUNTS_THREADS threads."
echo "featureCounts will count fragments (-p) with strandedness STRANDED=$STRANDED"

# Run featureCounts with fragment counting and user-set strandedness
"$FEATURECOUNTS" -T "$FEATURECOUNTS_THREADS" -p -s "$STRANDED" -t CDS -g gene_id -a "$GTF" -o "$OUT_FC" "${PKS_BAMS[@]}"

echo "featureCounts finished -> $OUT_FC"

# done: optional cleanup of tmp dir if empty
if [[ "$CLEANUP" == "true" ]]; then
  rmdir --ignore-fail-on-non-empty "$TMPDIR" || true
fi

echo "ALL DONE."
