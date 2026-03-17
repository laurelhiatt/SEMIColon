import gzip
import pysam
import argparse
import os

PRIMARY_CHR = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def compute_gc(file_path, fasta, depth_threshold=8):

    gc = 0
    total = 0
    lines_seen = 0

    with gzip.open(file_path, "rt") as f:
        for line in f:
            chrom, start, end, depth = line.strip().split()
            lines_seen += 1

            if chrom not in PRIMARY_CHR:
                continue

            depth = int(depth)
            if depth < depth_threshold:
                continue

            start = int(start)
            end = int(end)

            seq = fasta.fetch(chrom, start, end).upper()

            gc += seq.count("G") + seq.count("C")
            total += len(seq)
            if lines_seen % 5_000_000 == 0:
                print(f"  processed {lines_seen:,} lines...", flush=True)

    gc_content = gc / total if total > 0 else 0
    return gc_content, total


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--depth_file", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    fasta = pysam.FastaFile(args.reference)
    print(f"Processing {args.output} from {args.depth_file}.", flush=True)
    gc_content, total = compute_gc(args.depth_file, fasta)

    basename = os.path.basename(args.depth_file).replace(".per-base.bed.gz", "")
    print(f"GC content for {basename}: {gc_content:.6f} (total bases: {total})" )
    with open(args.output, "w") as out:
        out.write("file_name\tGC_content\ttotal_bases\n")
        out.write(f"{basename}\t{gc_content:.6f}\t{total}\n")

if __name__ == "__main__":
    main()