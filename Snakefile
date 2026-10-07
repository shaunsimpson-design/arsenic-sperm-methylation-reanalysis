# RRBS pipeline: gestational arsenic, mouse sperm (GSE150650 / SRP261822)
# Run from the project root:  snakemake --cores 14
import csv, os

with open("samples.csv") as f:
    RUNS = [row["run"] for row in csv.DictReader(f)]

GENOME = os.path.abspath("genome/mm10")
ALIGNED = os.path.abspath("data/aligned")
METH = os.path.abspath("data/methylation")

rule all:
    input:
        expand("data/methylation/{run}_1_val_1_bismark_bt2_pe.bismark.cov.gz", run=RUNS)

# 1. Download from SRA and convert to FASTQ (deleted after trimming)
rule download:
    output:
        r1=temp("data/raw/{run}_1.fastq"),
        r2=temp("data/raw/{run}_2.fastq")
    threads: 14
    priority: 1
    log: "logs/download/{run}.log"
    shell:
        "mkdir -p data/raw && cd data/raw && prefetch {wildcards.run} > ../../{log} 2>&1 "
        "&& fasterq-dump {wildcards.run} >> ../../{log} 2>&1 "
        "&& rm -rf {wildcards.run}"

# 2. Trim: RRBS mode, poly-G off (HiSeq X Ten is 4-colour)
rule trim:
    input:
        r1="data/raw/{run}_1.fastq",
        r2="data/raw/{run}_2.fastq"
    output:
        r1=temp("data/trimmed/{run}_1_val_1.fq"),
        r2=temp("data/trimmed/{run}_2_val_2.fq"),
        rep1="data/trimmed/{run}_1.fastq_trimming_report.txt",
        rep2="data/trimmed/{run}_2.fastq_trimming_report.txt"
    threads: 14
    priority: 2
    log: "logs/trim/{run}.log"
    shell:
        "trim_galore --rrbs --paired --no_poly_g --fastqc -o data/trimmed "
        "{input.r1} {input.r2} > {log} 2>&1"

# 3. Align: one Bismark instance, 6 Bowtie2 threads each (fits 24 GB RAM)
rule align:
    input:
        r1="data/trimmed/{run}_1_val_1.fq",
        r2="data/trimmed/{run}_2_val_2.fq"
    output:
        bam="data/aligned/{run}_1_val_1_bismark_bt2_pe.bam",
        report="data/aligned/{run}_1_val_1_bismark_bt2_PE_report.txt"
    threads: 14
    priority: 3
    log: "logs/align/{run}.log"
    shell:
        "cd data/trimmed && bismark -p 6 --genome {GENOME} -o {ALIGNED} "
        "-1 {wildcards.run}_1_val_1.fq -2 {wildcards.run}_2_val_2.fq "
        "> ../../{log} 2>&1"

# 4. Extract per-CpG calls, excluding M-bias positions found in C1;
#    then delete the large per-context files (only the .cov.gz is needed)
rule extract:
    input:
        "data/aligned/{run}_1_val_1_bismark_bt2_pe.bam"
    output:
        cov="data/methylation/{run}_1_val_1_bismark_bt2_pe.bismark.cov.gz",
        mbias="data/methylation/{run}_1_val_1_bismark_bt2_pe.M-bias.txt"
    threads: 14
    priority: 4
    log: "logs/extract/{run}.log"
    shell:
        "bismark_methylation_extractor -p --gzip --bedGraph "
        "--ignore_r2 2 --ignore_3prime 2 --ignore_3prime_r2 1 "
        "-o {METH} {input} > {log} 2>&1 "
        "&& rm -f {METH}/C??_*_{wildcards.run}_1_val_1_bismark_bt2_pe.txt.gz"
