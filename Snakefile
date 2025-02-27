SAMPLES = ["SRR2584857_1", "SRR2584857_2"]
RAW_DATA_SAMPLE_DIR = "raw_data/sample_data"
REF_FILE_NAME = "ecoli-rel606"
REFERENCE_FILE = "raw_data/ref_data/" + REF_FILE_NAME + ".fa.gz"
FASTQC_DIR = "results/fastqc"
MINIMAP_DIR = "results/minimap"

rule all:
    input:
        expand(f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz", sample=SAMPLES),
        REFERENCE_FILE,
        expand(f"{FASTQC_DIR}/{{sample}}_fastqc.html", sample=SAMPLES),
        expand(f"{MINIMAP_DIR}/{{sample}}.{REF_FILE_NAME}_align.sam", sample=SAMPLES)

rule download_sample_data:
    output: expand(f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz", sample=SAMPLES)
    shell: "wget -O {output} https://osf.io/4rdza/download"

rule download_ref_data:
    output: REFERENCE_FILE
    shell: "wget -O {output} https://osf.io/8sm92/download"

rule run_fastqc:
    input: expand(f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz", sample=SAMPLES)
    output: expand(f"{FASTQC_DIR}/{{sample}}_fastqc.html", sample=SAMPLES), 
            expand(f"{FASTQC_DIR}/{{sample}}_fastqc.zip", sample=SAMPLES)
    shell: f"fastqc -o {FASTQC_DIR} -t 4 {input}"

rule align_reads:
    input: f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz"
    output: f"{MINIMAP_DIR}/{{sample}}.{REF_FILE_NAME}_align.sam"
    shell: "minimap2 -ax sr {REFERENCE_FILE} {input} > {output}"

