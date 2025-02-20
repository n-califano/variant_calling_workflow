SAMPLES = ["SRR2584857_1"]
RAW_DATA_SAMPLE_DIR = "raw_data/sample_data"
REFERENCE_FILE = "raw_data/ref_data/ecoli-rel606.fa.gz"
FASTQC_DIR = "results/fastqc"

rule all:
    input:
        expand(f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz", sample=SAMPLES),
        REFERENCE_FILE,
        expand(f"{FASTQC_DIR}/{{sample}}_fastqc.html", sample=SAMPLES)

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



