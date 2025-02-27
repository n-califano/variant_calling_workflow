#TODO: move parameters in configfile

SAMPLES = ["SRR2584857_1", "SRR2584857_2"]
RAW_DATA_SAMPLE_DIR = "raw_data/sample_data"
REF_FILE_NAME = "ecoli-rel606"
REFERENCE_FILE = "raw_data/ref_data/" + REF_FILE_NAME + ".fa.gz"
FASTQC_DIR = "results/fastqc"
MINIMAP_DIR = "results/minimap"
PICARD_DIR = "results/picard"

rule all:
    input:
        expand(f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz", sample=SAMPLES),
        REFERENCE_FILE,
        expand(f"{FASTQC_DIR}/{{sample}}_fastqc.html", sample=SAMPLES),
        expand(f"{MINIMAP_DIR}/{{sample}}.{REF_FILE_NAME}_align.sam", sample=SAMPLES),
        expand(f"{PICARD_DIR}/{{sample}}_coord_sorted.bam", sample=SAMPLES)

rule download_sample_data:
    output: f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz"
    shell: "wget -O {output} https://osf.io/4rdza/download"

rule download_ref_data:
    output: REFERENCE_FILE
    shell: "wget -O {output} https://osf.io/8sm92/download"

rule run_fastqc:
    input: f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz"
    output: f"{FASTQC_DIR}/{{sample}}_fastqc.html",
            f"{FASTQC_DIR}/{{sample}}_fastqc.zip",
    shell: f"fastqc -o {FASTQC_DIR} -t 4 {input}"

rule align_reads:
    input: f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz"
    output: f"{MINIMAP_DIR}/{{sample}}.{REF_FILE_NAME}_align.sam"
    shell: "minimap2 -ax sr {REFERENCE_FILE} {input} > {output}"

rule process_sam:
    input:
        sam=f"{MINIMAP_DIR}/{{sample}}.{REF_FILE_NAME}_align.sam"
    output:
        query_sorted=f"{PICARD_DIR}/{{sample}}_query_sorted.bam",
        dedup=f"{PICARD_DIR}/{{sample}}_dedup.bam",
        metrics=f"{PICARD_DIR}/{{sample}}_metrics.txt",
        coord_sorted=f"{PICARD_DIR}/{{sample}}_coord_sorted.bam",
        coord_sorted_index=f"{PICARD_DIR}/{{sample}}_coord_sorted.bai"
    shell:
        """
        picard SortSam \
            --INPUT {input.sam} \
            --OUTPUT {output.query_sorted} \
            --SORT_ORDER queryname

        picard MarkDuplicates \
            --INPUT {output.query_sorted} \
            --OUTPUT {output.dedup} \
            --METRICS_FILE {output.metrics} \
            --REMOVE_DUPLICATES true

        picard SortSam \
            --INPUT {output.dedup} \
            --OUTPUT {output.coord_sorted} \
            --SORT_ORDER coordinate \
            --CREATE_INDEX true
        """