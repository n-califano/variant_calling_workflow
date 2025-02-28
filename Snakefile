import os.path
import scripts.python.utils as utils

configfile: "config/config.yml"

URL_SAMPLE = config['download_url_sample']
URL_REF = config['download_url_ref']
REFERENCE_FILE = config['reference_file']
REF_FILE_NAME = utils.remove_all_extensions(os.path.basename(REFERENCE_FILE))
RAW_DATA_SAMPLE_DIR = config['raw_data_sample_dir']
FASTQC_DIR = config['fastqc_dir']
MINIMAP_DIR = config['minimap_dir']
PICARD_DIR = config['picard_dir']

rule all:
    input:
        expand(f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz", sample=config['samples']),
        REFERENCE_FILE,
        expand(f"{FASTQC_DIR}/{{sample}}_fastqc.html", sample=config['samples']),
        expand(f"{MINIMAP_DIR}/{{sample}}.{REF_FILE_NAME}_align.sam", sample=config['samples']),
        expand(f"{PICARD_DIR}/{{sample}}_coord_sorted.bam", sample=config['samples'])

rule download_sample_data:
    output: f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz"
    shell: "wget -O {output} {URL_SAMPLE}"

rule download_ref_data:
    output: REFERENCE_FILE
    shell: "wget -O {output} {URL_REF}"

rule run_fastqc:
    input: f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz"
    output: f"{FASTQC_DIR}/{{sample}}_fastqc.html",
            f"{FASTQC_DIR}/{{sample}}_fastqc.zip",
    shell: "fastqc -o {FASTQC_DIR} -t 4 {input}"

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
        coord_sorted_index=f"{PICARD_DIR}/{{sample}}_coord_sorted.bai",
        align_summary_metrics=f"{PICARD_DIR}/{{sample}}_align_summary_metrics.txt"
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

        picard CollectAlignmentSummaryMetrics \
            --INPUT {output.coord_sorted} \
            --REFERENCE_SEQUENCE {REFERENCE_FILE} \
            --OUTPUT {output.align_summary_metrics}
        """