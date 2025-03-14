import os.path
import scripts.python.utils as utils

configfile: "config/config.yml"

URL_SAMPLE = config['download_url_sample']
URL_REF = config['download_url_ref']
URL_LCR_FILE = config['download_lcr_url']

TUMOR_NORMAL_PAIRS = config['tumor_normal_pairs']
ALL_SAMPLES = [sample for pair in TUMOR_NORMAL_PAIRS.items() for sample in pair]

RAW_DATA_REF_DIR = config['raw_data_ref_dir']
RAW_DATA_SAMPLE_DIR = config['raw_data_sample_dir']
FASTQC_DIR = config['fastqc_dir']
MINIMAP_DIR = config['minimap_dir']
PICARD_DIR = config['picard_dir']
MULTIQC_DIR = config['multiqc_dir']
GATK_DIR = config['gatk_dir']

REFERENCE_FILE = os.path.join(RAW_DATA_REF_DIR, os.path.basename(URL_REF))
REF_FILE_NAME = utils.remove_all_extensions(os.path.basename(REFERENCE_FILE))
LCR_FILE = os.path.join(RAW_DATA_REF_DIR, os.path.basename(URL_LCR_FILE))
LCR_FILE_NAME = utils.remove_all_extensions(os.path.basename(LCR_FILE))
UNZIPPED_REF_FILE = f"{RAW_DATA_REF_DIR}/{REF_FILE_NAME}.fa"
REF_DICTIONARY = f"{RAW_DATA_REF_DIR}/{REF_FILE_NAME}.dict"


rule all:
    input:
        #expand(f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz", sample=config['samples']),
        REFERENCE_FILE,
        #expand(f"{FASTQC_DIR}/{{sample}}_fastqc.html", sample=config['samples']),
        #expand(f"{MINIMAP_DIR}/{{sample}}.{REF_FILE_NAME}_align.sam", sample=config['samples']),
        #expand(f"{PICARD_DIR}/{{sample}}_coord_sorted.bam", sample=config['samples'])
        #expand(f"{PICARD_DIR}/{{sample}}_align_summary_metrics.txt", sample=config['samples'])
        expand(f"{MULTIQC_DIR}/multiqc_report.html", sample=ALL_SAMPLES),
        f"{RAW_DATA_REF_DIR}/{REF_FILE_NAME}.dict",
        #expand(f"{GATK_DIR}/{{sample}}.vcf", sample=ALL_SAMPLES)
        expand(f"{GATK_DIR}/{{tumor}}_vs_{{normal}}_passing_lcr_filtered.vcf", tumor=TUMOR_NORMAL_PAIRS.keys(), normal=TUMOR_NORMAL_PAIRS.values()),
        f"{RAW_DATA_REF_DIR}/{LCR_FILE_NAME}.bed"
    #output:
     #   expand(f"{PICARD_DIR}/{{sample}}_align_summary_metrics.txt", sample=config['samples'])

rule download_sample_data:
    output: f"{RAW_DATA_SAMPLE_DIR}/{{sample}}.fastq.gz"
    params: 
        sample_url=f"{URL_SAMPLE}/{{sample}}.fastq.gz"
    shell: "wget -O {output} {params.sample_url}"

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
    shell: "minimap2 -ax sr -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' {REFERENCE_FILE} {input} > {output}"
    #old version: -R '@RG\\tID:{params.sample}\\tSM:{params.sample}\\tPL:illumina'  the PL part seems to not be necessary

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

rule aggregate_qc:
    input:
        align_summary_metrics=expand(f"{PICARD_DIR}/{{sample}}_align_summary_metrics.txt", sample=ALL_SAMPLES),
        fastqc=expand(f"{FASTQC_DIR}/{{sample}}_fastqc.zip", sample=ALL_SAMPLES)
    output: f"{MULTIQC_DIR}/multiqc_report.html"
    shell: "multiqc {input} --outdir {MULTIQC_DIR}"

rule create_reference_index:
    input: REFERENCE_FILE
    output: 
        UNZIPPED_REF_FILE,
        ref_index_file=f"{RAW_DATA_REF_DIR}/{REF_FILE_NAME}.fa.fai"
    shell: 
        """
        gunzip -c {input} > {UNZIPPED_REF_FILE}
        samtools faidx {UNZIPPED_REF_FILE}
        """

rule create_reference_dict:
    output: REF_DICTIONARY
    shell: "picard CreateSequenceDictionary \
                --REFERENCE {UNZIPPED_REF_FILE} \
                --OUTPUT {REF_DICTIONARY}"

rule variant_calling:
    input:
        coord_sorted_normal_bam=f"{PICARD_DIR}/{{normal}}_coord_sorted.bam",
        coord_sorted_tumor_bam=f"{PICARD_DIR}/{{tumor}}_coord_sorted.bam",
    output:
        vcf_file=f"{GATK_DIR}/{{tumor}}_vs_{{normal}}_raw.vcf"
    shell:
        """
        gatk Mutect2 \
            --sequence-dictionary {REF_DICTIONARY} \
            --reference {UNZIPPED_REF_FILE} \
            --input {input.coord_sorted_normal_bam} \
            --normal-sample {wildcards.normal} \
            --input {input.coord_sorted_tumor_bam} \
            --tumor-sample {wildcards.tumor} \
            --annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample --annotation Coverage \
            --output {output.vcf_file}
        """

rule download_LCR_file:
    output: 
        lcr_file_with_chr=f"{RAW_DATA_REF_DIR}/{LCR_FILE_NAME}_with_chr.bed",
        lcr_file=f"{RAW_DATA_REF_DIR}/{LCR_FILE_NAME}.bed"
    shell:
        """
        wget -O {output.lcr_file_with_chr}.gz {URL_LCR_FILE} 
        gunzip -c {output.lcr_file_with_chr}.gz > {output.lcr_file_with_chr}
        sed 's/^chr//g' {output.lcr_file_with_chr} > {output.lcr_file}
        """

rule variant_filtering:
    input: 
        UNZIPPED_REF_FILE,
        vcf_file=f"{GATK_DIR}/{{tumor}}_vs_{{normal}}_raw.vcf",
        lcr_file=f"{RAW_DATA_REF_DIR}/{LCR_FILE_NAME}.bed"
    output:
        mutect_filtered_vcf_file=f"{GATK_DIR}/{{tumor}}_vs_{{normal}}_mutect_filtered.vcf",
        pass_filtered_vcf_file=f"{GATK_DIR}/{{tumor}}_vs_{{normal}}_passing_filtered.vcf",
        lcr_filtered_vcf_file=f"{GATK_DIR}/{{tumor}}_vs_{{normal}}_passing_lcr_filtered.vcf"
    shell:
        """
        gatk FilterMutectCalls \
            --reference {UNZIPPED_REF_FILE} \
            --variant {input.vcf_file} \
            --output {output.mutect_filtered_vcf_file}

        # Keep only the SNP marked as PASS
        SnpSift filter \
            -noLog \
            "( FILTER = 'PASS' )" \
            {output.mutect_filtered_vcf_file} > {output.pass_filtered_vcf_file}

        # Filter out LCR (Low Complexity Regions)
        SnpSift intervals \
            -noLog \
            -x \
            -i {output.pass_filtered_vcf_file} \
            {input.lcr_file} > {output.lcr_filtered_vcf_file}
        """

