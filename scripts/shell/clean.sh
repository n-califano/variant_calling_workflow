SCRIPT_DIR="$(dirname $0)"
rm -rf "$SCRIPT_DIR/../../results/fastqc/"*.{html,zip} \
   "$SCRIPT_DIR/../../results/minimap/"*.sam \
   "$SCRIPT_DIR/../../results/picard/"*.{bam,bai,txt} \
   "$SCRIPT_DIR/../../results/multiqc/"{*.html,multiqc_data*/} \
   "$SCRIPT_DIR/../../raw_data/ref_data/"*.{dict,fa,fai,bed} \
   "$SCRIPT_DIR/../../results/gatk/"*.vcf* \