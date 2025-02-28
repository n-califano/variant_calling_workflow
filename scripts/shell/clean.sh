SCRIPT_DIR="$(dirname $0)"
rm "$SCRIPT_DIR/../../results/fastqc/"*.{html,zip} \
   "$SCRIPT_DIR/../../results/minimap/"*.sam \
   "$SCRIPT_DIR/../../results/picard/"*.{bam,bai,txt}