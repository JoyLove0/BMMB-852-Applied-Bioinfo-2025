# Variable definitions
REF=refs/${CHR}.fasta
BAM=bam/KRAS-N-P.bam
VCF=vcf/KRAS-N-P.vcf.gz

# Run the variant calling for the region of interest.
make -f src/run/bcftools.mk \
        REF=${REF} \
        BAM=${BAM} \
        VCF=${VCF} \
        run