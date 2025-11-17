# The statistics for the gold standard.
VCF=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz
bcftools view $VCF chr12:24000000-26000000 | bcftools stats

# The statistics for the currentvariant calls.
bcftools view vcf/both.vcf.gz | bcftools stats