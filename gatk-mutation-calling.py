# Setting up filesnames here:
from os.path import join
configfile: "Melanoma.config.json"

# Tools
SAMTOOLS = "samtools"

# Files
SAMPLES = ["WT-1", "YM-1"]
normal_skin = ["WT-1"]
tumor = ["YM-1"]

#Reference File
FASTA = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/reference/my_dir/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna",

# path to where bams are
bams_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/sorted_bam/"
intermediate_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/gatk_intermediates/"
fasta_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/reference/my_dir/GCF_000001635.26_GRCm38.p6/"
peptide_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/peptides/"

# Path to packages not in conda
VarScan_path = "/scratch/external_scripts/VarScan.v2.3.9.jar"
perl_filter_path = "/scratch/external_scripts/fpfilter-tool-master/fpfilter-2.pl"
perl5lib_path = "/packages/6x/vcftools/0.1.12b/lib/per15/site_perl/"
gatk_path = "~/gatk-4.1.7.0/"

rule all:
    input:
        expand(fasta_path + "GCF_000001635.26_GRCm38.p6_genomic.dict", fasta_path=fasta_path),
        expand(intermediate_path + "{tumor}-somatic.vcf.gz", intermediate_path=intermediate_path,tumor=tumor),
        expand(intermediate_path + "{tumor}_filtered.vcf", intermediate_path=intermediate_path,tumor=tumor),
        expand(intermediate_path + "{tumor}_filtered_PASS.noheader.vcf", intermediate_path=intermediate_path, tumor=tumor),
        expand(intermediate_path + "{tumor}_filtered_noheader_vep.vcf", intermediate_path=intermediate_path, tumor=tumor),
        expand(peptide_path + "{tumor}_Gatk_vep.15.peptide", peptide_path=peptide_path, tumor=tumor),
        expand(peptide_path + "{tumor}_Gatk_vep.17.peptide", peptide_path=peptide_path, tumor=tumor),
        expand(peptide_path + "{tumor}_Gatk_vep.19.peptide", peptide_path=peptide_path, tumor=tumor),
        expand(peptide_path + "{tumor}_Gatk_vep.21.peptide", peptide_path=peptide_path, tumor=tumor)


rule fasta_dict:
    input:
        fa = FASTA
    output:
        dict = os.path.join(fasta_path, "GCF_000001635.26_GRCm38.p6_genomic.dict")
    shell:
        """
        java -jar $PICARD CreateSequenceDictionary REFERENCE={input.fa} OUTPUT={output.dict}
        """

rule gatk_mutations:
    input:
        fa = FASTA,
        tumor_bam = os.path.join(bams_path, "{tumor}_WES_HISAT2_aligned_sortedbycoord_RG.bam"),
        normal_bam = os.path.join(bams_path, "WT-1_WES_HISAT2_aligned_sortedbycoord_RG.bam"),
    output:
        somatic = os.path.join(intermediate_path, "{tumor}-somatic.vcf.gz")
    params:
        gatk = os.path.join(gatk_path, "gatk Mutect2")
    shell:
        """
        {params.gatk} -R {input.fa} -I {input.tumor_bam} -I {input.normal_bam} -normal WT -O {output.somatic}
        """

rule filter:
    input:
        fa = FASTA,
        unfiltered = os.path.join(intermediate_path, "{tumor}-somatic.vcf.gz")
    output:
        filtered= os.path.join(intermediate_path, "{tumor}_filtered.vcf")
    shell:
        """
        gatk FilterMutectCalls -R {input.fa} -V {input.unfiltered} -O {output.filtered}
        """

rule grep_pass:
    input:
        os.path.join(intermediate_path, "{tumor}_filtered.vcf")
    output:
        os.path.join(intermediate_path, "{tumor}_filtered_PASS.noheader.vcf")
    shell:
        """
        grep "PASS" {input} > {output}
        """

rule run_vep:
    input:
        os.path.join(intermediate_path, "{tumor}_filtered_PASS.noheader.vcf")
    output:
        os.path.join(intermediate_path, "{tumor}_filtered_noheader_vep.vcf")
    shell:
        """
        vep -i {input} --format vcf --species mus_musculus --cache --offline --vcf -o {output} --force_overwrite  --symbol --plugin Wildtype --terms SO --plugin Downstream
        """ 

rule generate_fasta:
    input:
        os.path.join(intermediate_path, "{tumor}_filtered_noheader_vep.vcf")
    output:
        len_15 = os.path.join(peptide_path, "{tumor}_Gatk_vep.15.peptide"),
        len_17 = os.path.join(peptide_path, "{tumor}_Gatk_vep.17.peptide"),
        len_19 = os.path.join(peptide_path, "{tumor}_Gatk_vep.19.peptide"),
        len_21 = os.path.join(peptide_path, "{tumor}_Gatk_vep.21.peptide"),
    shell:
       """
        pvacseq generate_protein_fasta {input} 15 {output.len_15};
        pvacseq generate_protein_fasta {input} 17 {output.len_17};
        pvacseq generate_protein_fasta {input} 19 {output.len_19};
        pvacseq generate_protein_fasta {input} 21 {output.len_21};
        """
