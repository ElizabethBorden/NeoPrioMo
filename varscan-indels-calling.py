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
bams_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/bam/"
pileup_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/pileup/"
bai_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/bam/"
intermediate_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/intermediate/"
bamreadcount_path = "bam-readcount"
peptide_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/peptides/"

# Path to packages not in conda
VarScan_path = "~/external_scripts/VarScan.v2.3.9.jar"
perl_filter_path = "~/external_scripts/fpfilter-tool-master/fpfilter-2.pl"
perl5lib_path = "/packages/6x/vcftools/0.1.12b/lib/per15/site_perl/"

rule all:
    input:
        expand(intermediate_path +"{filename}.Varscan.indel.Somatic.hc", intermediate_path=intermediate_path,  filename=tumor),
        expand(intermediate_path +"{filename}.Varscan.indel.Somatic.vcf", intermediate_path=intermediate_path,  filename=tumor),
        expand(intermediate_path + "{filename}.Varscan.indel.Somatic.vcf_TP_after_filtered_sorted.vcf", intermediate_path=intermediate_path,  filename=tumor),
        expand(intermediate_path +"{filename}_Varscan_indel_variants_noheader.vcf",intermediate_path=intermediate_path, filename=tumor),
        expand(intermediate_path +"{filename}_Varscan_indel_vep.vcf",intermediate_path=intermediate_path, filename=tumor),
#        expand(intermediate_path +"{filename}_Varscan_indel_vep.tsv", intermediate_path=intermediate_path, filename=tumor), 
        expand(peptide_path +"{filename}_VarScan_indel_vep.15.peptide", intermediate_path=intermediate_path, filename=tumor),
        expand(peptide_path +"{filename}_VarScan_indel_vep.17.peptide", intermediate_path=intermediate_path, filename=tumor),
        expand(peptide_path +"{filename}_VarScan_indel_vep.19.peptide", intermediate_path=intermediate_path, filename=tumor),
        expand(peptide_path +"{filename}_VarScan_indel_vep.21.peptide", intermediate_path=intermediate_path, filename=tumor)

rule isolate_calls_by_type_and_confidence:
    input:
        VarScan_indel = os.path.join(intermediate_path, "{filename}.Varscan.indel") 
    output:
        VarScan_indel = os.path.join(intermediate_path, "{filename}.Varscan.indel.Somatic.hc")
    params:
        VarScan = VarScan_path
    shell:
        """
        java -jar {params.VarScan} processSomatic {input.VarScan_indel}
        """

rule convert_to_vcf:
    input:
        os.path.join(intermediate_path, "{filename}.Varscan.indel.Somatic.hc")
    output:
        os.path.join(intermediate_path, "{filename}.Varscan.indel.Somatic.vcf")
    shell:                                                                                                      """
        python ~/external_scripts/vscan2vcf.py {input} > {output}
        """

rule FP_filter:
    input:
        fa = FASTA,
        unfiltered = os.path.join(intermediate_path, "{filename}.Varscan.indel.Somatic.vcf")
    output:
        filtered= os.path.join(intermediate_path, "{filename}.Varscan.indel.Somatic.vcf_TP_after_filtered_sorted.vcf")
    params:
        FP = "~/miniconda3/envs/Filtering/lib/FPfilter/"
    shell:
        """
        FPfilter -v {input.unfiltered} -p {params.FP}
        """

rule rm_header_from_vcf:
    input:
        os.path.join(intermediate_path, "{filename}.Varscan.indel.Somatic.vcf_TP_after_filtered_sorted.vcf")
    output:
        os.path.join(intermediate_path, "{filename}_Varscan_indel_variants_noheader.vcf")
    shell:
        """
        egrep -v "^#" {input} > {output}
        """

rule run_vep:
    input:
        os.path.join(intermediate_path, "{filename}_Varscan_indel_variants_noheader.vcf"),
        fa=FASTA
    output:
        os.path.join(intermediate_path, "{filename}_Varscan_indel_vep.vcf")
    shell:
        """
        vep -i {input} --format vcf --check_ref --dont_skip --fasta {input.fa} --species mus_musculus --cache --offline --vcf -o {output} --force_overwrite  --symbol --plugin Wildtype --terms SO --plugin Downstream
        """

rule generate_fasta:
    input:
        os.path.join(intermediate_path, "{filename}_Varscan_indel_vep.vcf")
    output:
        len_15 = os.path.join(peptide_path, "{filename}_VarScan_indel_vep.15.peptide"),
        len_17 = os.path.join(peptide_path, "{filename}_VarScan_indel_vep.17.peptide"),
        len_19 = os.path.join(peptide_path, "{filename}_VarScan_indel_vep.19.peptide"),
        len_21 = os.path.join(peptide_path, "{filename}_VarScan_indel_vep.21.peptide"),
    shell:
       """
        pvacseq generate_protein_fasta {input} 15 {output.len_15};
        pvacseq generate_protein_fasta {input} 17 {output.len_17};
        pvacseq generate_protein_fasta {input} 19 {output.len_19};
        pvacseq generate_protein_fasta {input} 21 {output.len_21};
        """
