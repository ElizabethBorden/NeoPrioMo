# Setting up filesnames here:
from os.path import join
#configfile: "Melanoma.config.json"

#Reference File
FASTA = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/reference/my_dir/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna",

results_path="Melanoma_somatic/results/variants/"
bam_path="/data/storage/DATASETS/DOWNLOADS/MELANOMA/sorted_bam/"
intermediate_path="/data/storage/DATASETS/DOWNLOADS/MELANOMA/intermediates_strelka/"
peptide_path = "/data/storage/DATASETS/DOWNLOADS/MELANOMA/peptides_strelka/"

bamreadcount_path = "bam-readcount"
perl_filter_path = "~/external_scripts/fpfilter-tool-master/fpfilter-2.pl"
perl5lib_path = "/packages/6x/vcftools/0.1.12b/lib/per15/site_perl/"

rule all:
    input:
        os.path.join(results_path, "somatic.snvs.vcf.gz"),
        os.path.join(results_path, "somatic.indels.vcf.gz"),
        "Melanoma_somatic/runWorkflow.py",
        os.path.join(results_path, "somatic.snvs.vcf"),
        os.path.join(results_path, "somatic.indels.vcf"),
        os.path.join(results_path, "somatic.snvs_PASS.vcf"),
        os.path.join(results_path, "somatic.indels_PASS.vcf"),
        os.path.join(results_path, "somatic.snvs_PASS.vcf_TP_after_filtered_sorted.vcf"),
        os.path.join(results_path, "somatic.indels_PASS.vcf_TP_after_filtered_sorted.vcf"),
        #os.path.join(intermediate_path, "YM-1_snvs.noheader.vcf"),
        #os.path.join(intermediate_path, "YM-1_indels.noheader.vcf"),
        os.path.join(intermediate_path, "YM-1_snvs_vep.vcf"),
        os.path.join(intermediate_path, "YM-1_indels_vep.vcf"),
        #os.path.join(intermediate_path, "YM-1_snvs_vep.tsv"),
        #os.path.join(intermediate_path, "YM-1_indels_vep.tsv"),
        os.path.join(peptide_path, "YM-1_snvs_vep.15.peptide"),
        os.path.join(peptide_path, "YM-1_snvs_vep.17.peptide"),
        os.path.join(peptide_path, "YM-1_snvs_vep.19.peptide"),
        os.path.join(peptide_path, "YM-1_snvs_vep.21.peptide"),
        os.path.join(peptide_path, "YM-1_indels_vep.15.peptide"),
        os.path.join(peptide_path, "YM-1_indels_vep.17.peptide"),
        os.path.join(peptide_path, "YM-1_indels_vep.19.peptide"),
        os.path.join(peptide_path, "YM-1_indels_vep.21.peptide")

rule config:
    input: 
        normal = os.path.join(bam_path, "WT-1_WES_HISAT2_aligned_sortedbycoord.bam"),
        tumor = os.path.join(bam_path, "YM-1_WES_HISAT2_aligned_sortedbycoord.bam"),
        fa = FASTA
    output:
        "Melanoma_somatic/runWorkflow.py"
    params:     
        strelka = "~/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py"
    shell:
        """
        {params.strelka} --normalBam {input.normal} --tumorBam {input.tumor} --referenceFasta {input.fa} --runDir Melanoma_somatic
        """
rule run:
    input:
        "Melanoma_somatic/runWorkflow.py"
    output: 
        snvs = os.path.join(results_path, "somatic.snvs.vcf.gz"),
        indels = os.path.join(results_path, "somatic.indels.vcf.gz")
    params:
        run = "~/Melanoma_dataset/Mutation_calling_strelka/Melanoma_somatic/runWorkflow.py"
    shell:
        """
        {params.run} -m local -j 20;
        """

rule gunzip:
    input:
        gz_snvs = os.path.join(results_path, "somatic.snvs.vcf.gz"),
        gz_indels = os.path.join(results_path, "somatic.indels.vcf.gz")
    output:
        snvs = os.path.join(results_path, "somatic.snvs.vcf"),
        indels = os.path.join(results_path, "somatic.indels.vcf")
    shell:
        """
        gunzip -c {input.gz_snvs} > {output.snvs};
        gunzip -c {input.gz_indels} > {output.indels}
        """
rule grep_pass:
    input:
        snvs = os.path.join(results_path, "somatic.snvs.vcf"),
        indels = os.path.join(results_path, "somatic.indels.vcf")
    output:
        snvs = os.path.join(results_path, "somatic.snvs_PASS.vcf"),
        indels = os.path.join(results_path, "somatic.indels_PASS.vcf")
    shell:
        """
        grep "PASS" {input.snvs} > {output.snvs};
        grep "PASS" {input.indels} > {output.indels}
        """

rule FP_filter:
    input:
        fa = FASTA,
        unfiltered = os.path.join(results_path, "somatic.indels_PASS.vcf"),
        unfiltered_snvs = os.path.join(results_path, "somatic.snvs_PASS.vcf")
    output:
        filtered= os.path.join(results_path, "somatic.indels_PASS.vcf_TP_after_filtered_sorted.vcf"),
        filtered_snvs= os.path.join(results_path, "somatic.snvs_PASS.vcf_TP_after_filtered_sorted.vcf")
    params:
        FP = "~/miniconda3/envs/Filtering/lib/FPfilter/"
    shell:
        """
        FPfilter -v {input.unfiltered} -p {params.FP};
        FPfilter -v {input.unfiltered_snvs} -p {params.FP}
        """

rule run_vep:
    input:
        snvs = os.path.join(results_path, "somatic.snvs_PASS.vcf_TP_after_filtered_sorted.vcf"),
        indels = os.path.join(results_path, "somatic.indels_PASS.vcf_TP_after_filtered_sorted.vcf")
    output:
        snvs = os.path.join(intermediate_path, "YM-1_snvs_vep.vcf"),
        indels = os.path.join(intermediate_path, "YM-1_indels_vep.vcf")
    shell:
        """
        vep -i {input.snvs} --format vcf --species mus_musculus --cache --offline --vcf -o {output.snvs} --force_overwrite  --symbol --plugin Wildtype --terms SO --plugin Downstream;
        vep -i {input.indels} --format vcf --species mus_musculus --cache --offline --vcf -o {output.indels} --force_overwrite  --symbol --plugin Wildtype --terms SO --plugin Downstream
        """

rule generate_fasta:
    input:
        snvs = os.path.join(intermediate_path, "YM-1_snvs_vep.vcf"),
        indels = os.path.join(intermediate_path, "YM-1_indels_vep.vcf")
    output:
        snvs_len_15 = os.path.join(peptide_path, "YM-1_snvs_vep.15.peptide"),
        snvs_len_17 = os.path.join(peptide_path, "YM-1_snvs_vep.17.peptide"),
        snvs_len_19 = os.path.join(peptide_path, "YM-1_snvs_vep.19.peptide"),
        snvs_len_21 = os.path.join(peptide_path, "YM-1_snvs_vep.21.peptide"),
        indels_len_15 = os.path.join(peptide_path, "YM-1_indels_vep.15.peptide"),
        indels_len_17 = os.path.join(peptide_path, "YM-1_indels_vep.17.peptide"),
        indels_len_19 = os.path.join(peptide_path, "YM-1_indels_vep.19.peptide"),
        indels_len_21 = os.path.join(peptide_path, "YM-1_indels_vep.21.peptide")
    shell:
        """
        pvacseq generate_protein_fasta {input.snvs} 15 {output.snvs_len_15};
        pvacseq generate_protein_fasta {input.snvs} 17 {output.snvs_len_17};
        pvacseq generate_protein_fasta {input.snvs} 19 {output.snvs_len_19};
        pvacseq generate_protein_fasta {input.snvs} 21 {output.snvs_len_21};
        pvacseq generate_protein_fasta {input.indels} 15 {output.indels_len_15};
        pvacseq generate_protein_fasta {input.indels} 17 {output.indels_len_17};
        pvacseq generate_protein_fasta {input.indels} 19 {output.indels_len_19};
        pvacseq generate_protein_fasta {input.indels} 21 {output.indels_len_21};
        """
