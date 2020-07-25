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
VarScan_path = "/scratch/external_scripts/VarScan.v2.3.9.jar"
perl_filter_path = "/scratch/external_scripts/fpfilter-tool-master/fpfilter-2.pl"
perl5lib_path = "/packages/6x/vcftools/0.1.12b/lib/per15/site_perl/"

rule all:
    input:
        expand(pileup_path + "{sample}.pileup", pileup_path=pileup_path, sample=SAMPLES),
        expand(bams_path + "{sample}_WES_HISAT2_aligned_sorted.bam", bams_path=bams_path, sample=SAMPLES),
        expand(bai_path + "{sample}_WES_HISAT2_aligned_sorted.bam.bai", bai_path=bai_path, sample=SAMPLES),
        expand(intermediate_path +"{filename}.Varscan.snp",intermediate_path=intermediate_path, filename=tumor),
        expand(intermediate_path +"{filename}.Varscan.indel", intermediate_path=intermediate_path, filename=tumor),
        expand(intermediate_path +"{filename}.Varscan.snp.Somatic.hc", intermediate_path=intermediate_path,  filename=tumor),
        expand(intermediate_path +"{filename}.Varscan.snp.Somatic.hc.filter", intermediate_path=intermediate_path, filename=tumor),
        expand(intermediate_path +"{filename}.Varscan.snp.Somatic.hc.filter.bed", intermediate_path=intermediate_path, filename=tumor),
        expand(intermediate_path +"{filename}.readcounts", intermediate_path=intermediate_path, filename=tumor),
        expand(intermediate_path +"{filename}_Varscan_variants_filter.pass", intermediate_path=intermediate_path, filename=tumor),
        expand(intermediate_path +"{filename}_Varscan_variants_filter.pass.vcf", intermediate_path=intermediate_path, filename=tumor), 
        expand(intermediate_path +"{filename}_Varscan_variants_filter.pass.noheader.vcf",intermediate_path=intermediate_path, filename=tumor),
        expand(intermediate_path +"{filename}_Varscan_vep.vcf",intermediate_path=intermediate_path, filename=tumor),
        expand(intermediate_path +"{filename}_Varscan_vep.tsv", intermediate_path=intermediate_path, filename=tumor), 
        expand(peptide_path +"{filename}_VarScan_vep.15.peptide", intermediate_path=intermediate_path, filename=tumor),
        expand(peptide_path +"{filename}_VarScan_vep.17.peptide", intermediate_path=intermediate_path, filename=tumor),
        expand(peptide_path +"{filename}_VarScan_vep.19.peptide", intermediate_path=intermediate_path, filename=tumor),
        expand(peptide_path +"{filename}_VarScan_vep.21.peptide", intermediate_path=intermediate_path, filename=tumor)

rule sort_bam:
    input:
        bam = os.path.join(bams_path, "{sample}_WES_HISAT2_aligned.bam")
    output:
        bam = os.path.join(bams_path, "{sample}_WES_HISAT2_aligned_sorted.bam")
    shell:
        """
        samtools sort {input.bam} -o {output.bam}
        """

rule bam_pileup:
    input:
        fa = FASTA,
        bam = os.path.join(bams_path, "{sample}_WES_HISAT2_aligned_sorted.bam")
    output:
        pileup = os.path.join(pileup_path, "{sample}.pileup")
    threads: 6
    shell:
        """
        samtools mpileup -f {input.fa} {input.bam} > {output.pileup}
        """

rule index_bam:
    input:
        bam = os.path.join(bams_path, "{sample}_WES_HISAT2_aligned_sorted.bam")
    output:
        bai = os.path.join(bai_path, "{sample}_WES_HISAT2_aligned_sorted.bam.bai")
    shell:
        """
        samtools index {input.bam} {output.bai}
        """

rule run_VarScan:
    input:
        pileup = expand("/data/storage/DATASETS/DOWNLOADS/MELANOMA/pileup/{sample}.pileup", sample=SAMPLES)
    output:
        snp = os.path.join(intermediate_path, "{filename}.Varscan.snp"),
        indel = os.path.join(intermediate_path, "{filename}.Varscan.indel")
    params:
        VarScan = VarScan_path,
        basename = os.path.join(intermediate_path, "{filename}.Varscan")
    threads: 6
    shell:
        "java -jar {params.VarScan} somatic {input.pileup} {params.basename} min-coverage 10 min-var-freq 0.08 somatic-p-value 0.05"

rule isolate_calls_by_type_and_confidence:
    input:
        VarScan_snp = os.path.join(intermediate_path, "{filename}.Varscan.snp") 
    output:
        VarScan_snp = os.path.join(intermediate_path, "{filename}.Varscan.snp.Somatic.hc")
    params:
        VarScan = VarScan_path
    shell:
        """
        java -jar {params.VarScan} processSomatic {input.VarScan_snp}
        """
rule somatic_filter:
    input:
        snp_somatic_hc = os.path.join(intermediate_path, "{filename}.Varscan.snp.Somatic.hc"),
        indel = os.path.join(intermediate_path, "{filename}.Varscan.indel")
    output:
        snp_somatic_hc_filter = os.path.join(intermediate_path, "{filename}.Varscan.snp.Somatic.hc.filter")
    params:
        VarScan = VarScan_path
    shell:
        """
        java -jar {params.VarScan} somaticFilter {input.snp_somatic_hc} -indel-file {input.indel} -output-file {output.snp_somatic_hc_filter}
        """

rule convert_to_bed_fmt:
    input:
        snp_somatic_hc_filter = os.path.join(intermediate_path, "{filename}.Varscan.snp.Somatic.hc.filter")
    output:
        snp_somatic_hc_filter_bed = os.path.join(intermediate_path, "{filename}.Varscan.snp.Somatic.hc.filter.bed")
    shell:
        """
        awk -F "\t" '{{print $1 "\t" $2 "\t" $2 }}' {input.snp_somatic_hc_filter} | tail -n+2 > {output.snp_somatic_hc_filter_bed}
        """

rule readcount:
    input:
        fa = FASTA,
        snp_somatic_hc_filter_bed = os.path.join(intermediate_path, "{filename}.Varscan.snp.Somatic.hc.filter.bed"),
        tumor = expand(os.path.join(bams_path, "{tumor}_WES_HISAT2_aligned_sorted.bam"), tumor=tumor),
        bai = expand(os.path.join(bai_path, "{sample}_WES_HISAT2_aligned_sorted.bam.bai"), sample=SAMPLES)
    output:
        readcounts = os.path.join(intermediate_path, "{filename}.readcounts")
    params:
        bamreadcount = bamreadcount_path
    shell:
        """
        {params.bamreadcount} -q 1 -b 20 -f {input.fa} -l {input.snp_somatic_hc_filter_bed} {input.tumor} > {output.readcounts}
        """

rule perl_filter:
    input:
        snp_somatic_hc_filter = os.path.join(intermediate_path, "{filename}.Varscan.snp.Somatic.hc.filter"),
        readcounts = os.path.join(intermediate_path, "{filename}.readcounts")
    output:
        out = os.path.join(intermediate_path, "{filename}_Varscan_variants_filter.pass")
    params:
        basename = os.path.join(intermediate_path, "{filename}_Varscan_variants_filter"),
        perlfilter = perl_filter_path
    shell:
        """
        perl {params.perlfilter} {input.snp_somatic_hc_filter} {input.readcounts} --output-basename {params.basename}
        """

rule convert_to_vcf:
    input:
        os.path.join(intermediate_path, "{filename}_Varscan_variants_filter.pass")
    output:
        os.path.join(intermediate_path, "{filename}_Varscan_variants_filter.pass.vcf")
    shell:
        """
        python ~/external_scripts/vscan2vcf.py {input} > {output}
        """

rule rm_header_from_vcf:
    input:
        os.path.join(intermediate_path, "{filename}_Varscan_variants_filter.pass.vcf")
    output:
        os.path.join(intermediate_path, "{filename}_Varscan_variants_filter.pass.noheader.vcf")
    shell:
        """
        egrep -v "^#" {input} > {output}
        """

rule run_vep:
    input:
        os.path.join(intermediate_path, "{filename}_Varscan_variants_filter.pass.noheader.vcf")
    output:
        os.path.join(intermediate_path, "{filename}_Varscan_vep.vcf")
    params:
        #vep = vep_path,
        #plugins = plugin_path,
        #perl5lib = perl5lib_path
    shell:
        """
        vep -i {input} --format vcf --species mus_musculus --cache --offline --vcf -o {output} --force_overwrite  --symbol --plugin Wildtype --terms SO --plugin Downstream
        """

rule pvacseq_convert_vcf:
    input:
        os.path.join(intermediate_path, "{filename}_Varscan_vep.vcf")
    output:
        os.path.join(intermediate_path, "{filename}_Varscan_vep.tsv")
    shell:
        """
        pvacseq convert_vcf {input} {output}
        """

rule generate_fasta:
    input:
        os.path.join(intermediate_path, "{filename}_Varscan_vep.tsv")
    output:
        len_15 = os.path.join(peptide_path, "{filename}_VarScan_vep.15.peptide"),
        len_17 = os.path.join(peptide_path, "{filename}_VarScan_vep.17.peptide"),
        len_19 = os.path.join(peptide_path, "{filename}_VarScan_vep.19.peptide"),
        len_21 = os.path.join(peptide_path, "{filename}_VarScan_vep.21.peptide"),
    shell:
       """
        pvacseq generate_fasta {input} 15 {output.len_15};
        pvacseq generate_fasta {input} 17 {output.len_17};
        pvacseq generate_fasta {input} 19 {output.len_19};
        pvacseq generate_fasta {input} 21 {output.len_21};
        """
