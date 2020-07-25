# Setting up filesnames here:
from os.path import join

# Files
sample = ["SRR3877300", "SRR3877303", "SRR3877307"]

# Path to files
fastq_path = "/data/CEM/wilsonlab/from_collaborators/Tsai_SCC/RNAseq/"
output_path = "/scratch/eknodel/Tsai_analysis/Fusions_arriba/"

rule all:
    input:
        expand(output_path + "{sample}/{sample}_fusions.out", sample=sample),
        expand(output_path + "{sample}/{sample}_fusions_pass.out", sample=sample),
        expand(output_path + "{sample}/{sample}_peptides.out",sample=sample)

rule sort_bam:
    input:
        fastq_1 = os.path.join(fastq_path, "{sample}_1.fastq.gz"),
        fastq_2 = os.path.join(fastq_path, "{sample}_2.fastq.gz")
    output:
        fusions = os.path.join(output_path, "{sample}/{sample}_fusions.out"),
        fusions_discarded = os.path.join(output_path, "{sample}/{sample}_fusions_discarded.out"),
        directory = os.path.join(output_path, "{sample}")
    params:
        arriba = "~/bin/arriba_v1.2.0/run_arriba.sh",
    shell:
        """
        cd {output.directory};
        {params.arriba} ~/bin/arriba_v1.2.0/STAR_index_GRCh38_GENCODE28 \
        ~/bin/arriba_v1.2.0/GENCODE28.gtf \
        ~/bin/arriba_v1.2.0/GRCh38.fa \
        ~/bin/arriba_v1.2.0/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz \
        {input.fastq_1} \
        {input.fastq_2} \
        8 \
        {output.fusions} \
        {output.fusions_discarded} \
        """

rule loose_filter:
    input:
        fusions = os.path.join(output_path, "{sample}/{sample}_fusions.out")
    output:
        fusions_filtered = os.path.join(output_path, "{sample}/{sample}_fusions_pass.out")
    shell:
        """
        grep "high" {input.fusions} > {output.fusions_filtered};
        grep "medium" {input.fusions} >> {output.fusions_filtered}
        """

rule generate_peptides:
    input:
        fusions = os.path.join(output_path, "{sample}/{sample}_fusions_pass.out")
    output:
        peptides = os.path.join(output_path, "{sample}/{sample}_peptides.out")
    shell:
        r"""
        cat {input.fusions} | awk '{{print $1"/"$2, $23}}' | sed 's/*//' | sed 's/|//' | sed 's/^/>/' | sed '/\./d' | sed 's/ /\n/' > {output.peptides}
        """
