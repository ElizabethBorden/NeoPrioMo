# Setting up filenames here:
from os.path import join

# Files
#sample = ["3466","3702", "3713", "3919", "3903", "3703h2.1", "3879h2.1", "3926h2.1", "3879h1.1", "3926h1.1", "3703h11.1"]
sample = ["fusion"]

# Path to files
#path = "/scratch/eknodel/cohen_melanoma/validated_peptides/"
path = "/scratch/eknodel/cohen_melanoma/RNAseq/fusions_arriba/"

rule all:
    input:
        expand(path + "{sample}/{sample}_lukszablast.in", sample=sample, path=path),
        expand(path + "{sample}/{sample}_lukszablast.out", sample=sample, path=path),
        expand(path + "{sample}/neoantigen_fitness_{sample}.txt", sample=sample, path=path)

rule format_blast_input:
    input:
        os.path.join(path, "{sample}/{sample}_lukszaprogram.in")
    output:
        os.path.join(path, "{sample}/{sample}_lukszablast.in")
    shell:
        """
        cat {input} | awk '{{print ">sample|"$1"|MUT|"$2 "\\n" $5 "\\n" ">sample|"$1"|WT|"$2 "\\n" $4}}' > {output};
        sed -i -e 1,4d {output} 
        """

rule blast:
    input:
        os.path.join(path, "{sample}/{sample}_lukszablast.in")
    output:
        os.path.join(path, "{sample}/{sample}_lukszablast.out")
    shell:
        """
        blastp -query {input} -db ~/Cohen_melanoma/luksza_program/iedb.fasta -outfmt 5 -evalue 100000000  -gapopen 11 -gapextend 1 > {output}
        """

rule luksza:
    input:
        blast = os.path.join(path, "{sample}/{sample}_lukszablast.out"),
        neoantigens = os.path.join(path, "{sample}/{sample}_lukszaprogram.in")
    output:
        os.path.join(path, "{sample}/neoantigen_fitness_{sample}.txt")
    params:
        directory = "/scratch/eknodel/cohen_melanoma/validated_peptides/{sample}/"
    shell:
        """
        cd ~/Luksza_programs;
        python src/main.py {input.neoantigens} {params.directory} 26 4.86936 {output} {input.blast}
        """
