# Setting up filesnames here:
#sample=["3919"]
sample = ["fusion"]

# path to where bams are
blast_path = "/scratch/eknodel/cohen_melanoma/RNAseq/fusions_arriba/"

rule all:
    input:
        expand(blast_path + "{sample}/{sample}.epitopes.annotated.tsv", blast_path=blast_path, sample=sample),
        expand(blast_path + "{sample}/{sample}.WTmatch.in", blast_path=blast_path, sample=sample),
        expand(blast_path + "{sample}/{sample}.WTmatch_9_netCTL.out", blast_path=blast_path, sample=sample),
        expand(blast_path + "{sample}/{sample}.WTmatch_10_netCTL.out", blast_path=blast_path, sample=sample),
        expand(blast_path + "{sample}/{sample}.WTmatch_netCTL.out", blast_path=blast_path, sample=sample)

rule filter_blast:
    input:
        os.path.join(blast_path, "{sample}/{sample}_peptides.in")
    output:
        os.path.join(blast_path, "{sample}/{sample}.epitopes.annotated.tsv")
    params:
        sample = "{sample}",
        outdir = os.path.join(blast_path,"{sample}")
    shell:
        """
        python sequence_similarity.py -o {params.outdir} -s {params.sample} -i {input}
        """

rule prep_WT_netMHC:
    input:
        os.path.join(blast_path, "{sample}/{sample}.epitopes.annotated.tsv")
    output:
        os.path.join(blast_path, "{sample}/{sample}.WTmatch.in")
    shell:
        """
        cat {input} | awk '{{print ">"$1 "\\n" $2}}' > {output}
        """

rule run_WT_netMHC:
    input:
        peptides = os.path.join(blast_path, "{sample}/{sample}.WTmatch.in")
    output:
        netCTL_9 = os.path.join(blast_path, "{sample}/{sample}.WTmatch_9_netCTL.out"),
        netCTL_10 = os.path.join(blast_path, "{sample}/{sample}.WTmatch_10_netCTL.out")
    params:
        hla = "HLA-A01:01"
    shell:
        """
        netCTLpan -a {params.hla} -f {input.peptides} -s -l 9 -xls -xlsfile {output.netCTL_9};
        netCTLpan -a {params.hla} -f {input.peptides} -s -l 10 -xls -xlsfile {output.netCTL_10}
        """

rule combine_outputs:
    input:
        netCTL_9 = os.path.join(blast_path, "{sample}/{sample}.WTmatch_9_netCTL.out"),
        netCTL_10 = os.path.join(blast_path, "{sample}/{sample}.WTmatch_10_netCTL.out")
    output:
        netCTL = os.path.join(blast_path, "{sample}/{sample}.WTmatch_netCTL.out")
    shell:
        """
        echo "Peptide  WTmatch MHC TAP Cle Comb    %Rank" >> {output.netCTL};
        cat {input.netCTL_10} | sed '/Allele/d' | awk {{'print $2, $3, $5, $6, $7, $8, $9}}' >> {output.netCTL};
        cat {input.netCTL_9} | sed '/Allele/d' | awk {{'print $2, $3, $5, $6, $7, $8, $9}}' >> {output.netCTL}
        """
