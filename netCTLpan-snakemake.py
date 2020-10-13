# Setting up filenames here:
from os.path import join

# Files
sample = ["fusion"]

# Path to files
peptide_path = "/scratch/eknodel/cohen_melanoma/RNAseq/fusions_arriba/"

rule all:
    input:
        expand(peptide_path + "{sample}/{sample}_9_netctl.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_10_netctl.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_netCTL.out", sample=sample)

rule netCTL:
    input:
        peptides = os.path.join(peptide_path, "{sample}/{sample}_peptides.in")
    output:
        netCTL_9  = os.path.join(peptide_path, "{sample}/{sample}_9_netctl.xsl"),
        netCTL_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netctl.xsl")
    params:
        hla = "HLA-A01:01"
    shell:
        """
        netCTLpan -a {params.hla} -f {input.peptides} -s -xls -l 9  -xlsfile {output.netCTL_9};
        netCTLpan -a {params.hla} -f {input.peptides} -s -xls -l 10 -xlsfile {output.netCTL_10};
        """

rule combine_netCTL:
    input:
        netCTL_9 = os.path.join(peptide_path, "{sample}/{sample}_9_netctl.xsl"),
        netCTL_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netctl.xsl")
    output:
        netCTL = os.path.join(peptide_path, "{sample}/{sample}_netCTL.out")
    shell:
        """
        echo "Peptide  WTmatch MHC TAP Cle Comb    %Rank" >> {output.netCTL};
        cat {input.netCTL_10} | sed '/Allele/d' | awk {{'print $2, $5, $6, $7, $8, $9}}' >> {output.netCTL};
        cat {input.netCTL_9} | sed '/Allele/d' | awk {{'print $2, $5, $6, $7, $8, $9}}' >> {output.netCTL} 
        """
