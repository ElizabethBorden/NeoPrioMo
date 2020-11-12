# Setting up filenames here:
from os.path import join

# Files
sample = ["3466"]

# Path to files
peptide_path = "/scratch/eknodel/cohen_melanoma/validated_peptides/"

rule all:
    input:
        expand(peptide_path + "{sample}/{sample}_9_netmhc.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_10_netmhc.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_netMHC.out", sample=sample)

rule netMHC:
    input:
        peptides = os.path.join(peptide_path, "{sample}/{sample}_validated_peptides.in")
    output:
        netMHC_9  = os.path.join(peptide_path, "{sample}/{sample}_9_netmhc.xsl"),
        netMHC_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netmhc.xsl")
    params:
        hla = "HLA-A01:01"
    shell:
        """
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -xls -l 9  -xlsfile {output.netMHC_9};
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -xls -l 10 -xlsfile {output.netMHC_10};
        """

rule combine_netMHC:
    input:
        netMHC_9 = os.path.join(peptide_path, "{sample}/{sample}_9_netmhc.xsl"),
        netMHC_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netmhc.xsl")
    output:
        netMHC = os.path.join(peptide_path, "{sample}/{sample}_netMHC.out")
    shell:
        """
        cat {input.netMHC_10} | sed '/Allele/d' | awk {{'print $2, $5, $6, $7, $8, $9}}' >> {output.netMHC};
        cat {input.netMHC_9} | sed '/Allele/d' | awk {{'print $2, $5, $6, $7, $8, $9}}' >> {output.netMHC} 
        """
