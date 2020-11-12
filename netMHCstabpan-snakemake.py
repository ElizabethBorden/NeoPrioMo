# Setting up filenames here:
from os.path import join

# Files
sample = ["3466"]

# Path to files
peptide_path = "/scratch/eknodel/cohen_melanoma/validated_peptides/"

rule all:
    input:
        expand(peptide_path + "{sample}/{sample}_9_netmhcstab.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_10_netmhcstab.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_netMHCstab.out", sample=sample)

rule netMHCstab:
    input:
        peptides = os.path.join(peptide_path, "{sample}/{sample}_validated_peptides.in")
    output:
        netMHCstab_9  = os.path.join(peptide_path, "{sample}/{sample}_9_netmhcstab.xsl"),
        netMHCstab_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netmhcstab.xsl")
    params:
        hla = "HLA-A01:01"
    shell:
        """
        netMHCstabpan -a {params.hla} -f {input.peptides} -s -xls -l 9  -xlsfile {output.netMHCstab_9};
        netMHCstabpan -a {params.hla} -f {input.peptides} -s -xls -l 10 -xlsfile {output.netMHCstab_10};
        """

rule combine_netMHCstab:
    input:
        netMHCstab_9 = os.path.join(peptide_path, "{sample}/{sample}_9_netmhcstab.xsl"),
        netMHCstab_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netmhcstab.xsl")
    output:
        netMHCstab = os.path.join(peptide_path, "{sample}/{sample}_netMHCstab.out")
    shell:
        """
        cat {input.netMHCstab_10} | sed '/Allele/d' | awk {{'print $2, $5, $6, $7, $8, $9}}' >> {output.netMHCstab};
        cat {input.netMHCstab_9} | sed '/Allele/d' | awk {{'print $2, $5, $6, $7, $8, $9}}' >> {output.netMHCstab} 
        """
