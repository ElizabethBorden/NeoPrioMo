# Setting up filesnames here:
from os.path import join

# Files
sample = ["SRR3877300", "SRR3877303", "SRR3877307", "SRR3877308", "SRR3877311", "SRR3877314", "SRR3877317", "SRR3877297", "SRR3877299", "SRR3877302", "SRR3877305", "SRR3877310", "SRR3877319", "SRR3877321"]

# Path to files
output_path = "/scratch/eknodel/Tsai_analysis/Fusions_arriba/"

rule all:
    input:
        expand(output_path + "{sample}/{sample}_netctl.xsl", sample=sample),
        expand(output_path + "{sample}/{sample}_netctl.out", sample=sample),
        expand(output_path + "{sample}/{sample}_netctl_filtered.out", sample=sample)

rule netCTL:
    input:
        peptides = os.path.join(output_path, "{sample}/{sample}_peptides.out")
    output:
        netCTL = os.path.join(output_path, "{sample}/{sample}_netctl.xsl")
    shell:
        """
        netCTLpan -a HLA-A03:01,HLA-B14:02,HLA-B14:03,HLA-C08:02 -f {input.peptides} -s -xls -xlsfile {output.netCTL}
        """

rule remove_header:
    input:
        netCTL = os.path.join(output_path, "{sample}/{sample}_netctl.xsl")
    output:
        netCTL_noheader = os.path.join(output_path, "{sample}/{sample}_netctl.out")
    shell:
        """
        cat {input.netCTL} | sed '/Allele/d' > {output.netCTL_noheader}
        """

rule filter_netCTL:
    input: 
        netCTL_noheader = os.path.join(output_path, "{sample}/{sample}_netctl.out")
    output: 
        netCTL_filter = os.path.join(output_path, "{sample}/{sample}_netctl_filtered.out")
    shell:
        """
        python ~/Novel_model/filter_MHCrank.py {input.netCTL_noheader} {output.netCTL_filter}
        """
