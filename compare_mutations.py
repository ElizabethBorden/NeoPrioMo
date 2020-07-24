#### Input files ####
gatk = set(line.strip() for line in open('/data/storage/DATASETS/DOWNLOADS/MELANOMA/peptides/YM-1_Gatk_vep.21.collapsed'))
varscan = set(line.strip() for line in open('/data/storage/DATASETS/DOWNLOADS/MELANOMA/peptides/YM-1_VarScan_vep.21.collapsed'))
strelka = set(line.strip() for line in open('/data/storage/DATASETS/DOWNLOADS/MELANOMA/peptides_strelka/YM-1_Strelka_vep.21.collapsed'))
varscan_indels = set(line.strip() for line in open('/data/storage/DATASETS/DOWNLOADS/MELANOMA/peptides/YM-1_VarScan_indel_vep.21.collapsed'))

#### Counters #####
a=0
b=0
c=0
d=0

#### Get data for Venn diagram ####
with open ('venn_diagram_numbers.out', 'a') as file:
    for line in gatk:
        if line in varscan:
            a=a+1
    file.write('gatk_varscan:')
    file.write(str(a))
    file.write('\n')
    for line in strelka:
        if line in varscan:
            b=b+1
    file.write('strelka_varscan:')
    file.write(str(b))
    file.write('\n')
    for line in strelka:
        if line in gatk:
            c=c+1
    file.write('strelka_gatk:')
    file.write(str(c))
    file.write('\n') 
    for line in strelka:
        if line in gatk & varscan:
            d=d+1
    file.write('strelka_gatk_varscan:')
    file.write(str(d))
    file.write('\n')

a=0
#### Get data for Venn diagram indels####
with open ('indel_venn_diagram_numbers.out', 'a') as file:
    for line in gatk:
        if line in varscan_indels:
            a=a+1
    file.write('gatk_varscan:')
    file.write(str(a))
    file.write('\n')
with open ('indel_gatk_varscan.out', 'a') as file:
    for line in gatk:
        if line in varscan_indels:
            file.write(line)
            file.write('\n')

#### Summarize epitopes from any two programs ####
with open ('summary_anytwo.out', 'a') as file:
        for line in strelka:
            if line in gatk:
                file.write(line)
                file.write('\n')
            else:
                if line in varscan:
                    file.write(line)
                    file.write('\n')
        for line in gatk:
            if line in strelka:
                print('Already accounted for - skipping')
            else:
                if line in varscan:
                    file.write(line)
                    file.write('\n')

