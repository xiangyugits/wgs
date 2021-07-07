include: "rules/common.smk"


##### Target rules #####


rule all:
    input:
        expand("recal/{sample}-{unit}.bam",sample=units.loc[:,'sample'],unit=units.unit.drop_duplicates()),

        #"annotated/all.vcf.gz",
        #"qc/multiqc.html",
        #"plots/depths.svg",
        #"plots/allele-freqs.svg",



##### Modules #####
    

include: "rules/ref.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
