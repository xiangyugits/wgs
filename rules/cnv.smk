SAMPLES, = glob_wildcards('/path/to/BAMs/{sample}_merged_markDupl_BQSR.bam')
CHRN = list(range(1, 22))
CHRN.append('X','Y')
CHR = CHRN
MAP = 'canFam3_mappability_150.merged.bed.gz'
SEGDUP = 'segmental_duplication.bed.gz'


rule all:
    input:
        expand('chr{j}_{sample}_intervals_cohort.vcf.gz', j=CHR, sample=SAMPLES),
        expand('chr{j}_{sample}_segments_cohort.vcf.gz', j=CHR, sample=SAMPLES)


rule make_intervals:
    input:
        REF="resources/genome.fasta",
    params:
        'chr{j}'
    output:
        'resources/interval_chr{j}.interval_list'
    cache: True
    shell:
        '''
        gatk --java-options "-Xmx8G" PreprocessIntervals \
        -R {input} \
        --padding 0 \
        -L {params} \
        -imr OVERLAPPING_ONLY \
        -O {output}
        '''

rule annotate:
    input:
        ref = "resources/genome.fasta",
        interval = 'interval_chr{j}.interval_list',
        mappability = MAP,
        segduplication = SEGDUP
    output:
        'resources/annotated_intervals_chr{j}.tsv'
    cache:True
    shell:
        '''
        gatk --java-options "-Xmx8G" AnnotateIntervals \
        -R {input.ref} \
        -L {input.interval} \
        --mappability-track {input.mappability} \
        --segmental-duplication-track {input.segduplication} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O {output}
        '''

rule count_reads:
    input:
        ref = "resources/genome.fasta",
        bam = get_sample_bams,
        interval = 'resources/interval_chr{j}.interval_list'
    output:
        'cnv/{sample}_chr{j}.hdf5'
    shell:
        '''
        gatk --java-options "-Xmx8G" CollectReadCounts \
        -R {input.ref} \
        -imr OVERLAPPING_ONLY \
        -L {input.interval} \
        -I {input.bam} \
        -O {output}
        '''
rule filter_intervals:
    input:
        interval = 'resources/interval_chr{j}.interval_list',
        annotated = 'resources/annotated_intervals_chr{j}.tsv',
        samples = expand('{sample}_{chromosome}.hdf5', sample=SAMPLES, chromosome='chr{j}'),
    output:
        'cnv/gcfiltered_chr{j}.interval_list'
    params:
        files = lambda wildcards, input: ' -I '.join(input.samples)
    shell:
        '''
        gatk --java-options "-Xmx8G" FilterIntervals \
        -L {input.interval} \
        --annotated-intervals {input.annotated} \
        -I {params.files} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O {output}
        '''


#This step is needed to generate global baseline coverage and noise data for the subsequent steps:
rule determine_ploidy:
    input:
        interval = 'cnv/gcfiltered_chr{j}.interval_list',
        samples = expand('cnv/{sample}_{chromosome}.hdf5', sample=SAMPLES, chromosome='chr{j}'),
        prior = 'ploidy_priors.tsv',
    params:
        prefix = 'dogs',#!!!!!!!!!!!!#
        files = lambda wildcards, input: ' -I '.join(input.samples)
    output:
        'ploidy-calls_chr{j}'
    shell:
        '''
        gatk --java-options "-Xmx8G" DetermineGermlineContigPloidy \
        -L {input.interval} \
        -I {params.files} \
        --contig-ploidy-priors {input.prior} \
        --output-prefix  {params.prefix} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O {output}
        '''

rule scattering:
    input:
        interval = 'cnv/gcfiltered_chr{j}.interval_list'
    output:
        dynamic('cnv/scatter_chr{j}/{fragment}/scattered.interval_list')
    params:
        'scatter_chr{j}'
    shell:
        '''
        mkdir -p {params} # needed because Snakemake fails creating this directory automatically
        gatk --java-options "-Xmx8G" IntervalListTools \
        --INPUT {input.interval} \
        --SUBDIVISION_MODE INTERVAL_COUNT \
        --SCATTER_CONTENT 15000 \
        --OUTPUT {params}
        '''

rule cnvcall:
    input:
        interval = 'cnv/scatter_chr{j}/{fragment}/scattered.interval_list',
        sample = expand("{sample}_{chromosome}.hdf5", sample=SAMPLES, chromosome='chr{j}'),
        annotated = 'annotated_intervals_chr{j}.tsv',
        ploidy = 'ploidy-calls_chr{j}'
    output:
        modelf = "cohort-calls_chr{j}/frag_{fragment}-model",
        callsf = "cohort-calls_chr{j}/frag_{fragment}-calls"
    params:
        outdir = 'cohort-calls_chr{j}',
        outpref = 'frag_{fragment}',
        files = lambda wildcards, input: " -I ".join(input.sample)
    shell:
        '''
        gatk --java-options "-Xmx8G" GermlineCNVCaller  \
        --run-mode COHORT \
        -L {input.interval} \
        -I {params.files} \
        --contig-ploidy-calls {input.ploidy}/dogs-calls \
        --annotated-intervals {input.annotated} \
        --output-prefix {params.outpref} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O {params.outdir}
        '''

def sampleindex(sample):
    index = SAMPLES.index(sample)
    return index

rule process_cnvcalls:
    input:
        model = dynamic("cohort-calls_chr{j}/frag_{fragment}-model"),
        calls = dynamic("cohort-calls_chr{j}/frag_{fragment}-calls"),
        dict  = "resources/genome.dict",
        ploidy = 'ploidy-calls_chr{j}'
    output:
        intervals = 'chr{j}_{sample}_intervals_cohort.vcf.gz',
        segments = 'chr{j}_{sample}_segments_cohort.vcf.gz'
    params:
        index = lambda wildcards: sampleindex(wildcards.sample),
        modelfiles = lambda wildcards, input: " --model-shard-path ".join(input.model),
        callsfiles = lambda wildcards, input: " --calls-shard-path ".join(input.calls)
    shell:
        '''
        gatk --java-options "-Xmx8G" PostprocessGermlineCNVCalls \
        --model-shard-path {params.modelfiles} \
        --calls-shard-path {params.callsfiles} \
        --sequence-dictionary {input.dict} \
        --allosomal-contig chrX \
        --contig-ploidy-calls {input.ploidy}/dogs-calls \
        --sample-index {params.index} \
        --output-genotyped-intervals  {output.intervals} \
        --output-genotyped-segments  {output.segments}
        '''