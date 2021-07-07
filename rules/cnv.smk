SAMPLES = glob_wildcards('recal/{sample}.bam')
CHRN = list(range(1, 22))
CHRN.append('X','Y')
CHR = CHRN
MAP = 'canFam3_mappability_150.merged.bed.gz'
SEGDUP = 'segmental_duplication.bed.gz'


rule all:
    input:
        expand('chr{j}_{sample}_intervals_cohort.vcf.gz', j=CHR, sample=SAMPLES),
        expand('chr{j}_{sample}_segments_cohort.vcf.gz', j=CHR, sample=SAMPLES)


# 对bins进行前期处理以用来计算reads coverage，
# 首先检查输入的interval是否有overlap，有则合并；然后根据指定参数扩充interval，分成bins,按指定bin长切割bins，最后过滤掉都是N的bins
#对WGS数据分析，需要使用窗口分割（--bin-length）
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

#注释每个目标区段，GC、mappability等信息
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

#计算指定的intervals的reads数，即计算起始位点落入intervals的reads数。
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

#对AnnotateIntervals生成的带注释的intervals和/或CollectReadCounts输出的reads计数信息(sample*.counts.hdf5)，
#指定过滤区域，输出一个过滤后的Picard intervals列表。
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

#给定CollectReadCounts生成的HDF5或TSV计数文件，确定种系样品的染色体倍数。
#This step is needed to generate global baseline coverage and noise data for the subsequent steps:
#使用reads count进行种系核型分析需要校准（“建模”）各染色体的coverage bias和variantion。
rule determine_ploidy:
    input:
        interval = 'cnv/gcfiltered_chr{j}.interval_list',
        samples = expand('cnv/{sample}_{chromosome}.hdf5', sample=SAMPLES, chromosome='chr{j}'),
        prior = 'ploidy_priors.tsv',
    params:
        prefix = 'prefix',#!!!!!!!!!!!!#
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

# 根据样本的CollectReadsCounts生成的reads count计数HDF5/TSV文件和DetermineGermlineContigPloidy的相应输出，在种系样品中检测拷贝数变异。
# 建议使用此工具时将intervals文件切割成几个子集，以免内存溢出。
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

#处理GermlineCNVCaller的输出并生成相应的VCF文件。
#此工具生成“intervals”和“segments”VCF文件，用于补充上一步GermlineCNVCaller的信息。
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