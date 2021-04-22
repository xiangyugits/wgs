rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        temp("trimmed/{sample}-{unit}.fastq.gz"),
    params:
        **config["params"]["trimmomatic"]["se"],
        extra="",
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.60.7/bio/trimmomatic/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp("trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp("trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp("trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog="trimmed/{sample}-{unit}.trimlog.txt",
    threads:thread
    params:
        **config["params"]["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.60.7/bio/trimmomatic/pe"

#bwa 比对
rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        temp("mapped/{sample}-{unit}.sorted.bam"),
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: thread
    wrapper:
        "0.60.7/bio/bwa/mem"

#去除重复序列
rule mark_duplicates:
    input:
        "mapped/{sample}-{unit}.sorted.bam",
    output:
        bam=temp("dedup/{sample}-{unit}.bam"),
        metrics="qc/dedup/{sample}-{unit}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}-{unit}.log",
    params:
        config["params"]["picard"]["MarkDuplicates"],
    wrapper:
        "0.60.7/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        bam=protected("recal/{sample}-{unit}.bam"),
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log",
    wrapper:
        "0.57.0/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools/index/{prefix}.log",
    wrapper:
        "0.60.7/bio/samtools/index"
