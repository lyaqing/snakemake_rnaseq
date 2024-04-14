rule rsem_prepare_reference:
    input:
        gtf = lambda wc: config["ref"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["ref"]["refeq_gtf"],
        fa = config["ref"]["genome"]
    output:
        index = directory(config["index"]["rsem"])
    message:
        "Generating RSEM genome index"
    conda:
        "../envs/rsem.yaml"
    shell:
        "mkdir -p {output.index};"
        "rsem-prepare-reference --gtf {input.gtf} --bowtie {input.fa} {output.index}/rsem"


rule rsem_calculate_expression:
    input:
        index = config["index"]["rsem"],
        fq1 = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        genes = RESULT_DIR + "rsem/{sample}.genes.results",
        isoforms = RESULT_DIR + "rsem/{sample}.isoforms.results",
        stat_cnt = RESULT_DIR + "rsem/{sample}.stat/{sample}.cnt",
        stat_model = RESULT_DIR + "rsem/{sample}.stat/{sample}.model",
        stat_theta = RESULT_DIR + "rsem/{sample}.stat/{sample}.theta"
    message:
        "Quantifying {wildcards.sample} using RSEM"
    conda:
        "../envs/rsem.yaml"
    params:
        ss = lambda wildcards: '--strandedness reverse' if samples.loc[wildcards.sample, 'ss'] == 'true' else '--strandedness none',
        output_dir = RESULT_DIR + "rsem",
        prefix = RESULT_DIR + "rsem/{sample}"
    shell:
        "mkdir -p {params.output_dir};"
        "rsem-calculate-expression --no-bam-output --paired-end {input.fq1} {input.fq2} {input.index}/rsem {params.prefix}"
