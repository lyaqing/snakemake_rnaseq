rule subread_index:
    input:
        fa = config["ref"]["genome"]
    output:
        index = directory(config["index"]["subread"])
    conda:
        "../envs/main.yaml"
    shell:
        "subread-buildindex -o {output.index} {input.fa}"


rule subread_alignment:
    input:
        fq1 = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",
        index = config["index"]["subread"],
        gtf = lambda wc: config["ref"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["ref"]["refeq_gtf"]
    output:
        bam = RESULT_DIR + "subread/{sample}.bam",
        bai = RESULT_DIR + "subread/{sample}.bam.bai"
    log:
        RESULT_DIR + "logs/subread/{sample}.log.txt"
    params:
        output_dir = RESULT_DIR + "subread/"
    message:
        "Mapping {wildcards.sample} reads to the genome using Subread."
    threads: 16
    conda:
        "../envs/main.yaml"
    shell:
        "mkdir -p {params.output_dir};"
        "subread-align -t 0 -T {threads} --multiMapping -B 2 -a {input.gtf} -i {input.index} -r {input.fq2} -R {input.fq1} | samtools sort -@ {threads} -o {output.bam}"
        "samtools index {output.bam}"
