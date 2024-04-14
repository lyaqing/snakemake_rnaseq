# Rule to generate the STAR genome index
rule star_index:
    input:
        fa = config["ref"]["genome"],
        gtf = lambda wc: config["ref"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["ref"]["refeq_gtf"]
    output:
        index = directory(config["index"]["star"])
    message:
        "Generating STAR genome index"
    threads: 10
    resources: mem_mb=100000
    shell:
        """
        mkdir -p {output.index};" # if directory not created STAR will ask for it
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.index} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf}
        """


# Rule to map reads using STAR
rule star_alignment:
    input:
        index = config["index"]["star"],
        fq1 = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        bam = RESULT_DIR + "star/{sample}.Aligned.sortedByCoord.out.bam",
        bai = RESULT_DIR + "star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        sj = RESULT_DIR + "star/{sample}.SJ.out.tab",
        gene_reads = RESULT_DIR + "star/{sample}.ReadsPerGene.out.tab"
    log:
        log = RESULT_DIR + "star/{sample}.Log.out",
        log_final = RESULT_DIR + "star/{sample}.Log.final.out"
    message:
        "Mapping {wildcards.sample} reads to the genome using STAR."
    threads: 20
    resources:
        cpus=20,
    params:
        prefix = RESULT_DIR + "star/{sample}."
    shell:
        """
        STAR --genomeDir {input.index} --readFilesIn {input.fq1} {input.fq2} --runThreadN {threads} --readFilesCommand zcat --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --quantMode GeneCounts --twopassMode Basic --outSAMunmapped Within
        samtools index {output.bam}
        """