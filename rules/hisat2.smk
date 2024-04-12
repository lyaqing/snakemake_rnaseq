rule hisat2_index:
    input:
        fa = config["refs"]["genome"],
        gtf = lambda wc: config["refs"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["refs"]["refseq_gtf"]
    output:
        ss = config["indexes"]["hisat2"] + "/genome.ss",
        exon = config["indexes"]["hisat2"] + "/genome.exon",
        index = directory(config["indexes"]["hisat2"])
    message:
        "Generating HISAT2 genome index"
    shell:
        "mkdir -p {output.index}; "
        "hisat2_extract_splice_sites.py {input.gtf} > {output.ss} &"
        "hisat2_extract_exons.py {input.gtf} > {output.exon} &"
        "hisat2-build -p 30 --ss {output.ss} --exon {output.exon} {input.fa} {output.index}/genome"


rule hisat2_alignment:
    input:
        index = config["indexes"]["hisat2"],
        fq1 = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        bam = RESULT_DIR + "hisat2/{sample}.bam",
        bai = RESULT_DIR + "hisat2/{sample}.bai",
        novel_splicesite = RESULT_DIR + "hisat2/{sample}.novel_splicesite.txt"
    log:
        RESULT_DIR + "logs/hisat2/{sample}.log"
    params:
        output_dir = RESULT_DIR + "hisat2/"
    message:
        "Mapping {wildcards.sample} reads to the genome using HISAT2."
    threads: 30
    shell:
        "mkdir -p {params.output_dir}; "
        "hisat2 -p {threads} -x {input.index}/genome --rna-strandness RF -1 {input.fq1} -2 {input.fq2} --novel-splicesite-outfile {output.novel_splicesite} | samtools view -Sb | samtools sort -o {output.bam}"
        "samtools index {output.bam}"
