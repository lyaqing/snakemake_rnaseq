rule htseq:
    input:
        bam = lambda wc: get_bam_paths(wc),
        bai = lambda wc: bam_to_bai(get_bam_paths(wc)),
        gtf = lambda wc: config["ref"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["ref"]["refeq_gtf"]
    output:
        counts = RESULT_DIR + "htseq/{sample}.count"
    threads: 30
    params:
         output_dir = RESULT_DIR + "htseq"
    shell:
        "mkdir -p {params.output_dir};"
        "htseq-count -n {threads} -s no -r pos {input.bam} {input.gtf} > {output.counts}"
