rule featurecounts:
    input:
        bam = lambda wc: get_bam_paths(wc),
        gtf = lambda wc: config["ref"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["ref"]["refeq_gtf"]
    output:
        counts = RESULT_DIR + "featurecounts/{sample}",
        jcounts = RESULT_DIR + "featurecounts/{sample}.jcounts",
        summary = RESULT_DIR + "featurecounts/{sample}.summary"
    params:
        output_dir = RESULT_DIR + "featurecounts",
        prefix = RESULT_DIR + "featurecounts/{sample}",
        ss = lambda wildcards: '-s 2' if samples.loc[wildcards.sample, 'ss'] == 'true' else '-s 0'
    shell:
        "mkdir -p {params.output_dir};"
        "featureCounts -a {input.gtf} -B -C -g gene_id -J -o {params.prefix} -p --countReadPairs {params.ss} -t exon {input.bam}"
