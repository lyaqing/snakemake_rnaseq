rule stringtie:
    input:
        bam = lambda wc: get_bam_paths(wc),
        bai = lambda wc: bam_to_bai(get_bam_paths(wc)),
        gtf = lambda wc: config["ref"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["ref"]["refeq_gtf"]
    output:
        gene_abundance = RESULT_DIR + "stringtie/{sample}.gene.abundance.txt",
        cov_ref_gtf = RESULT_DIR + "stringtie/{sample}.cov.ref.gtf",
        ballgown_gtf = RESULT_DIR + "stringtie/{sample}.ballgown.gtf"
    params:
        output_dir = RESULT_DIR + "stringtie/",
        ss = lambda wildcards: '--rf' if samples.loc[wildcards.sample, 'ss'] == 'true' else ''
    threads: 60
    conda:
        "../envs/main.yaml"
    shell:
        "mkdir -p {params.output_dir};"
        "stringtie {input.bam} {params.ss} -e -p {threads} -G {input.gtf} -o {output.ballgown_gtf} -C {output.cov_ref_gtf} -A {output.gene_abundance}"