rule fastqc:
    input:
        fq1 = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        fq1_html = RESULT_DIR + "fastqc/{sample}_R1_trimmed_fastqc.html",
        fq1_zip = RESULT_DIR + "fastqc/{sample}_R1_trimmed_fastqc.zip",
        fq2_html = RESULT_DIR + "fastqc/{sample}_R2_trimmed_fastqc.html",
        fq2_zip = RESULT_DIR + "fastqc/{sample}_R2_trimmed_fastqc.zip"
    params:
        outdir = RESULT_DIR + "fastqc"
    threads: 10
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        fastqc {input.fq1} {input.fq2} -t {threads} -o {params.outdir}
        """
