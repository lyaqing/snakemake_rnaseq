rule kallisto_index:
    input:
        fa = config["refs"]["transcriptome"]
    output:
        index = config['indexes']['kallisto']
    message:
        "Indexing {input.fa} with Kallisto"
    shell:
        "kallisto index -i {output.index} ${fa}"


rule kallisto_quant:
    input:
        index = config['indexes']['kallisto'],
        fq1 = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        tsv = RESULT_DIR + "kallisto/{sample}/abundance.tsv",
        h5 = RESULT_DIR + "kallisto/{sample}/abundance.h5",
        json = RESULT_DIR + "kallisto/{sample}/run_info.json",
        fusion = RESULT_DIR + "kallisto/{sample}/fusion.txt"
    params:
        output_dir = RESULT_DIR + "kallisto/{sample}/",
        ss = lambda wildcards: '--rf-stranded' if samples.loc[wildcards.sample, 'ss'] == 'true' else ''
    message:
        "Quantifying {wildcards.sample} with Kallisto."
    shell:
        "kallisto quant -i {input.index} -o {params.output_dir} {params.ss} --bias {input.fq1} {input.fq2}"
