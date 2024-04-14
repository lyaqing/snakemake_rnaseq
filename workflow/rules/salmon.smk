rule salmon_index:
    input:
        transcriptome = config['refs']['transcriptome'],
        genome = config['refs']['genome']
    output:
        index = directory(config['indexes']['salmon'])
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        mkdir -p {output.index};
        grep "^>" {input.genome} | cut -d "" -f 1 > {output.index}/decoys.txt;
        sed -i.bak -e 's/>//g' {output.index}/decoys.txt;
        cut -d "" -f 1 {input.transcriptome} > {output.index}/salmon.cdna.fa;
        cat {output.index}/salmon.cdna.fa {input.genome} > {output.index}/gentrome.fa.gz;
        salmon index -t {output.index}/gentrome.fa.gz -d {output.index}/decoys.txt -i {output.index}
        """


rule salmon_quant:
    input:
        fq1 = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",
        index = config['indexes']['salmon'],
        gtf = lambda wc: config["refs"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["refs"]["refseq_gtf"]
    output:
        quant_genes = RESULT_DIR + "salmon/{sample}/quant.genes.sf",
        quant = RESULT_DIR + "salmon/{sample}/quant.sf"
    log:
        log = RESULT_DIR + "salmon/{sample}/logs/salmon_quant.log"
    threads: 30
    params:
        lib_type = lambda wildcards: '-l ISR' if samples.loc[wildcards.sample, 'ss'] == 'true' else '-l IU',
        output_dir = RESULT_DIR + "salmon/{sample}"
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        mkdir -p {params.output_dir};
        salmon quant -p {threads} -i {input.index} {params.lib_type} -1 {input.fq1} -2 {input.fq2} --validateMappings --gcBias --seqBias -g {input.gtf} -o {params.output_dir}
        """