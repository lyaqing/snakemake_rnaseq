rule sailfish_index:
    input:
        fa = config["refs"]["transcriptome"]
    output:
        index = directory(config["indexes"]["sailfish"])
    conda:
        "envs/sailfish.yaml"
    shell:
        "sailfish index -t {input.fa} -o {output.index}"

rule sailfish_quant:
    input:
        index = config["indexes"]["sailfish"],
        fq1 = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        quant = directory(RESULT_DIR + "sailfish")
    params:
        lib_type = lambda wildcards: '-l ISR' if samples.loc[wildcards.sample, 'ss'] == 'true' else '-l IU'
    conda:
        "envs/sailfish.yaml"
    threads: 30
    shell:
        "sailfish quant -p {threads} -i {input.index} {params.lib_type} "
        "-1 <(gzip -dc {input.fq1}) -2 <(gzip -dc {input.fq2}) -o {output.quant}"
