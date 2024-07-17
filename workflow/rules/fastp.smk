rule fastp_trim:
    input:
        fq1=lambda wc: samples.loc[wc.sample, "fq1"],
        fq2=lambda wc: samples.loc[wc.sample, "fq2"] #None if pd.isnull(samples.loc[wc.sample, "fq2"]) else samples.loc[wc.sample, "fq2"]
    output:
        trimmed_fq1=WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        trimmed_fq2=WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",
        html=WORKING_DIR + "fastp/{sample}_fastp.html",
        json=WORKING_DIR + "fastp/{sample}_fastp.json"
    message:
        "Trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        phred = config["fastp"]["qualified_quality_phred"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastp --thread {threads} --in1 {input.fq1} --out1 {output.trimmed_fq1} --in2 {input.fq2} --out2 {output.trimmed_fq2} --html {output.html} --json {output.json} --qualified_quality_phred {params.phred}
        """