rule qualimap:
    input:
        bam = lambda wc: get_bam_paths(wc),
        bai = lambda wc: bam_to_bai(get_bam_paths(wc)),
        gtf = lambda wc: config["ref"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["ref"]["refeq_gtf"]
    output:
        rnaseq = RESULT_DIR + "qualimap/{sample}_qualimap_rnaseq.tar.gz",
        bamqc = RESULT_DIR + "qualimap/{sample}_qualimap_bamqc.tar.gz"
    log:
        stdout = RESULT_DIR + "logs/qualimap/{sample}.o",
        stderr = RESULT_DIR + "logs/qualimap/{sample}.e"
    params:
        outdir = RESULT_DIR + "qualimap",
        lib_type = lambda wildcards: '-p strand-specific-reverse' if samples.loc[wildcards.sample, 'ss'] == 'true' else '-p non-strand-specific',
    resources:
        mem_gb = 100,
    threads: 10,
    conda:
        "../envs/qc.yaml"
    shell:
        """
        qualimap bamqc -bam {input.bam} \
        -pe {params.lib_type} \
        -outdir {params.outdir}/{wildcards.sample}/bamqc -outformat PDF:HTML \
        --java-mem-size={resources.mem_gb}G -nt {threads}
        tar -zcvf {params.outdir}/{wildcards.sample}_qualimap_bamqc.tar.gz {params.outdir}/{wildcards.sample}/bamqc/*
        
        qualimap rnaseq -bam {input.bam} -gtf {input.gtf} \
        -pe {params.lib_type} \
        -outdir {params.outdir}/{wildcards.sample}/rnaseq -outformat HTML \
        --java-mem-size={resources.mem_gb}G
        tar -zcvf {params.outdir}/{wildcards.sample}_qualimap_rnaseq.tar.gz {params.outdir}/{wildcards.sample}/rnaseq/*
        """