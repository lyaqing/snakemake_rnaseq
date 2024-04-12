# Rule to generate the STAR genome index
rule star_index:
    input:
        fa = config["refs"]["genome"],
        gtf = lambda wc: config["refs"]["ensembl_gtf"] if config["annotation"]["ensembl"] else config["refs"]["refseq_gtf"]
    output:
        index = directory(config["indexes"]["star"])
    message:
        "Generating STAR genome index"
    params:
        sjdb_overhang = config["star_index"]["sjdbOverhang"],
        limit_genome_generate_ram = config["star_index"]["limitGenomeGenerateRAM"],
        genome_sa = config["star_index"]["genomeSAindexNbases"],
        genome_chr_bin_n_bits = config["star_index"]["genomeChrBinNbits"]
    threads: 10
    resources: mem_mb=100000
    shell:
        "mkdir -p {output.index};" # if directory not created STAR will ask for it
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.index} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.sjdb_overhang} --limitGenomeGenerateRAM {params.limit_genome_generate_ram} --genomeSAindexNbases {params.genome_sa} --genomeChrBinNbits {params.genome_chr_bin_n_bits}"


# Rule to map reads using STAR
rule star_alignment:
    input:
        index = config["indexes"]["star"],
        fq1 = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        bam = RESULT_DIR + "star/{sample}.Aligned.sortedByCoord.out.bam",
        bai = RESULT_DIR + "star/{sample}.Aligned.sortedByCoord.out.bai",
        sj = RESULT_DIR + "star/{sample}.SJ.out.tab",
        gene_reads = RESULT_DIR + "star/{sample}.ReadsPerGene.out.tab"
    log:
        log = RESULT_DIR + "star/{sample}.Log.out",
        log_final = RESULT_DIR + "star/{sample}.Log.final.out"
    message:
        "Mapping {wildcards.sample} reads to the genome using STAR."
    threads: 20
    resources:
        cpus=20,
    params:
        prefix = RESULT_DIR + "star/{sample}."
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fq1} {input.fq2} --runThreadN {threads} --readFilesCommand zcat --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --quantMode GeneCounts --twopassMode Basic --outSAMunmapped Within"
        "samtools index {output.bam}"
