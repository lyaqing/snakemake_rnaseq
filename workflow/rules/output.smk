##################################
# Produce table of raw gene counts
##################################

if config["count_multimappers"]:
  rule create_counts_table:
    input:
      bams = expand(RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES),
      gff  = config["ref"]["gtf"]
    output:
      WORKING_DIR + "raw_counts.tsv"
    message: "Producing the table of raw counts (counting read multimappers)"
    threads: 10
    shell:
      "featureCounts -T {threads} -M -t exon -g transcript_id -F 'gtf' -a {input.gff} -o {output} {input.bams}"
else:
  rule create_counts_table:
    input:
      bams = expand(RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES),
      gff  = config["ref"]["gtf"]
    output:
      WORKING_DIR + "raw_counts.tsv"
    message: "Producing the table of raw counts (not counting multimappers)"
    threads: 10
    shell:
      "featureCounts -T {threads} -t exon -g transcript_id -F 'gtf' -a {input.gff} -o {output} {input.bams}"


rule parse_raw_counts:
  input:
    WORKING_DIR + "raw_counts.tsv"
  output:
    RESULT_DIR + "raw_counts.parsed.tsv"
  message: 
    "Parsing the raw counts file for scaling (removal of comment and header renaming)"
  params:
    star_result_dir_name = RESULT_DIR + "star/"
  shell:
    "python scripts/parse_raw_counts.py {input} {params.star_result_dir_name} {output}"

#########################################
# Produce table of normalised gene counts
#########################################

rule normalise_raw_counts:
  input:
    raw = RESULT_DIR + "raw_counts.parsed.tsv"
  output:
    norm = RESULT_DIR + "scaled_counts.tsv"
  message: 
    "Normalising raw counts the DESeq2 way"
  shell:
    "Rscript --vanilla scripts/deseq2_normalization.R {input.raw} {output.norm}"

#############################
# Plots of mapping statistics
#############################
