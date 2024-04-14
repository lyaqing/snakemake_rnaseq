rule multiqc:
  input:
    expand(WORKING_DIR + "fastp/{sample}_fastp.json", sample = SAMPLES)
  output:
    RESULT_DIR + "multiqc_report.html"
  params:
    fastp_directory = WORKING_DIR + "fastp/",
    outdir = RESULT_DIR
  message: "Summarising fastp reports with multiqc"
  shell:
    "multiqc --force "
    "--outdir {params.outdir} "
    "{params.fastp_directory} "