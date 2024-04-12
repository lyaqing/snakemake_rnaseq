########### Libraries ###########
import pandas as pd
import datrie

####### Configuration #######
configfile: "config/config.yaml"
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]

# Load sample data
samples = pd.read_table(config["units"], dtype=str).set_index("sample", drop=False)
SAMPLES=samples.index.unique()

# Include rules
include: "rules/fastp.smk"
include: "rules/star.smk"
include: "rules/subread.smk"
include: "rules/hisat2.smk"
include: "rules/stringtie.smk"
include: "rules/featurecounts.smk"
include: "rules/htseq.smk"
include: "rules/salmon.smk"
include: "rules/sailfish.smk"
include: "rules/kallisto.smk"
include: "rules/rsem.smk"
# include: "rules/multiqc.smk"


##### Define a function to determine which BAM files to use
def get_bam_paths(wildcards):
    paths = []
    if config["alignment"]["star"]:
        paths.append(RESULT_DIR + "star/{sample}.Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample))
    if config["alignment"]["hisat2"]:
        paths.append(RESULT_DIR + "hisat2/{sample}.bam".format(sample=wildcards.sample))
    if config["alignment"]["subread"]:
        paths.append(RESULT_DIR + "subread/{sample}.bam".format(sample=wildcards.sample))
    return paths

# Utility function to convert BAM paths to BAI paths
def bam_to_bai(bam_paths):
    return [path + ".bai" for path in bam_paths]


##### Define final outputs
if config["quantification"]['star']:
  rule all:
    input:
      expand(RESULT_DIR + "star/{sample}.ReadsPerGene.out.tab",sample=SAMPLES)
elif config["quantification"]['stringtie']:
  rule all:
    input:
      expand(RESULT_DIR + "stringtie/{sample}.gene.abundance.txt", sample=SAMPLES),
      expand(RESULT_DIR + "stringtie/{sample}.cov.ref.gtf", sample=SAMPLES),
      expand(RESULT_DIR + "stringtie/{sample}.ballgown.gtf", sample=SAMPLES)
elif config["quantification"]['featurecounts']:
  rule all:
    input:
      expand(RESULT_DIR + "featurecounts/{sample}", sample=SAMPLES),
      expand(RESULT_DIR + "featurecounts/{sample}.jcounts", sample=SAMPLES),
      expand(RESULT_DIR + "featurecounts/{sample}.summary", sample=SAMPLES)
elif config["quantification"]['htseq']:
  rule all:
    input:
      expand(RESULT_DIR + "htseq/{sample}.count", sample=SAMPLES)
elif config["quantification"]['rsem']:
  rule all:
    input:
      expand(RESULT_DIR + "rsem/{sample}.genes.results", sample=SAMPLES),
      expand(RESULT_DIR + "rsem/{sample}.isoforms.results", sample=SAMPLES),
      expand(RESULT_DIR + "rsem/{sample}.stat/{sample}.cnt", sample=SAMPLES),
      expand(RESULT_DIR + "rsem/{sample}.stat/{sample}.model", sample=SAMPLES),
      expand(RESULT_DIR + "rsem/{sample}.stat/{sample}.theta", sample=SAMPLES)
elif config["quantification"]["kallisto"]:
  rule all:
    input:
      expand(RESULT_DIR + "kallisto/{sample}/abundance.tsv", sample=SAMPLES),
      expand(RESULT_DIR + "kallisto/{sample}/abundance.h5", sample=SAMPLES),
      expand(RESULT_DIR + "kallisto/{sample}/run_info.json", sample=SAMPLES),
      expand(RESULT_DIR + "kallisto/{sample}/fusion.txt", sample=SAMPLES)
elif config["quantification"]["sailfish"]:
  rule all:
    input:
      
elif config["quantification"]["salmon"]:
  rule all:
    input:
      expand(RESULT_DIR + "salmon/{sample}/quant.genes.sf", sample=SAMPLES),
      expand(RESULT_DIR + "salmon/{sample}/quant.sf", sample=SAMPLES)
else:
  raise ValueError('please specify only "True" or "False" for the "quantification" parameter in the config file.')

