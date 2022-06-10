import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts('config/DESeq2comparisons.yaml')
print(BASES)
print(CONTRASTS)

SALMON_DIR = get_output_dir(config['project_top_level'], config['salmon_output_folder'])

DESEQ2_DIR = get_output_dir(config['project_top_level'], config['DESeq2_output'])
DESEQ2_DIR = DESEQ2_DIR + "salmon/"


rule deseqOutput:
    input:
        expand(os.path.join(DESEQ2_DIR,"{bse}-{contrast}" + ".DESeq2_normalized_counts.csv"),zip, bse = BASES,contrast = CONTRASTS)



rule run_standard_deseq_salmon:
    input:
        base_group = lambda wildcards: salmon_files_from_contrast(wildcards.bse),
        contrast_group = lambda wildcards: salmon_files_from_contrast(wildcards.contrast),
        tx2gene = 
        
    wildcard_constraints:
        bse="|".join(BASES),
        contrast="|".join(CONTRASTS)
    output:
        os.path.join(DESEQ2_DIR,"{bse}_{contrast}" + "normed_counts.csv.gz")
    params:
        metadata_path = config['sampleCSVpath']
        salmon_output = SALMON_DIR,
        baseName = "{bse}",
        contrastName = "{contrast}",
        out = DESEQ2_DIR + "{bse}_{contrast}",
        base_grep = lambda wildcards: sample_names_from_contrast(wildcards.bse),
        contrast_grep = lambda wildcards: sample_names_from_contrast(wildcards.contrast)
    shell:
        """
        Rscript scripts/salmon_deseq2.R \
        --salmondirectory {params.salmon_output} \
        --metadata {params.metadata_path} \
        --tx2gene '{params.contrast_grep}' \
        --outputdir '{params.bam_suffix}' \
        --column_name '{params.bam_suffix}' \
        --baseline '{params.bam_suffix}' \
        --contrast '{params.bam_suffix}' \
        --controls_name '{params.bam_suffix}' \
        --contrast_name '{params.bam_suffix}' 
        """
