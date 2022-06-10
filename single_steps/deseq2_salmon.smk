

rule majiq_delta_psi:
    input:
        base_group_deseq2 = lambda wildcards: salmon_files_from_metadata(wildcards.bse),
        contrast_group_deseq2 = lambda wildcards: salmon_files_from_metadata(wildcards.contrast)
    output:
        os.path.join(MAJIQ_DIR,"delta_psi","{bse}-{contrast}" + "deseq2_results.csv")
    params:
        majiq_path = config['majiq_path'],
        delta_psi_output_folder = os.path.join(MAJIQ_DIR,"delta_psi"),
        majiq_psi_extra_parameters = return_parsed_extra_params(config['majiq_psi_extra_parameters'])
    threads:
        8
    shell:
        """
        mkdir -p {params.delta_psi_output_folder}
        salmon_deseq2.R \
        --basefiles {input.base_group_deseq2} \
        --contrastfiles {input.contrast_group_deseq2} 
        -o {params.delta_psi_output_folder} \
        --basename {wildcards.bse}\
        --contrastname {wildcards.contrast}
        """