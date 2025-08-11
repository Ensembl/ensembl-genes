def build_annotation_commands(core_adaptor: dict, output_params: dict, anno_settings: dict, settings: dict) -> None:
    """
    Construct and store command-line strings for genome annotation and repeat analysis.

    This function builds two command-line strings based on the provided configuration:
    - `anno_commandline`: full annotation pipeline
    - `anno_red_commandline`: repeat and masking pipeline

    The resulting command strings are stored in the `output_params` dictionary under
    the keys "anno_commandline" and "anno_red_commandline".

    Args:
        core_adaptor (dict): Contains database connection info with keys:
            'dbname', 'host', 'port', 'user', 'pass'.
        output_params (dict): Dictionary where input/output file paths and
            directories are specified. The final command strings are added here.
        anno_settings (dict): Annotation settings such as number of threads and
            diamond validation DB path.
        settings (dict): General pipeline settings, including repeat library usage.

    Returns:
        None: Modifies `output_params` in place by adding two command strings.
    """

    get = lambda k: output_params.get(k, "")  # Short helper

    anno_commandline = (
        f" --genome_file {get('reheadered_toplevel_genome_file')}"
        f" --db_details {core_adaptor['dbname']},{core_adaptor['host']},{core_adaptor['port']},{core_adaptor['user']},{core_adaptor['pass']}"
        f" --output_dir {get('output_path')}"
        f" --short_read_fastq_dir {get('short_read_dir')}"
        f" --long_read_fastq_dir {get('long_read_dir')}"
        f" --max_intron_length {get('max_intron_length')}"
        f" --protein_file {get('protein_file')}"
        f" --busco_protein_file {get('busco_protein_file')}"
        f" --rfam_accessions_file {get('rfam_accessions_file')}"
        f" --num_threads {anno_settings['num_threads']}"
    )

    if settings.get("use_existing_repeatmodeler_library"):
        anno_commandline += (
            f" --repeatmasker_library {settings['use_existing_repeatmodeler_library']}"
        )

    if anno_settings.get("diamond_validation_db"):
        anno_commandline += (
            f" --diamond_validation_db {anno_settings['diamond_validation_db']}"
        )

    if output_params.get("validation_type"):
        anno_commandline += f" --validation_type {output_params['validation_type']}"

    anno_commandline += " --run_full_annotation --load_to_ensembl_db"

    anno_red_commandline = (
        f" --genome_file {get('reheadered_toplevel_genome_file')}"
        f" --db_details {core_adaptor['dbname']},{core_adaptor['host']},{core_adaptor['port']},{core_adaptor['user']},{core_adaptor['pass']}"
        f" --output_dir {get('output_path')}"
        f" --num_threads {get('num_threads')}"
        " --run_masking --run_repeats --run_simple_features --load_to_ensembl_db"
    )

    # Add both commands back into the params dict
    output_params["anno_commandline"] = anno_commandline
    output_params["anno_red_commandline"] = anno_red_commandline
