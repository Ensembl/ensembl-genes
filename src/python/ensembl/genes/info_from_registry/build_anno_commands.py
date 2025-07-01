def build_annotation_commands(adaptors, output_params):
    core_db = adaptors["core_string"]
    get = lambda k: output_params.get(k, "")  # Short helper

    anno_commandline = (
        f" --genome_file {get('reheadered_toplevel_genome_file')}"
        f" --db_details {core_db['dbname']},{core_db['host']},{core_db['port']},{core_db['user']},{core_db['pass']}"
        f" --output_dir {get('output_path')}"
        f" --short_read_fastq_dir {get('short_read_dir')}"
        f" --long_read_fastq_dir {get('long_read_dir')}"
        f" --max_intron_length {get('max_intron_length')}"
        f" --protein_file {get('protein_file')}"
        f" --busco_protein_file {get('busco_protein_file')}"
        f" --rfam_accessions_file {get('rfam_accessions_file')}"
        f" --num_threads {get('num_threads')}"
    )

    if output_params.get("repeatmodeler_library"):
        anno_commandline += (
            f" --repeatmasker_library {output_params['repeatmodeler_library']}"
        )

    if output_params.get("diamond_validation_db"):
        anno_commandline += (
            f" --diamond_validation_db {output_params['diamond_validation_db']}"
        )

    if output_params.get("validation_type"):
        anno_commandline += f" --validation_type {output_params['validation_type']}"

    anno_commandline += " --run_full_annotation --load_to_ensembl_db"

    anno_red_commandline = (
        f" --genome_file {get('reheadered_toplevel_genome_file')}"
        f" --db_details {core_db['dbname']},{core_db['host']},{core_db['port']},{core_db['user']},{core_db['pass']}"
        f" --output_dir {get('output_path')}"
        f" --num_threads {get('num_threads')}"
        " --run_masking --run_repeats --run_simple_features --load_to_ensembl_db"
    )

    # Add both commands back into the params dict
    output_params["anno_commandline"] = anno_commandline
    output_params["anno_red_commandline"] = anno_red_commandline
