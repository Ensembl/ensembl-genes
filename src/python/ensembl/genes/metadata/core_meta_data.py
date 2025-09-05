import requests
import xmltodict
import json
from datetime import datetime
import argparse
import pymysql
import re
import logging
import logging.config
from pathlib import Path


def mysql_fetch_data(query, database, host, port, user):
    try:
        conn = pymysql.connect(host=host, user=user, port=port, database=database.strip())

        cursor = conn.cursor()
        cursor.execute(query)
        info = cursor.fetchall()

    except pymysql.Error as err:
        print(err)

    cursor.close()
    conn.close()
    return info


def get_ena_metadata(accession, truth_dict):
    assembly_url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
    assembly_xml = requests.get(assembly_url)
    assembly_dict = xmltodict.parse(assembly_xml.text)

    assembly_attribs = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"]["ASSEMBLY_ATTRIBUTES"]["ASSEMBLY_ATTRIBUTE"]

    return_dict = {}

    # assembly meta data
    # assembly name
    try:
        return_dict["assembly.name"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"]["NAME"].replace(" ", "_")
    except KeyError:
        return_dict["assembly.name"] = ""
    # assembly level
    try:
        return_dict["assembly.level"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"]["ASSEMBLY_LEVEL"]
    except KeyError:
        return_dict["assembly.level"] = ""

    for attrib_set in assembly_attribs:
        # assembly date
        if attrib_set["TAG"] == "ENA-LAST-UPDATED":
            return_dict["assembly.date"] = attrib_set["VALUE"]

        # organism meta data
        # sample id
        if "organism.biosample_id" not in truth_dict:
            try:
                return_dict["organism.biosample_id"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"]["SAMPLE_REF"]["IDENTIFIERS"]["PRIMARY_ID"]
            except KeyError:
                logger.critical(
                    " | BIOSAMPLE_ID | organism.biosample_id could not be found in the ENA metadata, this is a required key!"
                )
        # taxonomy id
        # found a species that has incorrect taxon id in INSDC records, hardcoding the fix
        if accession == "GCA_944452655.1" or accession == "GCA_944452715.1":
            return_dict["organism.taxonomy_id"] = "1539398"
        else:
            try:
                return_dict["organism.taxonomy_id"] = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"]["TAXON"]["TAXON_ID"]
            except KeyError:
                logger.warning(" | TAXONOMY_ID | organism.taxonomy_id could not be found in the ENA metadata")

        return return_dict


def get_ncbi_metadata(accession, assembly_name, scientific_name, search):
    organism = f"{accession}_{assembly_name}"
    ncbi_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession[0:3]}/{accession[4:7]}/{accession[7:10]}/{accession[10:13]}/{organism}/{organism}_assembly_report.txt"
    ncbi_return = requests.get(ncbi_url).text.split("\n")

    return_dict = {}
    return_dict["assembly.ucsc_alias"] = ""
    return_dict["organism.strain"] = ""
    return_dict["organism.strain_type"] = ""

    if search == "ucsc":
        for line in ncbi_return:
            ucscREGEX = re.search("\# Synonyms\: +([a-zA-Z0-9]+)", line)
            if ucscREGEX:
                ucsc_alias = ucscREGEX.group(1)
                return_dict["assembly.ucsc_alias"] = ucsc_alias

    elif search == "biosample":
        for line in ncbi_return:
            strainREGEX = re.search("\# Infraspecific name: +([A-Za-z]+)=([a-zA-Z0-9 -\.\/]+)", line)
            if strainREGEX:
                strain_type = strainREGEX.group(1)
                strain = strainREGEX.group(2)
                return_dict["organism.strain"] = strain
                return_dict["organism.strain_type"] = strain_type

    return return_dict


def get_biosample_metadata(biosample_id, assembly_accession, assembly_name, scientific_name):
    biosample_url = f"https://www.ebi.ac.uk/biosamples/samples/{biosample_id}"
    biosample_return = requests.get(biosample_url).text
    return_dict = {}
    return_dict["organism.strain"] = ""
    return_dict["organism.strain_type"] = ""

    try:
        biosample_data = json.loads(biosample_return)

        strain_types = {"population", "race", "ecotype", "breed", "strain", "cultivar"}

        for strain_type in strain_types:
            try:
                return_dict["organism.strain"] = biosample_data["characteristics"][strain_type][0]["text"]
                return_dict["organism.strain_type"] = strain_type
                break
            except KeyError:
                continue

        if return_dict["organism.strain"] == "Caucasian":
            return_dict["organism.strain"] = "European"

        try:
            return_dict["assembly.tol_id"] = biosample_data["characteristics"]["tolid"][0]["text"]
        except KeyError:
            return_dict["assembly.tol_id"] = ""

    except json.decoder.JSONDecodeError:
        # I couldn't get the biosample info from ENA because metadata is so ridiculously poorly recorded so I'm going to try and grab it from NCBI... FML
        biosample_dict = get_ncbi_metadata(assembly_accession, assembly_name, scientific_name, "biosample")
        return_dict["organism.strain"] = biosample_dict["organism.strain"]
        return_dict["organism.strain_type"] = biosample_dict["organism.strain_type"]
        return_dict["assembly.tol_id"] = ""

    return return_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare SQL updates for core dbs")
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help="Path where the output and temp files will write to. \
        Uses current dir by default",
    )
    parser.add_argument(
        "-d",
        "--db_name",
        help="Database name",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--host",
        help="Host server",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--port",
        help="Host server port",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--production_name",
        help="species.production_name",
        required=False,
    )
    parser.add_argument(
        "-t", "--team", required=True, type=lambda x: x.capitalize(), help="Team responsible for the database"
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose output (check that all required keys are not NULL/ empty)"
    )
    
    args = parser.parse_args()

    server_info = {
        "staging": {
            "db_host": args.host,
            "db_port": int(args.port),
            "db_user": "ensro",
        },
        "meta": {
            "db_host": "mysql-ens-meta-prod-1.ebi.ac.uk",
            "db_port": 4483,
            "db_user": "ensro",
            "db_pass": "",
        },
    }

    metadata_dir = Path(__file__).parent
    provider_static_file = metadata_dir / "provider_static.txt"
    ref_static_file = metadata_dir / "ref_static.txt"
    snp_static_file = metadata_dir / "snp_static.txt"
    url_static_file = metadata_dir / "url_static.txt"

    db = args.db_name
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if (args.production_name):
        sql_out = open(output_dir / f"{args.production_name}.sql", "w")
    else:
        sql_out = open(output_dir / f"{db}.sql", "w")
        
    print(f"Working on database: {db}")
    print(f"USE {db};", file=sql_out)

    # set up logger
    if (args.production_name):
        log_file_path = output_dir / f"{args.production_name}_metadata.log"
    else:
        log_file_path = output_dir / f"{db}_metadata.log"
    log_ini_path = metadata_dir / "logging.conf"
    logging.config.fileConfig(
        log_ini_path,
        defaults={"logfilename": log_file_path},
        disable_existing_loggers=True,
    )
    logger = logging.getLogger()
    logger.propagate = False

    core_dict = {}
    if (args.production_name):
        core_dict['species.production_name'] = args.production_name
        #getting species_id via the production_name
        species_query = f"SELECT species_id FROM meta WHERE meta_key='species.production_name' AND meta_value='{core_dict['species.production_name']}';"
        species_meta = mysql_fetch_data(
            species_query,
            host=server_info["staging"]["db_host"],
            user=server_info["staging"]["db_user"],
            port=server_info["staging"]["db_port"],
            database=db,
        )
        species_id = species_meta[0][0]
        print("species ID: "+str(species_id))
    else:
        species_id = 1
        print("WARNING: no production_name provided, using default species ID ="+str(species_id))
        
    # get all existing assembly, species and genebuild metadata from the core db
    core_query = f"SELECT meta_key,meta_value FROM meta WHERE species_id = {species_id} AND (meta_key LIKE 'assembly%' OR meta_key LIKE 'species%' OR meta_key LIKE 'genebuild%' OR meta_key LIKE 'organism%' OR meta_key LIKE 'sample%' OR meta_key LIKE 'annotation%' OR meta_key LIKE 'gencode%');"
    print(core_query)
    core_meta = mysql_fetch_data(
        core_query,
        host=server_info["staging"]["db_host"],
        user=server_info["staging"]["db_user"],
        port=server_info["staging"]["db_port"],
        database=db,
    )
    for meta_pair in core_meta:
        core_dict[meta_pair[0]] = meta_pair[1]

        
    # get the assembly metadata from the sources of truth (sources of truth in parentheses)
    # expected assembly meta_keys: assembly.accession (from core), assembly.date (from ena), assembly.is_reference (static), assembly.name (ena), assembly.provider_name (core or default), assembly.provider_url (core or default), assembly.level (ena), assembly.tolid (biosample), assembly.ucsc_alias (ncbi), assembly.long_name, assembly.url_name (static)
    # expected organism meta_keys: organism.taxonomy_id (ena), organism.species_taxonomy_id (taxonomy db), organism.common_name (taxonomy db), organism.strain (biosample), organism.scientific_name (taxonomy db), organism.scientific_parlance_name (static), organism.strain_type (biosample), organism.sample_accession (ena)
    # expected genebuild meta_keys: genebuild.initial_release_date, genebuild.last_geneset_update, genebuild.level, genebuild.method, genebuild.method_display, genebuild.start_date, genebuild.version (create and check and required), genebuild.sample_gene (core), genebuild.sample_location (core), genebuild.id, genebuild.projection_source_db, genebuild.havana_datafreeze_date, genebuild.provider_name (static or core or default), genebuild.provider_url (static or core or default), genebuild.annotation_source (core or default)
    truth_dict = {}

    # now some assembly.accession values will be GCFs - that breaks finding things based on a GCA
    gca_accession = core_dict["assembly.accession"]
    if core_dict["assembly.accession"].startswith("GCF"):
        gca_accession = core_dict["assembly.alt_accession"]
        truth_dict["assembly.alt_accession"] = core_dict["assembly.alt_accession"]
        truth_dict["assembly.accession_refseq"] = core_dict["assembly.accession"]
        truth_dict["assembly.accession_body"] = "RefSeq"
    elif "assembly.accession_refseq" in core_dict:
        # if the core comes from the main site, we need to put the GCF back as the main accession
        truth_dict["assembly.alt_accession"] = core_dict["assembly.alt_accession"]
        core_dict["assembly.accession"] = core_dict["assembly.accession_refseq"]
        truth_dict["assembly.accession_body"] =	"RefSeq"

    # Some goddamn ENA assembly records do not have the BioSample ID, so I'm hardcoding them, I swear to jaysus, I'm so done with metadata!!!
    if db == "caenorhabditis_elegans_core_57_110_282":
        truth_dict["organism.biosample_id"] = "SAMN04256190"
    elif db == "ciona_intestinalis_core_110_3":
        truth_dict["organism.biosample_id"] = "SAMD00414333"
    elif db == "homo_sapiens_37_core_110_37":
        truth_dict["organism.biosample_id"] = "SAMN12121739"
    elif db == "homo_sapiens_core_110_38":
        truth_dict["organism.biosample_id"] = "SAMN12121739"
    elif db == "saccharomyces_cerevisiae_core_57_110_4":
        truth_dict["organism.biosample_id"] = "SAMEA3184125"
    elif db == "mus_musculus_core_110_39":
        truth_dict["organism.biosample_id"] = "SAMN26853311"
    elif db == "bos_taurus_core_110_1":
        truth_dict["organism.biosample_id"] = "SAMN03145444"

    # get metadata from ENA records
    try:
        truth_dict.update(get_ena_metadata(gca_accession, truth_dict))
    except KeyError:
        logger.critical("No assembly accession found, cannot process any further")

    # get common and scientific names from NCBI taxonomy (seems inefficient to query the db twice, but I don't think it returns ordered results and I don't want to risk mixing them up)
    try:
        name_query = f"SELECT name FROM ncbi_taxa_name WHERE taxon_id={truth_dict['organism.taxonomy_id']} AND name_class='genbank common name';"
        name_info = mysql_fetch_data(
            name_query,
            host=server_info["meta"]["db_host"],
            user=server_info["meta"]["db_user"],
            port=server_info["meta"]["db_port"],
            database="ncbi_taxonomy",
        )
        truth_dict["organism.common_name"] = (name_info[0][0]).capitalize()
    except IndexError:  # not everything has a genbank common name
        truth_dict["organism.common_name"] = ""
    s_name_query = f"SELECT name FROM ncbi_taxa_name WHERE taxon_id={truth_dict['organism.taxonomy_id']} AND name_class='scientific name';"
    s_name_info = mysql_fetch_data(
        s_name_query,
        host=server_info["meta"]["db_host"],
        user=server_info["meta"]["db_user"],
        port=server_info["meta"]["db_port"],
        database="ncbi_taxonomy",
    )
    truth_dict["organism.scientific_name"] = (s_name_info[0][0]).capitalize()

    # get metadata from NCBI taxonomy
    taxonomy_query = f"SELECT rank, parent_id FROM ncbi_taxa_node WHERE taxon_id={truth_dict['organism.taxonomy_id']};"
    rank_info = mysql_fetch_data(
        taxonomy_query,
        host=server_info["meta"]["db_host"],
        user=server_info["meta"]["db_user"],
        port=server_info["meta"]["db_port"],
        database="ncbi_taxonomy",
    )
    rank_dict = {}
    for pair in rank_info:
        rank_dict["rank"] = pair[0]
        rank_dict["parent_id"] = str(pair[1])

    if rank_dict["rank"] == "species":
        truth_dict["organism.species_taxonomy_id"] = truth_dict["organism.taxonomy_id"]
    else:
        while rank_dict["rank"] != "species":
            current_id = str(rank_dict["parent_id"])
            parent_taxonomy_query = f"SELECT rank, parent_id FROM ncbi_taxa_node WHERE taxon_id={rank_dict['parent_id']};"
            parent_rank_info = mysql_fetch_data(
                parent_taxonomy_query,
                host=server_info["meta"]["db_host"],
                user=server_info["meta"]["db_user"],
                port=server_info["meta"]["db_port"],
                database="ncbi_taxonomy",
            )
            for pair in parent_rank_info:
                rank_dict["rank"] = pair[0]
                rank_dict["parent_id"] = str(pair[1])

            truth_dict["organism.species_taxonomy_id"] = current_id

    # get metadata from NCBI records
    truth_dict.update(
        get_ncbi_metadata(
            core_dict["assembly.accession"],
            truth_dict["assembly.name"],
            truth_dict["organism.scientific_name"],
            "ucsc",
        )
    )

    # get metadata from BioSample records
    if truth_dict["organism.biosample_id"] != "":
        truth_dict.update(
            get_biosample_metadata(
                truth_dict["organism.biosample_id"],
                gca_accession,
                truth_dict["assembly.name"],
                truth_dict["organism.scientific_name"],
            )
        )
    else:
        # there's probably a better source for ToLIDs - DToL portal?
        truth_dict["assembly.tol_id"] = ""

    # now to get some metadata that can only come from static files
    # scientific_parlance_name
    snp_list = open(snp_static_file).readlines()
    for line in snp_list:
        if truth_dict["organism.scientific_name"] in line:
            snp = line.split("\t")[1].strip()
            truth_dict["organism.scientific_parlance_name"] = snp

    # assembly.url_name
    url_list = open(url_static_file).readlines()
    for line in url_list:
        if gca_accession in line:
            url = line.split("\t")[1].strip()
            truth_dict["assembly.url_name"] = url

    # assembly.is_reference
    ref_list = open(ref_static_file).readlines()
    for line in ref_list:
        if truth_dict["organism.scientific_name"] in line:
            ref_accession = line.split("\t")[1].strip()
            if gca_accession == ref_accession:
                truth_dict["assembly.is_reference"] = 1

    # now to create some values
    # assembly provider and url - if not in core already, set to default provider,"ENA"
    if "assembly.provider_name" not in core_dict:
        truth_dict["assembly.provider_name"] = "ENA"
        truth_dict["assembly.provider_url"] = "https://www.ebi.ac.uk/ena/browser/home"
        logger.warning(" | ASSEMBLY_PROVIDER | No assembly provider information found, using default, ENA.")

    provider_dict = {}
    with open(provider_static_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 2: 
                provider_dict[parts[0]] = {
                    'name': parts[1],
                    'url': parts[2]
                }

    # genebuild.provider_name and _url keys: check if annotation.provider_name and _url keys exist, if not check spreadsheet, else set to default "Ensembl", "Ensembl url" (full_genebuild, anno, braker, hprc)
    # genebuild.version key is now being used for metadata loading, we are setting it to the first genebuild.version for everything (unless there is already a value set for this key), data teams need to be aware when handing over an updated annotation, i.e. assembly is same as an existing genome, but the gene set has been updated (new data, or a fix)
    # maybe add a check on the metadata database here - does the label GCA_XXXX_ENSXX match an existing dataset, then warn! (could be part of the check for required keys!)
    if "gencode.version" in core_dict:
        truth_dict["genebuild.version"] = core_dict["gencode.version"].replace(" ", "")
        truth_dict["genebuild.method_display"] = "Manual annotation"
    elif "genebuild.method" in core_dict:
        # a quick check where sample genes have BRAKER stable id prefixes, because some braker annotations had incorrect value for genebuild.method (I'm doing this with the sample gene text as it will contain the BRAKER stable id prefix)
        if "BRAKER" in core_dict["sample.gene_text"]:
            truth_dict["genebuild.version"] = "BRK01"
            truth_dict["genebuild.method"] = "braker"
            truth_dict["genebuild.method_display"] = "BRAKER2"
            truth_dict["genebuild.annotation_source"] = "braker"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            truth_dict["genebuild.provider_url"] = "https://beta.ensembl.org/help/articles/braker-2-genome-annotation"
        # otherwise check the genebuild.method for annotation key updates/additions
        if "HELIXER" in core_dict["sample.gene_text"]:
            truth_dict["genebuild.version"] = "HLX01"
            truth_dict["genebuild.method"] = "helixer"
            truth_dict["genebuild.method_display"] = "Helixer"
            truth_dict["genebuild.annotation_source"] = "helixer"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            truth_dict["genebuild.provider_url"] = "https://beta.ensembl.org/info/genome/helixer-genome-annotation"
        elif core_dict["genebuild.method"] == "full_genebuild":
            truth_dict["genebuild.version"] = "ENS01"
            truth_dict["genebuild.method_display"] = "Ensembl Genebuild"
            truth_dict["genebuild.annotation_source"] = "ensembl"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            truth_dict["genebuild.provider_url"] = "https://beta.ensembl.org/help/articles/vertebrate-genome-annotation"
        elif core_dict["genebuild.method"] == "anno":
            truth_dict["genebuild.version"] = "ENS01"
            truth_dict["genebuild.method_display"] = "Ensembl Genebuild"
            truth_dict["genebuild.annotation_source"] = "ensembl"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            truth_dict["genebuild.provider_url"] = "https://beta.ensembl.org/help/articles/non-vertebrate-genome-annotation"
        elif core_dict["genebuild.method"] == "projection_build":
            truth_dict["genebuild.version"] = "ENS01"
            truth_dict["genebuild.annotation_source"] = "ensembl"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            if core_dict["species.scientific_name"] == "homo sapiens":
                truth_dict["genebuild.provider_url"] = "https://beta.ensembl.org/help/articles/human-genome-automated-annotation"
                truth_dict["genebuild.method_display"] = "Mapping from GRCh38"
            else:
                truth_dict["genebuild.provider_url"] = "https://beta.ensembl.org/help/articles/vertebrate-genome-annotation"
                truth_dict["genebuild.method_display"] = "Mapping from reference"
        # for non-genebuilds, I don't have a set of rules (as above), so I rely on the core meta keys and the static
        elif core_dict["genebuild.method"] == "import":
            if "genebuild.method_display" not in core_dict:
                truth_dict["genebuild.method_display"] = "Import"
            if "genebuild.version" not in core_dict:
                truth_dict["genebuild.version"] = "EXT01"
            # try to get the annotation source
            if "species.annotation_source" in core_dict and "genebuild.annotation_source" not in core_dict:
                truth_dict["genebuild.annotation_source"] = core_dict["species.annotation_source"]
            elif "genebuild.annotation_source" not in core_dict:
                logger.critical(" | GENEBUILD_PROVIDER_URL | No genebuild.annotation_source could be found in the core db, this is a required key!")
            # try to get the annotation provider name from core
            if "annotation.provider_name" in core_dict and "genebuild.provider_name" not in core_dict:
                truth_dict["genebuild.provider_name"] = core_dict["annotation.provider_name"]
            #otherwise get it from provider_static.txt
            elif core_dict["species.production_name"] in provider_dict and "genebuild.provider_name" not in core_dict:
                truth_dict["genebuild.provider_name"] = provider_dict[core_dict["species.production_name"]["name"]]
            elif "genebuild.provider_name" not in core_dict:
                logger.critical(" | GENEBUILD_PROVIDER_NAME | No genebuild.provider_name could be found either in the core db or in provider_static.txt, this is a required key!")    
            # try to get the annotation provider url from core                    
            if "annotation.provider_url" in core_dict and "genebuild.provider_url" not in core_dict:
                truth_dict["genebuild.provider_url"] = core_dict["annotation.provider_url"]
            #otherwise get it from provider_static.txt
            elif core_dict["species.production_name"] in provider_dict and "genebuild.provider_url" not in core_dict:
                truth_dict["genebuild.provider_url"] = provider_dict[core_dict["species.production_name"]["url"]]
            elif "genebuild.provider_url" not in core_dict:
                logger.critical(" | GENEBUILD_PROVIDER_URL | No genebuild.provider_url could be found either in the core db or in provider_static.txt, this is a required key!")
            
                 
    else:
        logger.warning(" | GENEBUILD_METHOD | No genebuild.method, assuming this is an Ensembl annotation")
        truth_dict["genebuild.version"] = "ENS01"
        truth_dict["genebuild.annotation_source"] = "ensembl"
        truth_dict["genebuild.provider_name"] = "Ensembl"
        truth_dict["genebuild.provider_url"] = "https://rapid.ensembl.org/info/genome/genebuild/full_genebuild.html"

    # if the genebuild.version already exists in the core db, I'll just leave that value
    if "genebuild.version" in core_dict and core_dict["genebuild.version"] != "":
        truth_dict["genebuild.version"] = core_dict["genebuild.version"]

    #if the sample gene info is in the core db, update the meta key names
    try:
        truth_dict["genebuild.sample_gene"] = core_dict["sample.gene_param"]
    except KeyError:
        logger.critical(" | SAMPLE.GENE_PARAM | No sample.gene_param could be found in the core db, this is a required key!")
    try:
        truth_dict["genebuild.sample_location"] = core_dict["sample.location_param"]
    except KeyError:
        logger.critical(" | SAMPLE.LOCATION_PARAM | No sample.location_param could be found in the core db, this is a required key!")
        
    #let's do a check for the genebuild.last_geneset_update key because it's required by web
    try:
        truth_dict["genebuild.last_geneset_update"] = core_dict["genebuild.last_geneset_update"]
    except KeyError:
        logger.warning(" | GENEBUILD.LAST_GENESET_UPDATE | No genebuild.last_geneset_update could be found in the core db, this is a required key, I'm setting it to today.")
        truth_dict["genebuild.last_geneset_update"] = datetime.now().strftime("%Y-%m")

        
    #update the remaining species keys -> organisms keys
    try:
        truth_dict["organism.production_name"] = core_dict["species.production_name"]
    except KeyError:
        logger.critical(" | SPECIES.PRODUCTION_NAME | No species.production_name could be found in the core db, this is a required key!")

    # Set the team responsible for this genome
    truth_dict["genebuild.team_responsible"] = args.team
        
    # if the expected meta_key does not exist in the core metadata but there is a truth value, INSERT
    # if the expected meta_key exists in the core metadata but the meta_value does not match the truth value and truth value is not NULL, UPDATE
    # if the expected meta_key exists in the core metadata and matches the truth value, NO ACTION
    for meta_key in truth_dict:
        # Escape single quotes in truth_dict[meta_key]
        meta_value = truth_dict[meta_key].replace("'", "''") if truth_dict[meta_key] else None
        
        if meta_key not in core_dict:
            if truth_dict[meta_key]:
                print(
                    f"INSERT IGNORE INTO meta (species_id, meta_key, meta_value) VALUES({species_id}, '{meta_key}', '{meta_value}');",
                    file=sql_out,
                )
        elif truth_dict[meta_key] != core_dict[meta_key] and truth_dict[meta_key] != "":
            print(
                f"UPDATE meta SET meta_value='{meta_value}' WHERE meta_key='{meta_key}';",
                file=sql_out,
            )

    # do a check on required keys
    required_meta_keys = [
        "organism.production_name",
        "organism.taxonomy_id",
        "organism.scientific_name",
        "assembly.accession",
        "assembly.name",
        "genebuild.version",
        "genebuild.method",
        "genebuild.method_display",
        "genebuild.start_date",
        "genebuild.annotation_source",
        "genebuild.provider_name",
        "genebuild.provider_url",
        "genebuild.sample_gene",
        "genebuild.sample_location",
        "genebuild.last_geneset_update",
        "organism.biosample_id",
    ]
    for required_key in required_meta_keys:
        if required_key not in core_dict:
            if required_key not in truth_dict:
                logger.critical(f"You are missing required meta key: {required_key}")

    # report if there has been a change in common name
    if "species.common_name" in core_dict:
        if core_dict["species.common_name"].lower() != truth_dict["organism.common_name"].lower():
            logger.warning(
                ' | COMMON NAME | The value for species.common_name in your meta table: "'
                + core_dict["species.common_name"]
                + '" does not match the value that I am assigning to organism.common_name: "'
                + truth_dict["organism.common_name"]
                + '"'
            )

    meta_keys_to_remove = [
        "species.strain",
        "strain.type",
        ]
    for remove_key in meta_keys_to_remove:
        print(
            f"DELETE from meta WHERE meta_key='{remove_key}';",
            file=sql_out,
        )
            
    if args.verbose:
        # Print with the 'required' column if -v flag is present
        print("\n".join(f"{key:<28}   {val:<20}   {'required' if key in required_meta_keys else ''}" 
                    for key, val in truth_dict.items()))

