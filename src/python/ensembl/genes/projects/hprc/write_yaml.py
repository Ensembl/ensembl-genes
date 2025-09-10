import os.path, sys, getopt, re
import pymysql
import requests
from ftplib import FTP

def get_population(accession):
    assembly_link = "https://www.ncbi.nlm.nih.gov/assembly/"+accession
    assembly_html = requests.get(assembly_link).text.splitlines()

    population = ""
    for line in assembly_html:
        biosample_regex =  re.search("\/biosample\/([A-Z0-9]+)\/\"", line)
        if biosample_regex:
            biosample_id = (biosample_regex.group(1))
            biosample_link = "https://www.ncbi.nlm.nih.gov/biosample/"+biosample_id
            biosample_html = requests.get(biosample_link).text.splitlines()
            
            for line in biosample_html:
                pop_desc_regex = re.search("Population Description\<\/th\>\<td\>([A-Za-z0-9 ]+)", line)
                population_regex = re.search("population=([A-Za-z0-9 ,]+)", line)
                race_regex = re.search("race\<\/th\>\<td\>([A-Za-z0-9 ]+)", line)
                if pop_desc_regex:
                    population += pop_desc_regex.group(1)
                elif population_regex:
                    population += population_regex.group(1)
                elif race_regex:
                    population += race_regex.group(1)
    
    population = population.title()
    population = population.replace('Usa', 'USA')
    return(population)

def check_for_file(species_name, prod_name, accession, source, file_type):
#This checks for 2 different types of files in the FTP
#The repeatmodeler file, e.g. https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_021951015.1/ensembl/variation/2022_10/vcf/Homo_sapiens-GCA_021951015.1-2022_10-gnomad.vcf.gz
    if file_type == "gnomad":
        ftp = FTP("ftp.ensembl.org")
        ftp_path = "https://ftp.ensembl.org/"
        path = (
            "pub/rapid-release/species/"
            + species_name
            + "/"
            + accession
            + "/"
            + source
            + "/variation/2022_10/vcf/"
        )
        file_name = species_name + "-" + accession + "-2022_10-gnomad.vcf.gz"

    elif file_type == "clinvar":
        ftp = FTP("ftp.ensembl.org")
        ftp_path = "https://ftp.ensembl.org/"
        path = (
            "pub/rapid-release/species/"
            + species_name
            + "/"
            + accession
            + "/"
            + source
            + "/variation/2022_10/vcf/"
        )
        file_name = species_name + "-" + accession + "-2022_10-clinvar.vcf.gz"
        
    ftp.login()
    try:
        ftp.cwd(path)
        if file_name in ftp.nlst():
            return ftp_path + path + file_name
        else:
            return 0
    except:
        return 0
    ftp.close()

if __name__ == "__main__":

    server_info = {
        "db_host": "mysql-ens-sta-5.ebi.ac.uk",
        "db_port": 4684,
        "db_user": "ensro",
        "db_pass": "",
    }

    db_list = open("hprc_cores.txt").readlines()
    
    yaml_out = open("hprc_species.yaml", "w")

    for db in db_list:
        db = db.strip()

        conn = pymysql.connect(
            host=server_info["db_host"],
            user=server_info["db_user"],
            passwd=server_info["db_pass"],
            port=server_info["db_port"],
            database=db
        )

        query = "select meta_key, meta_value from meta where meta_key in ('assembly.name', 'assembly.accession', 'genebuild.last_geneset_update', 'species.production_name')"
    
        cursor = conn.cursor()
        cursor.execute(query)
        info = cursor.fetchall()

        info_dict = {}
        for tuple in info:
            info_dict[tuple[0]] = tuple[1]
    
        if (info_dict["assembly.accession"] == "GCA_009914755.4"):
            assembly_provider = "T2T Consortium"
        else:
            assembly_provider = "UCSC Genomics Institute"

        ftp_base = ftp_base = "https://ftp.ensembl.org/pub/rapid-release/species"
        species_name = "Homo_sapiens"
        date = info_dict["genebuild.last_geneset_update"].replace("-", "_")
        variation_date = "2022_10" #need better way to get this
        uc_prod_name = (info_dict["species.production_name"]).capitalize()

        population = get_population(info_dict["assembly.accession"])

        gnomad_file = check_for_file(
            species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            "ensembl",
            "gnomad",
        )
        clinvar_file = check_for_file(
      	    species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            "ensembl",
            "clinvar",
        )
        
        yaml = (
            "- assembly: "
            + info_dict["assembly.name"]
            + "\n"
            )
        if population:
            yaml += (
            "  population: "
            + population
            + "\n"
            )
        yaml += (
            "  assembly_accession: "
            + info_dict["assembly.accession"]
            + "\n"
            + "  assembly_link: "
            + "https://www.ebi.ac.uk/ena/browser/view/"
            + info_dict["assembly.accession"]
            + "\n"
            + "  assembly_submitter: "
            + assembly_provider
            + "\n"
            + "  annotation_gtf: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/ensembl/geneset/"
            + date
            + "/"
            + species_name
            + "-"
            + info_dict["assembly.accession"]
            + "-"
            + date
            + "-genes.gtf.gz\n"
            + "  annotation_gff3: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/ensembl/geneset/"
            + date
            + "/"
            + species_name
            + "-"
            + info_dict["assembly.accession"]
            + "-"
            + date
            + "-genes.gff3.gz\n"
            "  proteins: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/ensembl/geneset/"
            + date
            + "/"
            + species_name
            + "-"
            + info_dict["assembly.accession"]
            + "-"
            + date
            + "-pep.fa.gz\n"
            "  transcripts: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/ensembl/geneset/"
            + date
            + "/"
            + species_name
            + "-"
            + info_dict["assembly.accession"]
            + "-"
            + date
            + "-cdna.fa.gz\n"
            )
        if clinvar_file:
            yaml += ("  variants_clinvar: "
            + clinvar_file
            + "\n"
            )
        if gnomad_file:
            yaml += ("  variants_gnomad: "
            + gnomad_file
            + "\n"
            )
        yaml += ("  ftp_dumps: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "\n"
            "  rapid_link: https://rapid.ensembl.org/"
            + uc_prod_name
            + "/Info/Index\n"        
        )
        
        print(yaml, file=yaml_out)
