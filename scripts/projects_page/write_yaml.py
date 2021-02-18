import os.path, sys, getopt
import pymysql
import requests
import argparse
from ftplib import FTP
import re
from collections import OrderedDict

def check_for_repeatmodeler (species_name, accession):
  ftp = FTP('ftp.ebi.ac.uk')
  ftp.login()
  try:
    ftp.cwd("pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/"+species_name+"/")
    if accession+".repeatmodeler.fa" in ftp.nlst():
      return "http://ftp.ebi.ac.uk/pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/"+species_name+"/"+accession+".repeatmodeler.fa"
    else:
      return 0
  except:
    return 0
        
def mysql_fetch_data (query,database,host,port,user,password):
   try:
      conn = pymysql.connect(host=host, 
                             user=user, 
                             passwd=password, 
                             port=port,
                             database=database.strip())

      cursor = conn.cursor()
      cursor.execute(query)
      info = cursor.fetchall()

   except mysql.connector.Error as err:
      if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
         print("Something is wrong with your user name or password")
      elif err.errno == errorcode.ER_BAD_DB_ERROR:
         print("Database does not exist on the server")
      else:
         print(err)

   cursor.close()
   conn.close()
   return info

def write_yaml(info_dict,icon,yaml_out,project,use_server):
  #there are species on main for which the upper case production name is used in the url instead of the upper case species name
  prod_url_list = ['bos_taurus_hybrid', 'bos_indicus_hybrid']

  assembly_name = info_dict['assembly.name'].replace(" ", "_")
  date = info_dict['genebuild.last_geneset_update'].replace("-", "_")
  species_name = info_dict['species.scientific_name'].replace(" ", "_")
  if ('species.strain' in info_dict): 
    if (info_dict['species.strain'] != 'reference'):
      info_dict['species.scientific_name'] = info_dict['species.scientific_name']+" ("+info_dict['species.strain']+")"
  lc_species_name = (info_dict['species.scientific_name'].replace(" ", "_")).lower()
  uc_prod_name = (info_dict['species.production_name']).capitalize()
   
  assembly_report_url = "https://www.ncbi.nlm.nih.gov/assembly/"+info_dict['assembly.accession']
  assembly_report_response = requests.get(assembly_report_url)
  submitter_match = re.search("Submitter: </dt><dd>([^<]*)</dd><dt>", assembly_report_response.text)
  if submitter_match:
    submitter = submitter_match.group(1)
  else:
    submitter = 'unknown'

  if (use_server == 'rapid'): 
    ftp_base = "http://ftp.ensembl.org/pub/rapid-release/species"

    yaml = "- species: "+info_dict['species.scientific_name']+"\n"

    if project in ('vgp', 'dtol'):
      yaml += "  image: "+icon+"\n"
    else:
      yaml += "  submitted_by: "+submitter+"\n"

    yaml += "  accession: "+info_dict['assembly.accession']+"\n"
    yaml += "  annotation_gtf: "+ftp_base+"/"+species_name+"/"+info_dict['assembly.accession']+"/geneset/"+date+"/"+species_name+"-"+info_dict['assembly.accession']+"-"+date+"-genes.gff3.gz\n"
    yaml += "  annotation_gff3: "+ftp_base+"/"+species_name+"/"+info_dict['assembly.accession']+"/geneset/"+date+"/"+species_name+"-"+info_dict['assembly.accession']+"-"+date+"-genes.gtf.gz\n"
    yaml += "  proteins: "+ftp_base+"/"+species_name+"/"+info_dict['assembly.accession']+"/geneset/"+date+"/"+species_name+"-"+info_dict['assembly.accession']+"-"+date+"-pep.fa.gz\n"
    yaml += "  transcripts: "+ftp_base+"/"+species_name+"/"+info_dict['assembly.accession']+"/geneset/"+date+"/"+species_name+"-"+info_dict['assembly.accession']+"-"+date+"-cdna.fa.gz\n"
    yaml += "  softmasked_genome: "+ftp_base+"/"+species_name+"/"+info_dict['assembly.accession']+"/genome/"+species_name+"-"+info_dict['assembly.accession']+"-softmasked.fa.gz\n"
    
    rm_file = check_for_repeatmodeler(lc_species_name,info_dict['assembly.accession'])
    if (rm_file):
      yaml += "  repeat_library: "+rm_file+"\n"

    yaml += "  ftp_dumps: "+ftp_base+"/"+species_name+"/"+info_dict['assembly.accession']+"\n"
 
    main_species_url = "http://www.ensembl.org/info/about/species.html"
    main_species_response = requests.get(main_species_url)
    if info_dict['assembly.accession'] in main_species_response.text:
      yaml += "  ensembl_link: https://www.ensembl.org/"+species_name+"/Info/Index\n"
    else:
      yaml += "  rapid_link: https://rapid.ensembl.org/"+uc_prod_name+"/Info/Index\n"

  elif (use_server == 'main'):
    release = "release-"+info_dict['schema_version']
    ftp_base = "http://ftp.ensembl.org/pub/"+release

    yaml = "- species: "+info_dict['species.scientific_name']+"\n"
    if project in ('vgp', 'dtol'):
      yaml += "  image: "+icon+"\n"
    else:
      yaml += "  submitted_by: "+submitter+"\n"

    yaml += "  accession: "+info_dict['assembly.accession']+"\n"
    yaml += "  annotation_gtf: "+ftp_base+"/gtf/"+info_dict['species.production_name']+"/"+uc_prod_name+"."+assembly_name+".gtf.gz\n"
    yaml += "  annotation_gff3: "+ftp_base+"/gff3/"+info_dict['species.production_name']+"/"+uc_prod_name+"."+assembly_name+".gff3.gz\n"
    yaml += "  proteins: "+ftp_base+"/fasta/"+info_dict['species.production_name']+"/pep/"+uc_prod_name+"."+assembly_name+".pep.all.fa.gz\n"
    yaml += "  transcripts: "+ftp_base+"/fasta/"+info_dict['species.production_name']+"/cdna/"+uc_prod_name+"."+assembly_name+".cdna.all.fa.gz\n"
    yaml += "  softmasked_genome: "+ftp_base+"/fasta/"+info_dict['species.production_name']+"/dna/"+uc_prod_name+"."+assembly_name+".dna_sm.toplevel.fa.gz\n"

    rm_file = check_for_repeatmodeler(lc_species_name,info_dict['assembly.accession'])
    if (rm_file):
      yaml += "  repeat_library: "+rm_file+"\n"

    yaml += "  ftp_dumps: "+ftp_base+"\n"
    if (info_dict['species.production_name'] in prod_url_list):
      yaml += "  ensembl_link: https://www.ensembl.org/"+uc_prod_name+"/Info/Index\n"
    else:
      yaml += "  ensembl_link: https://www.ensembl.org/"+species_name+"/Info/Index\n"

  print (yaml, file=yaml_out)

if __name__ == '__main__':

  server_dict = {'rapid' : {'db_host' : 'mysql-ens-sta-5',
                            'db_port' : 4672,
                            'db_user' : 'ensro',
                            'db_pass' : ''},
                 'main' :  {'db_host' : 'mysql-ens-mirror-1',
                            'db_port' : 4240,
                            'db_user' : 'ensro',
                            'db_pass' : ''}
                 }
  
  icon_list = open('icons.txt').readlines()
  icon_dict = {}
  for line in icon_list:
    icon_dict[line.split()[0]] = line.split()[1]
    
  parser = argparse.ArgumentParser(description='Create species.yaml file for a given project page.')
  parser.add_argument('-f','--db_file', help='Name for file containing list of VGP databases on the Rapid Release or Main server', required=True)
  parser.add_argument('-p','--project', choices=['aquafaang','bovreg','dtol','geneswitch','vgp'], help='Name of the project this set of database belongs to', required=True)
  args = parser.parse_args()

  db_list = open(args.db_file).readlines()
  project = args.project
  yaml_out = open(project+"_species.yaml", 'w')
  
  sorted_db_list = sorted(db_list)
  
  if (project == 'aquafaang'):
    # move danio rerio reference dbs to end of the list
    for db in (sorted_db_list):
      if ('danio_rerio_core' in db):
        sorted_db_list.append(sorted_db_list.pop(sorted_db_list.index(db))) 
        
  if (project == 'bovreg'):
    # move bos taurus reference dbs to top of the list
    for db in (sorted_db_list):
      if ('bos_taurus_core' in db):
        sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))

  if (project == 'geneswitch'):
    # move sus scrofa or gallus gallus reference dbs to top of the list
    for db in (sorted_db_list):
      if ('sus_scrofa_core' in db or 'gallus_gallus_core' in db):
        sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))

  for db in (sorted_db_list):
    db = db.strip()
    #check if the db is on rapid
    conn = pymysql.connect(host=server_dict['rapid']['db_host'],
                           user=server_dict['rapid']['db_user'],
                           passwd=server_dict['rapid']['db_pass'],
                           port=server_dict['rapid']['db_port'])
      
    cur = conn.cursor()
    exists_query = "SHOW DATABASES LIKE '"+db.strip()+"'"
    exists_rapid = cur.execute(exists_query)
      
    if (exists_rapid):
      use_server = 'rapid'
    #if it's not on the rapid server check main
    else:
      conn = pymysql.connect(host=server_dict['main']['db_host'],
                           user=server_dict['main']['db_user'],
                           passwd=server_dict['main']['db_pass'],
                           port=server_dict['main']['db_port'])

      cur = conn.cursor()
      exists_query = "SHOW DATABASES LIKE '"+db.strip()+"'"
      exists_main = cur.execute(exists_query)
      if (exists_main):
        use_server = 'main'

    if (use_server):
      # retrieve the species name, assembly accession and assembly name from the database
      info_query = "SELECT meta_key,meta_value FROM meta WHERE meta_key in ('species.scientific_name','assembly.accession','assembly.name','species.production_name','species.strain','schema_version','genebuild.last_geneset_update')"
      info = mysql_fetch_data(info_query,db,server_dict[use_server]['db_host'],server_dict[use_server]['db_port'],server_dict[use_server]['db_user'],server_dict[use_server]['db_pass'])
        
      info_dict = {}
      for tuple in info:
        info_dict[tuple[0]] = tuple[1]
          
      # retrieve the species classification info from the database. in order to assign an icon 
      class_query = "SELECT meta_value FROM meta WHERE meta_key='species.classification'"
      classifications = mysql_fetch_data(class_query,db,server_dict[use_server]['db_host'],server_dict[use_server]['db_port'],server_dict[use_server]['db_user'],server_dict[use_server]['db_pass'])
      
      class_list = []
      for tuple in classifications:
        class_list.append(tuple[0])

      icon = "Metazoa.png"
      chordate = 0
      if "Chordata" in class_list:
        chordate = 1
      for classification in class_list:
        try:
          icon = icon_dict[classification]
          break
        except KeyError:
          continue
      if chordate and icon == "Metazoa.png":
        icon = "Chordates.png"

      write_yaml(info_dict,icon,yaml_out,project,use_server)
    else:
      print ("Could not find database "+db.strip()+" on mirror or rapid release servers!\n")
