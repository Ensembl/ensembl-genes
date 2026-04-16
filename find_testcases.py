import pymysql
import json

with open('src/python/ensembl/genes/projects/server_config.json') as f:
    conf = json.load(f)

# Connect to metadata
meta_conf = conf['meta_beta']
c_meta = pymysql.connect(host=meta_conf['db_host'], port=meta_conf['db_port'], user=meta_conf['db_user'], database='ensembl_metadata_qrp')
cur_meta = c_meta.cursor(pymysql.cursors.DictCursor)

# Connect to gb_schema
gb_conf = conf['gb1']
c_gb = pymysql.connect(host=gb_conf['db_host'], port=gb_conf['db_port'], user=gb_conf['db_user'], database='gb_assembly_metadata')
cur_gb = c_gb.cursor(pymysql.cursors.DictCursor)

# 1. GB-only pre_released genome
cur_gb.execute("SELECT gs.gca_accession FROM genebuild_status gs JOIN bioproject bp ON gs.assembly_id=bp.assembly_id JOIN main_bioproject mb ON bp.bioproject_id=mb.bioproject_id WHERE gs.gb_status='pre_released' AND mb.bioproject_name='CBP' LIMIT 1")
gb_only = cur_gb.fetchone()
print("GB-only pre_released (CBP):", gb_only)

# 2. Released-only or Both genome
# Let's find a CBP genome UUID in metadata
cur_meta.execute("SELECT g.genome_uuid, a.accession FROM genome g JOIN assembly a ON g.assembly_id=a.assembly_id JOIN dataset d ON g.genome_id=d.genome_id WHERE d.dataset_type_id=1 LIMIT 10")
metas = cur_meta.fetchall()

accs = [m['accession'] for m in metas]
if tuple(accs):
    cur_gb.execute(f"SELECT gca_accession, gb_status FROM genebuild_status WHERE gca_accession IN {tuple(accs)} AND gb_status='pre_released'")
    both = cur_gb.fetchall()
    print("In Both (pre-released in GB, released in Meta):", both)

    # find any other for released only
    cur_gb.execute(f"SELECT gca_accession FROM genebuild_status WHERE gca_accession IN {tuple(accs)}")
    in_gb = set(r['gca_accession'] for r in cur_gb.fetchall())
    released_only_acc = set(accs) - in_gb
    print("Released only accessions (not in GB at all):", released_only_acc)

print("UUIDs mapping:", metas)

