import pymysql
import json
with open('src/python/ensembl/genes/projects/server_config.json') as f:
    conf = json.load(f)['gb1']
conn = pymysql.connect(host=conf['db_host'], port=conf['db_port'], user=conf['db_user'], database='gb_assembly_metadata')
with conn.cursor(pymysql.cursors.DictCursor) as cur:
    cur.execute("SELECT COUNT(*) as c FROM genebuild_status WHERE gb_status='pre_released'")
    print('Total Pre-released:', cur.fetchone()['c'])
    
    cur.execute("SELECT gb_status, COUNT(*) as c from genebuild_status group by gb_status")
    print('Statuses:', cur.fetchall())

    cur.execute("SELECT bioproject_id, COUNT(*) as c FROM bioproject GROUP BY bioproject_id LIMIT 10")
    print('\nBioProjects Sample:', cur.fetchall())

    cur.execute("SELECT group_name, group_type, COUNT(*) as c FROM custom_group GROUP BY group_name, group_type")
    print('\nCustom Groups:', cur.fetchall())
    
    cur.execute("SELECT gca_accession FROM genebuild_status WHERE gb_status='pre_released' LIMIT 1")
    acc = cur.fetchone()
    if acc:
        print("\nExample accession:", acc['gca_accession'])
