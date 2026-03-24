import pymysql
import json
import os

with open("src/python/ensembl/genes/projects/server_config.json") as f:
    conf = json.load(f)["meta_beta"]

conn = pymysql.connect(host=conf["db_host"], port=conf["db_port"], user=conf["db_user"], database=conf["db_name"], connect_timeout=5)
try:
    with conn.cursor() as cur:
        # Search for HPRC UUIDs
        cur.execute("SELECT genome.genome_uuid, assembly.accession FROM genome JOIN assembly USING(assembly_id) WHERE assembly.accession IN ('GCA_018467165.1', 'GCA_009914755.4');")
        rows = cur.fetchall()
        
        with open("test_inputs_hprc.txt", "w") as out:
            for r in rows:
                print("Found target genome:", r[1], r[0])
                out.write(r[0] + "\n")
finally:
    conn.close()
