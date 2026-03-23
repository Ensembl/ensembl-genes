import pymysql
import json
import os

with open("src/python/ensembl/genes/projects/server_config.json") as f:
    conf = json.load(f)["meta_beta"]

conn = pymysql.connect(host=conf["db_host"], port=conf["db_port"], user=conf["db_user"], database=conf["db_name"], connect_timeout=5)
try:
    with conn.cursor() as cur:
        # Search for Aix sponsa and Anthonomus rubi UUIDs
        cur.execute("SELECT genome.genome_uuid, organism.scientific_name FROM genome JOIN organism USING(organism_id) WHERE organism.scientific_name IN ('Aix sponsa', 'Anthonomus rubi');")
        rows = cur.fetchall()
        
        with open("test_inputs.txt", "w") as out:
            for r in rows:
                print("Found target genome:", r[1], r[0])
                out.write(r[0] + "\n")
finally:
    conn.close()
