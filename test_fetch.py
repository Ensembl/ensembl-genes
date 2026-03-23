import pymysql
import json

with open("src/python/ensembl/genes/projects/server_config.json") as f:
    conf = json.load(f)["meta_beta"]

conn = pymysql.connect(host=conf["db_host"], port=conf["db_port"], user=conf["db_user"], database=conf["db_name"], connect_timeout=5)
try:
    with conn.cursor() as cur:
        cur.execute("SELECT genome_uuid FROM genome LIMIT 1")
        uuid = cur.fetchone()[0]
        with open("dummy_input.txt", "w") as out:
            out.write(uuid + "\n")
        print("Found UUID:", uuid)
finally:
    conn.close()
