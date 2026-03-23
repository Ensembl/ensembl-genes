import pymysql

conn = pymysql.connect(
    host="mysql-ens-meta-prod-1.ebi.ac.uk", user="ensro", port=4483, database="ensembl_metadata_qrp"
)
with conn.cursor() as cursor:
    cursor.execute("SHOW TABLES LIKE '%related%';")
    print("Tables like related:", cursor.fetchall())
    cursor.execute("SHOW TABLES LIKE '%assembly%';")
    print("Tables like assembly:", cursor.fetchall())
    cursor.execute("SHOW COLUMNS FROM assembly;")
    print("Assembly cols:", [row[0] for row in cursor.fetchall()])
conn.close()
