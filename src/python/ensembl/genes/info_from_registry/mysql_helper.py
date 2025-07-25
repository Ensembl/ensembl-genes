import pymysql
import logging


def mysql_fetch_data(query, database, host, port, user, password, params=None):
    try:
        conn = pymysql.connect(
            host=host,
            user=user,
            port=port,
            password=password,
            database=database.strip(),
            cursorclass=pymysql.cursors.DictCursor,
        )
        with conn.cursor() as cursor:
            cursor.execute(query, params or ())
            info = cursor.fetchall()
        conn.close()
        return info

    except pymysql.Error as err:
        print(f"MySQL error: {err}")
        return []


def mysql_update(query, database, host, port, user, password, params=None):
    try:
        conn = pymysql.connect(
            host=host,
            user=user,
            port=port,
            password=password,
            database=database.strip(),
            cursorclass=pymysql.cursors.DictCursor,
        )
        with conn.cursor() as cursor:
            cursor.execute(query, params or ())
            conn.commit()
        return True

    except pymysql.Error as err:
        print(f"MySQL error: {err}")
        return False
