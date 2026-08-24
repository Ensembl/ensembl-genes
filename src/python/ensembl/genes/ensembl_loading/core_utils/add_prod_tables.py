"""
Sync production-controlled tables from the production DB to a core DB.
"""

import argparse
import json
import logging
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path

import pymysql

logger = logging.getLogger(__name__)


TABLES_TO_SYNC = {
    "master_external_db": "external_db",
    "master_attrib_type": "attrib_type",
    "master_misc_set": "misc_set",
    "master_biotype": "biotype",
    "master_unmapped_reason": "unmapped_reason",
}
DEFAULT_CONFIG_PATH = (
    Path(__file__).resolve().parent.parent / "config" / "production_db_conf.json"
)


@dataclass(frozen=True)
class ProductionDbConfig:
    """Connection settings for the production database."""

    host: str
    port: int
    user: str
    password: str


def load_config(path: str | Path) -> ProductionDbConfig:
    """Load and validate production DB connection settings from JSON."""
    config_path = Path(path)
    try:
        with config_path.open("r", encoding="utf-8") as handle:
            values = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read production DB config {config_path}") from exc

    if not isinstance(values, dict):
        raise TypeError("Production DB config must contain a JSON object")

    required = ("host", "port", "user", "password")
    missing = [key for key in required if key not in values]
    if missing:
        raise ValueError("Production DB config is missing: " + ", ".join(missing))
    if not isinstance(values["host"], str) or not isinstance(values["user"], str):
        raise TypeError("Production DB config host and user must be strings")
    if not isinstance(values["password"], str):
        raise TypeError("Production DB config password must be a string")
    if isinstance(values["port"], bool) or not isinstance(values["port"], int):
        raise TypeError("Production DB config port must be an integer")
    if not 1 <= values["port"] <= 65535:
        raise ValueError("Production DB config port must be between 1 and 65535")

    return ProductionDbConfig(
        host=values["host"],
        port=values["port"],
        user=values["user"],
        password=values["password"],
    )


def load_default_config() -> ProductionDbConfig:
    """Load the bundled production database configuration."""
    return load_config(DEFAULT_CONFIG_PATH)


def get_connection(
    host: str,
    port: int,
    user: str,
    password: str,
    db: str,
) -> pymysql.Connection:
    """Connect to a MySQL database."""
    logger.info(
        "Connecting to database %s on %s:%s as %s",
        db,
        host,
        port,
        user,
    )

    return pymysql.connect(
        host=host,
        port=port,
        user=user,
        password=password,
        database=db,
        cursorclass=pymysql.cursors.DictCursor,
    )


def fetch_table_rows(
    prod_conn: pymysql.Connection,
    table: str,
) -> list[dict]:
    """Fetch current rows from a production table when supported."""
    logger.debug("Fetching rows from production table %s", table)

    with closing(prod_conn.cursor()) as cur:
        cur.execute(f"SHOW COLUMNS FROM `{table}`")
        production_columns = [row["Field"] for row in cur.fetchall()]
        current_filter = (
            " WHERE `is_current` = 1" if "is_current" in production_columns else ""
        )
        cur.execute(f"SELECT * FROM `{table}`{current_filter}")
        rows = cur.fetchall()

    logger.info(
        "Fetched %d rows from production table %s",
        len(rows),
        table,
    )

    return rows


def fetch_core_columns(
    core_conn: pymysql.Connection,
    table: str,
) -> list[str]:
    """Return column names present in a core table."""
    with closing(core_conn.cursor()) as cur:
        cur.execute(f"SHOW COLUMNS FROM `{table}`")
        rows = cur.fetchall()
        return [row["Field"] if isinstance(row, dict) else row[0] for row in rows]


def sync_table(
    prod_conn: pymysql.Connection,
    core_conn: pymysql.Connection,
    prod_table: str,
    core_table: str,
) -> None:
    """Sync a single production table to the core database."""
    logger.info(
        "Syncing %s -> %s",
        prod_table,
        core_table,
    )

    rows = fetch_table_rows(prod_conn, prod_table)
    core_columns = fetch_core_columns(core_conn, core_table)
    columns: list[str] = []
    excluded_columns: list[str] = []
    if rows:
        production_columns = list(rows[0].keys())
        columns = [column for column in core_columns if column in production_columns]
        excluded_columns = [
            column for column in production_columns if column not in columns
        ]
        if not columns:
            raise ValueError(f"No shared columns between {prod_table} and {core_table}")

    try:
        with closing(core_conn.cursor()) as cur:
            logger.debug("Truncating core table %s", core_table)
            cur.execute(f"TRUNCATE TABLE `{core_table}`")

            if rows:
                if excluded_columns:
                    logger.info(
                        "Ignoring production-only columns for %s: %s",
                        prod_table,
                        ", ".join(excluded_columns),
                    )

                col_list = ", ".join(f"`{column}`" for column in columns)
                placeholders = ", ".join(["%s"] * len(columns))

                values = [tuple(row[column] for column in columns) for row in rows]

                sql = (
                    f"INSERT INTO `{core_table}` "
                    f"({col_list}) "
                    f"VALUES ({placeholders})"
                )

                logger.debug(
                    "Inserting %d rows into core table %s",
                    len(rows),
                    core_table,
                )
                cur.executemany(sql, values)

        core_conn.commit()

        logger.info(
            "Successfully synced %s -> %s (%d rows)",
            prod_table,
            core_table,
            len(rows),
        )

    except Exception:
        core_conn.rollback()
        logger.exception(
            "Failed to sync %s -> %s",
            prod_table,
            core_table,
        )
        raise


def sync_tables(
    prod_conn: pymysql.Connection,
    core_conn: pymysql.Connection,
) -> None:
    """Sync all production tables to the core database."""
    logger.info(
        "Starting sync of %d production-controlled tables",
        len(TABLES_TO_SYNC),
    )

    for prod_table, core_table in TABLES_TO_SYNC.items():
        sync_table(
            prod_conn,
            core_conn,
            prod_table,
            core_table,
        )

    logger.info("All production-controlled tables synced successfully")


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Sync production tables to core database"
    )

    parser.add_argument(
        "--db_host",
        required=True,
        help="Production DB host",
    )
    parser.add_argument(
        "--db_port",
        type=int,
        required=True,
        help="Production DB port",
    )
    parser.add_argument(
        "--db_user",
        required=True,
        help="Production DB user",
    )
    parser.add_argument(
        "--db_password",
        required=True,
        help="Production DB password (read access)",
    )

    parser.add_argument(
        "--core_db_host",
        required=True,
        help="Core DB host",
    )
    parser.add_argument(
        "--core_db_port",
        type=int,
        required=True,
        help="Core DB port",
    )
    parser.add_argument(
        "--core_db_user",
        required=True,
        help="Core DB user",
    )
    parser.add_argument(
        "--core_db_password",
        required=True,
        help="Core DB password (write access)",
    )
    parser.add_argument(
        "--core_db_name",
        required=True,
        help="Core DB name",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    prod_cfg = {
        "host": args.db_host,
        "port": args.db_port,
        "user": args.db_user,
        "password": args.db_password,
        "db": "ensembl_production",
    }

    core_cfg = {
        "host": args.core_db_host,
        "port": args.core_db_port,
        "user": args.core_db_user,
        "password": args.core_db_password,
        "db": args.core_db_name,
    }

    prod_conn = None
    core_conn = None

    try:
        prod_conn = get_connection(**prod_cfg)
        core_conn = get_connection(**core_cfg)

        sync_tables(prod_conn, core_conn)

    except Exception:
        logger.exception("Production table sync failed")
        raise

    finally:
        if prod_conn is not None:
            prod_conn.close()
            logger.debug("Closed production database connection")

        if core_conn is not None:
            core_conn.close()
            logger.debug("Closed core database connection")


if __name__ == "__main__":
    main()
