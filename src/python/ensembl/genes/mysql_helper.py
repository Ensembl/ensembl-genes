"""Shared MySQL connection helpers for the Ensembl genes tools.

The helper deliberately returns PyMySQL's normal connection object and lets
the caller choose a cursor class.  This keeps tuple-cursor and dictionary-
cursor users compatible while giving all scripts one connection entry point.
"""

import logging
from typing import TYPE_CHECKING, Any, TypeAlias

import pymysql
from pymysql.connections import Connection
from pymysql.cursors import Cursor, DictCursor

if TYPE_CHECKING:
    MySQLConnection: TypeAlias = Connection[Any]
else:
    # PyMySQL versions before 1.0 expose Connection as a non-subscriptable class.
    MySQLConnection: TypeAlias = Connection

logger = logging.getLogger(__name__)


def get_mysql_connection(
    *,
    database: str | None = None,
    host: str,
    port: int = 3306,
    user: str,
    password: str | None = None,
    cursorclass: type[Cursor] | None = None,
    **kwargs: Any,
) -> MySQLConnection:
    """Create a PyMySQL connection with consistent connection options.

    ``cursorclass`` is optional because existing scripts use both PyMySQL's
    tuple cursors and ``DictCursor``.  Additional PyMySQL options can be
    passed through ``kwargs`` for future callers without changing this API.
    Connection errors are intentionally propagated to the caller so each
    script can apply its own error handling and logging policy.
    """

    options: dict[str, Any] = {
        "host": host,
        "port": int(port),
        "user": user,
        **kwargs,
    }
    if database is not None:
        options["database"] = database.strip()
    if password is not None:
        options["password"] = password
    if cursorclass is not None:
        options["cursorclass"] = cursorclass

    return pymysql.connect(**options)


def mysql_get_connection(
    database: str, host: str, port: int, user: str, password: str
) -> MySQLConnection | None:
    """Open a dictionary-cursor connection for legacy registry callers.

    New code should use :func:`get_mysql_connection` when it needs control
    over cursor type or connection error handling.
    """

    try:
        return get_mysql_connection(
            database=database,
            host=host,
            port=port,
            user=user,
            password=password,
            cursorclass=DictCursor,
        )
    except pymysql.Error as error:
        logger.error("MySQL connection failed: %s", error)
        return None


def mysql_fetch_data(  # pylint: disable=too-many-arguments
    query: str,
    database: str,
    host: str,
    port: int,
    user: str,
    password: str = "",
    params: tuple[Any, ...] | list[Any] | None = None,
) -> list[dict[str, Any]]:
    """Run a SELECT query and return rows as dictionaries.

    Database errors are logged and represented by an empty result, matching
    the behavior of the former registry-local helper.
    """

    connection: MySQLConnection | None = None
    try:
        connection = get_mysql_connection(
            database=database,
            host=host,
            port=port,
            user=user,
            password=password,
            cursorclass=DictCursor,
        )
        with connection.cursor() as cursor:
            cursor.execute(query, params or ())
            return list(cursor.fetchall())
    except pymysql.Error as error:
        logger.error("MySQL fetch failed: %s", error)
        return []
    finally:
        if connection is not None:
            connection.close()


def mysql_update(  # pylint: disable=too-many-arguments
    query: str,
    database: str,
    host: str,
    port: int,
    user: str,
    password: str = "",
    params: tuple[Any, ...] | list[Any] | None = None,
) -> bool:
    """Run an INSERT, UPDATE, or DELETE query and commit it."""

    connection: MySQLConnection | None = None
    try:
        connection = get_mysql_connection(
            database=database,
            host=host,
            port=port,
            user=user,
            password=password,
            cursorclass=DictCursor,
        )
        with connection.cursor() as cursor:
            cursor.execute(query, params or ())
        connection.commit()
        return True
    except pymysql.Error as error:
        logger.error("MySQL update failed: %s", error)
        return False
    finally:
        if connection is not None:
            connection.close()


__all__ = [
    "MySQLConnection",
    "get_mysql_connection",
    "mysql_fetch_data",
    "mysql_get_connection",
    "mysql_update",
]
