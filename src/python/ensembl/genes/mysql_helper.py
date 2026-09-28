"""Shared MySQL connection helpers for the Ensembl genes tools.

The helper deliberately returns PyMySQL's normal connection object and lets
the caller choose a cursor class.  This keeps tuple-cursor and dictionary-
cursor users compatible while giving all scripts one connection entry point.
"""

from typing import TYPE_CHECKING, Any, TypeAlias

import pymysql
from pymysql.connections import Connection
from pymysql.cursors import Cursor

if TYPE_CHECKING:
    MySQLConnection: TypeAlias = Connection[Any]
else:
    # PyMySQL versions before 1.0 expose Connection as a non-subscriptable class.
    MySQLConnection: TypeAlias = Connection


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


__all__ = ["MySQLConnection", "get_mysql_connection"]
