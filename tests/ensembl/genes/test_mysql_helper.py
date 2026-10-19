"""Tests for the shared MySQL connection helper."""

from unittest.mock import patch

import pymysql
from ensembl.genes.mysql_helper import get_mysql_connection
from pymysql.cursors import DictCursor


def test_get_mysql_connection_normalises_common_options():
    """Pass normalized connection options to PyMySQL."""
    expected_connection = object()

    with patch(
        "ensembl.genes.mysql_helper.pymysql.connect",
        return_value=expected_connection,
    ) as connect:
        result = get_mysql_connection(
            database=" registry ",
            host="db.example",
            port=3307,
            user="ensro",
            password="secret",
            cursorclass=DictCursor,
            connect_timeout=12,
        )

    assert result is expected_connection
    connect.assert_called_once_with(
        host="db.example",
        port=3307,
        user="ensro",
        database="registry",
        password="secret",
        cursorclass=DictCursor,
        connect_timeout=12,
    )


def test_get_mysql_connection_preserves_pymysql_errors():
    """Propagate connection errors to the caller."""
    with patch(
        "ensembl.genes.mysql_helper.pymysql.connect",
        side_effect=pymysql.MySQLError("connection failed"),
    ):
        try:
            get_mysql_connection(host="db.example", user="ensro")
        except pymysql.MySQLError as error:
            assert str(error) == "connection failed"
        else:
            raise AssertionError("MySQL errors should be propagated")
