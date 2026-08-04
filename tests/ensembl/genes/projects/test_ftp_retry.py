# Copyright 2026 EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for EnsemblFTP retry logic — all FTP and sleep calls are mocked."""

import socket
from ftplib import FTP, error_perm, error_temp
from unittest.mock import MagicMock, call, patch

import pytest

from ensembl.genes.projects.ftp_client import EnsemblFTP, _TRANSIENT_FTP_ERRORS


def make_client(max_retries: int = 2) -> EnsemblFTP:
    """Return an EnsemblFTP with FTP connections mocked out."""
    with (
        patch.object(EnsemblFTP, "_connect_ensembl"),
        patch.object(EnsemblFTP, "_connect_ebi"),
    ):
        client = EnsemblFTP(timeout=5, max_retries=max_retries)
    # Install mock FTP objects
    client.ensembl_ftp = MagicMock(spec=FTP)
    client.ebi_ftp = MagicMock(spec=FTP)
    # Prevent real sleeps
    client._sleep = MagicMock()
    return client


class TestGetConnection:
    """Group A: connection resolution tests."""

    def test_get_ebi_connection(self):
        c = make_client()
        assert c._get_connection("ebi") is c.ebi_ftp

    def test_get_ensembl_connection(self):
        c = make_client()
        assert c._get_connection("ensembl") is c.ensembl_ftp

    def test_unknown_identifier_raises_value_error(self):
        c = make_client()
        with pytest.raises(ValueError, match="Unknown FTP connection identifier"):
            c._get_connection("mystery_server")

    def test_ebi_and_ensembl_not_mixed(self):
        c = make_client()
        assert c._get_connection("ebi") is not c.ensembl_ftp
        assert c._get_connection("ensembl") is not c.ebi_ftp


class TestRetrySuccess:
    """Group B: correct retry and reconnect behaviour."""

    def test_successful_first_call_no_reconnect(self):
        """When the operation succeeds on the first attempt, no reconnect occurs."""
        c = make_client(max_retries=2)
        calls = []

        def _op(conn):
            calls.append(conn)
            return "ok"

        with patch.object(c, "_reconnect") as mock_reconnect:
            result = c._retry("ebi", _op)

        assert result == "ok"
        assert len(calls) == 1
        mock_reconnect.assert_not_called()
        c._sleep.assert_not_called()

    def test_transient_error_triggers_reconnect_and_retry(self):
        """A transient error on attempt 1 causes reconnect + sleep; attempt 2 succeeds."""
        c = make_client(max_retries=2)
        attempt = [0]

        def _op(conn):
            attempt[0] += 1
            if attempt[0] == 1:
                raise ConnectionResetError("broken pipe")
            return "recovered"

        with patch.object(c, "_reconnect") as mock_reconnect:
            result = c._retry("ebi", _op)

        assert result == "recovered"
        mock_reconnect.assert_called_once_with("ebi")
        c._sleep.assert_called_once()

    def test_replacement_connection_is_used_after_reconnect(self):
        """After reconnect, the operation runs on the NEW connection object."""
        c = make_client(max_retries=2)
        new_ftp = MagicMock(spec=FTP)
        connections_seen = []

        def _op(conn):
            connections_seen.append(conn)
            if len(connections_seen) == 1:
                raise EOFError("dropped")
            return "done"

        def replace_ebi():
            c.ebi_ftp = new_ftp

        with patch.object(c, "_reconnect", side_effect=lambda w: replace_ebi()):
            c._retry("ebi", _op)

        assert connections_seen[1] is new_ftp

    def test_original_connection_not_called_after_reconnect(self):
        """The original FTP object is never used after a successful reconnect."""
        c = make_client(max_retries=2)
        original_ftp = c.ebi_ftp
        new_ftp = MagicMock(spec=FTP)
        call_log = []

        def _op(conn):
            call_log.append(id(conn))
            if conn is original_ftp:
                raise ConnectionResetError("gone")
            return "ok"

        def replace_ebi():
            c.ebi_ftp = new_ftp

        with patch.object(c, "_reconnect", side_effect=lambda w: replace_ebi()):
            c._retry("ebi", _op)

        # Only the original connection was used on the failing attempt
        assert call_log[0] == id(original_ftp)
        # The second call used the new connection
        assert call_log[1] == id(new_ftp)

    def test_cwd_and_nlst_both_repeated_after_reconnect(self):
        """An operation that calls cwd then nlst repeats both calls on retry."""
        c = make_client(max_retries=2)
        first_ftp = c.ebi_ftp
        new_ftp = MagicMock(spec=FTP)
        new_ftp.nlst.return_value = ["file.txt"]
        call_sequence = []

        def _op(conn):
            conn.cwd("/")
            conn.cwd("/some/path")
            result = conn.nlst()
            return result

        first_ftp.cwd.side_effect = None
        first_ftp.nlst.side_effect = EOFError("dropped during nlst")

        def replace_ebi():
            c.ebi_ftp = new_ftp

        with patch.object(c, "_reconnect", side_effect=lambda w: replace_ebi()):
            files = c._retry("ebi", _op)

        assert files == ["file.txt"]
        # nlst on the new connection was called
        new_ftp.nlst.assert_called_once()
        # cwd was called on the new connection with both paths
        assert new_ftp.cwd.call_count == 2

    def test_all_attempts_exhausted_raises_last_exception(self):
        """After max_retries attempts, the last exception is re-raised."""
        c = make_client(max_retries=2)

        def _op(conn):
            raise ConnectionResetError("always fails")

        with patch.object(c, "_reconnect"):
            with pytest.raises(ConnectionResetError, match="always fails"):
                c._retry("ebi", _op)

    def test_error_perm_not_retried(self):
        """Permanent FTP errors (error_perm) propagate without retry."""
        c = make_client(max_retries=2)
        call_count = [0]

        def _op(conn):
            call_count[0] += 1
            raise error_perm("550 No such file")

        with patch.object(c, "_reconnect") as mock_reconnect:
            with pytest.raises(error_perm):
                c._retry("ebi", _op)

        # error_perm is not in _TRANSIENT_FTP_ERRORS, so only one attempt runs
        assert call_count[0] == 1
        mock_reconnect.assert_not_called()

    def test_arguments_forwarded_to_operation(self):
        """Extra args/kwargs are forwarded to the operation on each attempt."""
        c = make_client(max_retries=2)
        received = []

        def _op(conn, path, mode="r"):
            received.append((path, mode))
            if len(received) == 1:
                raise ConnectionResetError("oops")
            return received

        with patch.object(c, "_reconnect"):
            c._retry("ebi", _op, "/mypath", mode="w")

        assert all(r == ("/mypath", "w") for r in received)

    def test_reconnect_failure_becomes_last_exc(self):
        """If reconnect itself fails, the retry loop continues and the next
        operation attempt runs.  With max_retries=2, after attempt 1 fails
        (EOFError), reconnect fails (RuntimeError), attempt 2 runs the
        operation again and raises EOFError — which is the final exception."""
        c = make_client(max_retries=2)

        def _op(conn):
            raise EOFError("initial failure")

        reconnect_error = RuntimeError("reconnect completely failed")

        with patch.object(c, "_reconnect", side_effect=reconnect_error):
            # After reconnect failure last_exc becomes the reconnect error,
            # but then attempt 2 runs the operation and raises EOFError again,
            # which becomes the final last_exc that is raised.
            with pytest.raises(EOFError, match="initial failure"):
                c._retry("ebi", _op)


class TestRetrySleep:
    """Group C: sleep injection."""

    def test_sleep_called_once_on_single_retry(self):
        """_sleep is called exactly once when there is one retry after one failure."""
        c = make_client(max_retries=2)
        attempt = [0]

        def _op(conn):
            attempt[0] += 1
            if attempt[0] < 2:
                raise ConnectionResetError("first attempt")
            return "ok"

        with patch.object(c, "_reconnect"):
            c._retry("ebi", _op)

        c._sleep.assert_called_once_with(c.retry_sleep)

    def test_sleep_not_called_on_success(self):
        """_sleep is never called when the first attempt succeeds."""
        c = make_client(max_retries=2)
        c._retry("ebi", lambda conn: "ok")
        c._sleep.assert_not_called()

    def test_sleep_not_called_after_last_failed_attempt(self):
        """_sleep is not called on the final exhausted attempt (no more retries)."""
        c = make_client(max_retries=2)

        def _op(conn):
            raise EOFError("always fails")

        with patch.object(c, "_reconnect"):
            with pytest.raises(EOFError):
                c._retry("ebi", _op)

        # With max_retries=2 there is only 1 retry window (between attempt 1 and 2)
        assert c._sleep.call_count == 1


class TestMaxRetries:
    """Validate max_retries semantics."""

    def test_max_retries_one_means_one_attempt(self):
        """max_retries=1 means exactly one attempt (no retries)."""
        c = make_client(max_retries=1)
        call_count = [0]

        def _op(conn):
            call_count[0] += 1
            raise EOFError("fail")

        with pytest.raises(EOFError):
            c._retry("ebi", _op)

        assert call_count[0] == 1

    def test_max_retries_zero_raises_value_error(self):
        """max_retries=0 raises ValueError before any FTP call."""
        with pytest.raises(ValueError, match="max_retries must be at least 1"):
            make_client(max_retries=0)
