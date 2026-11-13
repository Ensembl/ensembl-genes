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
"""Format the final gene confidence metrics as a two-column CSV."""
import csv
from pathlib import Path


def write_gene_confidence_summary(
    metrics: dict[str, int | float], output: Path
) -> None:
    """Write ordered metric/value rows, with two decimal places for percentages."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        for metric, value in metrics.items():
            writer.writerow(
                [metric, f"{value:.2f}" if metric.endswith("_pct") else value]
            )
