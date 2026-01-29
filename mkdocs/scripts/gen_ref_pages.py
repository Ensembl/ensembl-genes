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

"""
Generate API reference pages under reference/ using mkdocs-gen-files.
"""

from pathlib import Path
import mkdocs_gen_files

SRC_ROOT = Path("src/python")
PKG_ROOT = SRC_ROOT / "ensembl/genes"
OUT_DIR = Path("reference")

nav = mkdocs_gen_files.Nav()

# ----------------------------
# 1. Add Home section
# ----------------------------
home_pages = [
    ("Overview", "index.md"),
    ("Install", "install.md"),
    ("Usage", "usage.md"),
]
for title, path in home_pages:
    nav["Home", title] = path

# ----------------------------
# 2. Add Development section
# ----------------------------
dev_pages = [
    ("Code of Conduct", "code_of_conduct.md"),
    ("Coverage report", "coverage.md"),
]
for title, path in dev_pages:
    nav["Development", title] = path

# ----------------------------
# 3. Add Reference section
# ----------------------------

# Find all python modules under src/python/ensembl
for path in sorted(PKG_ROOT.rglob("*.py")):
    if path.name == "__init__.py":
        continue

    # Example:
    #   src/python/ensembl/genes/automation/pre_release_ftp.py
    #
    # → module_parts = ("ensembl", "genes", "automation", "pre_release_ftp")
    module_parts = path.relative_to(SRC_ROOT).with_suffix("").parts
    module = ".".join(module_parts)

    # strip the top-level "ensembl" so literate-nav creates reference/ensembl/... correctly
    out_path = Path(*module_parts).with_suffix(".md")
    # Create nested directories
    mkdocs_gen_files.set_edit_path(out_path, path)

    with mkdocs_gen_files.open(out_path, "w") as f:
        f.write(f"# `{module}`\n\n::: {module}\n")

    # Add to navigation under "Code reference"
    nav["Code reference", *module_parts] = out_path.as_posix()

# ----------------------------
# 3a. Capture README.md files
# ----------------------------

for readme_path in sorted(PKG_ROOT.rglob("README.md")):
    # Example:
    # src/python/ensembl/genes/automation/README.md
    module_parts = readme_path.relative_to(SRC_ROOT).with_suffix("").parts
    out_path = Path(*module_parts).with_suffix(".md")

    # Copy the content of README.md to the mkdocs output
    mkdocs_gen_files.set_edit_path(out_path, readme_path)
    with mkdocs_gen_files.open(out_path, "w") as f:
        content = readme_path.read_text()
        # Optionally add a header with the folder name
        f.write(f"# {module_parts[-1].replace('_', ' ').title()}\n\n")
        f.write(content)

    # Add README to nav under "Code reference"
    nav["Code reference", *module_parts] = out_path.as_posix()
# Write summary.md for mkdocs-literate-nav
summary_path = OUT_DIR / "summary.md"
with mkdocs_gen_files.open(summary_path, "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
