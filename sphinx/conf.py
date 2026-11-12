# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations

import os
import sys
from importlib import metadata
from pathlib import Path
from shutil import rmtree

from sphinx.application import Sphinx


ROOT_DIR = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT_DIR / "src" / "python"
PACKAGE_DIR = SOURCE_DIR / "ensembl" / "genes"

sys.path.insert(0, str(SOURCE_DIR))

project = "Ensembl Genes"
author = "Ensembl"
copyright = "2016-2026, EMBL-European Bioinformatics Institute"

try:
    release = metadata.version("ensembl-genes")
except metadata.PackageNotFoundError:
    release = "0.1.0"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
root_doc = "index"
exclude_patterns = [
    "_build",
    "api/generated",
    "Thumbs.db",
    ".DS_Store",
]

myst_enable_extensions = ["colon_fence"]
myst_heading_anchors = 3

html_theme = "furo"
html_title = "Ensembl Genes"
html_logo = str(ROOT_DIR / "mkdocs" / "img" / "ebang.png")

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
autodoc_typehints = "description"
autodoc_mock_imports = [
    "Bio",
    "matplotlib",
    "numpy",
    "pandas",
    "pyfaidx",
    "pymysql",
    "pyranges",
    "pyranges1",
    "requests",
    "seaborn",
    "sqlalchemy",
    "yaml",
    "xmltodict",
]

LITERAL_ONLY_MODULES = {
    "ensembl.genes.metadata.beta_patcher": (
        "This module depends on the optional ensembl-metadata-api package at "
        "import time, so Sphinx shows the source instead of importing it."
    ),
    "ensembl.genes.tracking.bioproject_tracking": (
        "This module reads a local JSON configuration file at import time, so "
        "Sphinx shows the source instead of importing it."
    ),
}
SOURCE_FILE_SUFFIXES = {
    ".conf",
    ".csv",
    ".ipynb",
    ".json",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}
LITERAL_LANGUAGES = {
    ".conf": "text",
    ".csv": "text",
    ".ipynb": "json",
    ".json": "json",
    ".toml": "toml",
    ".txt": "text",
    ".yaml": "yaml",
    ".yml": "yaml",
}


def _relative_doc_path(source: Path, start: Path) -> str:
    return Path(os.path.relpath(source, start=start)).as_posix()


def _module_name(path: Path) -> str:
    return ".".join(path.relative_to(SOURCE_DIR).with_suffix("").parts)


def _package_name(path: Path) -> str:
    return ".".join(path.relative_to(SOURCE_DIR).parts)


def _markdown_doc_name(path: Path) -> str:
    return f"{_package_name(path.parent)}.{path.stem}"


def _source_file_doc_name(path: Path) -> str:
    safe_name = path.name.replace("-", "_").replace(".", "_")
    return f"{_package_name(path.parent)}.{safe_name}"


def _title(text: str, underline: str = "=") -> str:
    return f"{text}\n{underline * len(text)}\n\n"


def _write_module_page(output_dir: Path, module_path: Path) -> None:
    module = _module_name(module_path)
    output_path = output_dir / f"{module}.rst"
    source_path = _relative_doc_path(module_path, output_dir)

    lines = [
        _title(module),
        f"Source: ``{module_path.relative_to(ROOT_DIR).as_posix()}``\n\n",
    ]

    if module in LITERAL_ONLY_MODULES:
        lines.extend(
            [
                f".. note:: {LITERAL_ONLY_MODULES[module]}\n\n",
                f".. literalinclude:: {source_path}\n",
                "   :language: python\n",
                "   :linenos:\n",
            ]
        )
    else:
        lines.extend(
            [
                f".. automodule:: {module}\n",
                "   :members:\n",
                "   :show-inheritance:\n",
                "   :undoc-members:\n",
            ]
        )

    output_path.write_text("".join(lines))


def _write_markdown_page(output_dir: Path, markdown_path: Path) -> None:
    output_path = output_dir / f"{_markdown_doc_name(markdown_path)}.md"
    source_path = _relative_doc_path(markdown_path, output_dir)
    output_path.write_text(
        "\n".join(
            [
                f"```{{include}} {source_path}",
                ":parser: myst_parser.sphinx_",
                "```",
                "",
            ]
        )
    )


def _write_source_file_page(output_dir: Path, source_file_path: Path) -> None:
    doc_name = _source_file_doc_name(source_file_path)
    output_path = output_dir / f"{doc_name}.rst"
    source_path = _relative_doc_path(source_file_path, output_dir)
    language = LITERAL_LANGUAGES.get(source_file_path.suffix, "text")

    output_path.write_text(
        "".join(
            [
                _title(source_file_path.name),
                f"Source: ``{source_file_path.relative_to(ROOT_DIR).as_posix()}``\n\n",
                f".. literalinclude:: {source_path}\n",
                f"   :language: {language}\n",
                "   :linenos:\n",
            ]
        )
    )


def _write_package_page(
    output_dir: Path,
    package_dir: Path,
    package_dirs: set[Path],
    module_paths: set[Path],
    markdown_paths: set[Path],
    source_file_paths: set[Path],
) -> None:
    package = _package_name(package_dir)
    output_path = output_dir / f"{package}.rst"
    child_dirs = sorted(path for path in package_dirs if path.parent == package_dir)
    child_modules = sorted(path for path in module_paths if path.parent == package_dir)
    child_markdown = sorted(
        (path for path in markdown_paths if path.parent == package_dir),
        key=lambda path: (path.name != "README.md", path.name.lower()),
    )
    child_source_files = sorted(
        path for path in source_file_paths if path.parent == package_dir
    )

    toctree_entries = []
    toctree_entries.extend(
        f"{path.stem} <{_markdown_doc_name(path)}>" for path in child_markdown
    )
    toctree_entries.extend(
        f"{path.name} <{_package_name(path)}>" for path in child_dirs
    )
    toctree_entries.extend(
        f"{path.stem} <{_module_name(path)}>" for path in child_modules
    )
    toctree_entries.extend(
        f"{path.name} <{_source_file_doc_name(path)}>" for path in child_source_files
    )

    lines = [
        _title(package),
        f"Source directory: ``{package_dir.relative_to(ROOT_DIR).as_posix()}``\n\n",
    ]

    if (package_dir / "__init__.py").exists():
        lines.extend(
            [
                "Package Contents\n",
                "----------------\n\n",
                f".. automodule:: {package}\n",
                "   :members:\n",
                "   :show-inheritance:\n",
                "   :undoc-members:\n\n",
            ]
        )

    if toctree_entries:
        lines.extend(
            [
                "Folders, Scripts, and Files\n",
                "---------------------------\n\n",
                ".. toctree::\n",
                "   :maxdepth: 2\n\n",
                *[f"   {entry}\n" for entry in toctree_entries],
            ]
        )

    output_path.write_text("".join(lines))


def _write_reference_index(output_dir: Path) -> None:
    output_path = output_dir / "index.rst"
    output_path.write_text(
        "".join(
            [
                _title("Generated Code Reference"),
                "This tree is generated from ``src/python/ensembl/genes`` at "
                "Sphinx build time. It includes package folders, README files, "
                "Python scripts, and source/config files.\n\n",
                ".. toctree::\n",
                "   :maxdepth: 4\n\n",
                "   ensembl.genes\n",
            ]
        )
    )


def _generate_reference_pages(app: Sphinx) -> None:
    """Generate a package tree that includes folders, READMEs, and scripts."""
    output_dir = Path(app.srcdir) / "reference" / "generated"
    if output_dir.exists():
        rmtree(output_dir)
    output_dir.mkdir(parents=True)

    package_dirs = {
        path
        for path in [PACKAGE_DIR, *PACKAGE_DIR.rglob("*")]
        if path.is_dir()
        and "__pycache__" not in path.parts
        and not any(part.startswith(".") for part in path.parts)
    }
    module_paths = {
        path
        for path in PACKAGE_DIR.rglob("*.py")
        if path.name != "__init__.py" and "__pycache__" not in path.parts
    }
    markdown_paths = {
        path for path in PACKAGE_DIR.rglob("*.md") if "__pycache__" not in path.parts
    }
    source_file_paths = {
        path
        for path in PACKAGE_DIR.rglob("*")
        if path.is_file()
        and path.suffix in SOURCE_FILE_SUFFIXES
        and "__pycache__" not in path.parts
    }

    _write_reference_index(output_dir)

    for markdown_path in sorted(markdown_paths):
        _write_markdown_page(output_dir, markdown_path)

    for source_file_path in sorted(source_file_paths):
        _write_source_file_page(output_dir, source_file_path)

    for module_path in sorted(module_paths):
        _write_module_page(output_dir, module_path)

    for package_dir in sorted(package_dirs):
        _write_package_page(
            output_dir,
            package_dir,
            package_dirs,
            module_paths,
            markdown_paths,
            source_file_paths,
        )


def setup(app: Sphinx) -> dict[str, bool]:
    app.connect("builder-inited", _generate_reference_pages)
    return {
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
