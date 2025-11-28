"""Generate API reference pages for MkDocs using mkdocs-gen-files.

This script discovers Python modules under `src/python/ensembl` and
creates a `reference/` page for each module with a mkdocstrings directive.

It is intentionally simple and avoids importing project modules.
"""
from pathlib import Path
import mkdocs_gen_files

ROOT = Path("src/python")
PKG = Path("src/python/ensembl")
OUT_DIR = "reference"


def iter_modules(root: Path):
    for path in sorted(root.rglob("*.py")):
        # skip tests or private modules if needed
        if path.name == "__init__.py":
            continue
        rel = path.relative_to(ROOT).with_suffix("")
        parts = rel.parts
        # only include modules under the top-level 'ensembl' package
        if len(parts) == 0 or parts[0] != "ensembl":
            continue
        module = ".".join(parts)
        yield module


def build():
    nav = []
    for module in iter_modules(PKG):
        page_path = Path(OUT_DIR) / (module + ".md")
        # create parent dirs
        mkdocs_gen_files.set_edit_path(str(Path("src") / page_path))
        with mkdocs_gen_files.open(page_path, "w") as fh:
            fh.write(f"::: {module}\n")
        nav.append(page_path.as_posix())

    # write a simple summary index
    summary = Path(OUT_DIR) / "summary.md"
    with mkdocs_gen_files.open(summary, "w") as fh:
        fh.write("# API reference\n\n")
        for p in nav:
            fh.write(f"- [{p}]({p})\n")


if __name__ == "__main__":
    build()
