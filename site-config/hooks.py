"""Non-destructive MkDocs entry point and compatibility hook."""

from __future__ import annotations

import os
from pathlib import Path


# The repository itself is the documentation source, so a normal clean build
# would remove files from the generated site directory.  Project policy avoids
# bulk deletion.  Marking generated files as old makes a dirty build refresh
# every source page, including the complete search index, without deleting.
_OLD_TIMESTAMP = 946684800  # 2000-01-01 00:00:00 UTC

_MEDICAL_ETHICS_ASSIGNMENT_ORDER = {
    "第一次": 1,
    "第二次": 2,
    "第三次": 3,
    "第四次": 4,
    "第五次": 5,
    "第六次": 6,
    "第七次": 7,
}


def on_pre_build(*, config, **kwargs) -> None:
    site_dir = Path(config.site_dir)
    if not site_dir.exists():
        return

    for path in site_dir.rglob("*"):
        if path.is_file():
            os.utime(path, (_OLD_TIMESTAMP, _OLD_TIMESTAMP))


def _assignment_order(page) -> tuple[int, str]:
    """Return the intended order for the medical ethics online assignments."""
    source_path = getattr(getattr(page, "file", None), "src_uri", "")
    filename = Path(source_path).name

    for prefix, order in _MEDICAL_ETHICS_ASSIGNMENT_ORDER.items():
        if filename.startswith(prefix):
            return order, filename

    return len(_MEDICAL_ETHICS_ASSIGNMENT_ORDER) + 1, filename


def on_nav(nav, *, config, files, **kwargs):
    """Keep the medical ethics assignment pages in chronological order."""
    for course_section in nav.items:
        if course_section.title != "医学伦理学":
            continue

        for subsection in getattr(course_section, "children", []):
            if subsection.title == "课后作业答案":
                subsection.children.sort(key=_assignment_order)
                return nav

    return nav


def _preserve_directory(path: str) -> None:
    """Replace MkDocs' clean step with a non-destructive directory check."""
    Path(path).mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    # MkDocs normally clears site_dir before a full build.  The repository's
    # project policy prohibits bulk deletion, so the local scripts enter the
    # CLI through this small wrapper and preserve existing generated files.
    from mkdocs import utils

    utils.clean_directory = _preserve_directory

    from mkdocs.__main__ import cli

    cli(prog_name="mkdocs")
