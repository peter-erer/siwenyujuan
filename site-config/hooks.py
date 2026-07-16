"""Non-destructive MkDocs entry point and compatibility hook."""

from __future__ import annotations

import os
from pathlib import Path


# The repository itself is the documentation source, so a normal clean build
# would remove files from the generated site directory.  Project policy avoids
# bulk deletion.  Marking generated files as old makes a dirty build refresh
# every source page, including the complete search index, without deleting.
_OLD_TIMESTAMP = 946684800  # 2000-01-01 00:00:00 UTC


def on_pre_build(*, config, **kwargs) -> None:
    site_dir = Path(config.site_dir)
    if not site_dir.exists():
        return

    for path in site_dir.rglob("*"):
        if path.is_file():
            os.utime(path, (_OLD_TIMESTAMP, _OLD_TIMESTAMP))


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
