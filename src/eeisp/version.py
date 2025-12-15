from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version


def get_version() -> str:
    try:
        return version("eeisp")
    except PackageNotFoundError:
        return "0+unknown"


