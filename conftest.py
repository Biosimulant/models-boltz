from __future__ import annotations

import sys
from pathlib import Path

import pytest


def _model_dir_for_test(path: Path) -> Path | None:
    for parent in path.parents:
        if (parent / "model.yaml").exists() or (parent / "model.yml").exists():
            return parent
    return None


@pytest.fixture(autouse=True)
def _isolate_copied_src_packages(request):
    """Keep copied `src.*` model packages from leaking across lab tests."""

    test_path = Path(str(request.fspath)).resolve()
    model_dir = _model_dir_for_test(test_path)
    if model_dir is None:
        return

    for name in [key for key in sys.modules if key == "src" or key.startswith("src.")]:
        sys.modules.pop(name, None)

    model_dir_text = str(model_dir)
    sys.path[:] = [entry for entry in sys.path if entry != model_dir_text]
    sys.path.insert(0, model_dir_text)
