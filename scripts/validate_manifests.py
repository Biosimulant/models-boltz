#!/usr/bin/env python3
"""Validate model.yaml and space.yaml manifests in models."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import yaml


ROOT = Path(__file__).resolve().parents[1]


def _load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError("manifest root must be a mapping")
    return data


def _local_repo_aliases() -> set[str]:
    repo_name = ROOT.name
    return {
        repo_name,
        f"Biosimulant/{repo_name}",
    }


def _validate_model_manifest(path: Path) -> list[str]:
    errors: list[str] = []
    try:
        manifest = _load_yaml(path)
    except Exception as exc:  # noqa: BLE001
        return [f"{path}: {exc}"]

    standard = manifest.get("standard")
    if not isinstance(standard, str) or not standard.strip():
        errors.append(f"{path}: missing required 'standard'")

    biosim = manifest.get("biosim") if isinstance(manifest.get("biosim"), dict) else {}
    entrypoint = biosim.get("entrypoint") or manifest.get("entrypoint")
    if not isinstance(entrypoint, str) or not entrypoint.strip():
        errors.append(f"{path}: missing required 'biosim.entrypoint'")

    runtime = manifest.get("runtime") if isinstance(manifest.get("runtime"), dict) else {}
    deps = runtime.get("dependencies") if isinstance(runtime.get("dependencies"), dict) else {}
    packages = deps.get("packages") if isinstance(deps.get("packages"), list) else []
    for spec in packages:
        if not isinstance(spec, str) or "==" not in spec:
            errors.append(f"{path}: dependency package '{spec}' must be pinned with '=='")

    for key in ("requirements_file", "lockfile"):
        value = deps.get(key)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            errors.append(f"{path}: runtime.dependencies.{key} must be a non-empty string")

    return errors


def _validate_space_manifest(path: Path) -> list[str]:
    errors: list[str] = []
    try:
        manifest = _load_yaml(path)
    except Exception as exc:  # noqa: BLE001
        return [f"{path}: {exc}"]

    models = manifest.get("models")
    if not isinstance(models, list) or not models:
        return [f"{path}: missing non-empty 'models' list"]

    local_aliases = _local_repo_aliases()

    for idx, entry in enumerate(models):
        if not isinstance(entry, dict):
            errors.append(f"{path}: models[{idx}] must be a mapping")
            continue
        alias = entry.get("alias")
        if not isinstance(alias, str) or not alias.strip():
            errors.append(f"{path}: models[{idx}] missing required alias")
        repo = entry.get("repo") or entry.get("repo_full_name")
        if not isinstance(repo, str) or not repo.strip():
            errors.append(f"{path}: models[{idx}] missing repo/repo_full_name")
            continue

        manifest_path = entry.get("manifest_path")
        if not isinstance(manifest_path, str) or not manifest_path.strip():
            errors.append(f"{path}: models[{idx}] missing required manifest_path")
            continue

        # Only enforce local path existence for refs targeting this repo.
        if repo not in local_aliases:
            continue

        target = (ROOT / manifest_path).resolve()
        if ROOT not in target.parents and target != ROOT:
            errors.append(f"{path}: models[{idx}].manifest_path must be repo-relative: {manifest_path}")
        elif not target.exists():
            errors.append(f"{path}: models[{idx}].manifest_path does not exist: {manifest_path}")

    return errors


def main() -> int:
    errors: list[str] = []

    model_manifests = sorted(ROOT.rglob("model.yaml"))
    space_manifests = sorted(ROOT.rglob("space.yaml"))

    for path in model_manifests:
        if "/templates/" in str(path):
            continue
        errors.extend(_validate_model_manifest(path))
    for path in space_manifests:
        if "/templates/" in str(path):
            continue
        errors.extend(_validate_space_manifest(path))

    if errors:
        print("Manifest validation failed:")
        for err in errors:
            print(f" - {err}")
        return 1

    print(f"Validated {len(model_manifests)} model manifest(s) and {len(space_manifests)} space manifest(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
