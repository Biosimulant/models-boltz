#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
DESKTOP_SRC_TAURI = REPO_ROOT / "bsim-platform" / "biosimulant-desktop" / "src-tauri"
DEFAULT_DESKTOP_BIN = DESKTOP_SRC_TAURI / "target" / "debug" / "biosimulant-desktop"
MODEL_DIR = REPO_ROOT / "models" / "models-boltz" / "models" / "boltz-boltz2-affinity-predictor"
EXAMPLE_CONFIG = REPO_ROOT / "models" / "models-boltz" / "examples" / "boltz2-minimal" / "config.yaml"
DEFAULT_OUTPUT = Path.cwd() / "Boltz2_Remote_GPU_Example.bsispace"
HUB_API_BASE = "https://prod-api.biosimulant.com/api"


def parse_example_inputs(config_path: Path) -> tuple[str, str]:
    protein_sequence: str | None = None
    ligand_smiles: str | None = None
    in_inputs = False
    for raw_line in config_path.read_text().splitlines():
        line = raw_line.rstrip()
        if not line or line.lstrip().startswith("#"):
            continue
        if line.startswith("  inputs:") or line == "  inputs:":
            in_inputs = True
            continue
        if in_inputs and not raw_line.startswith("    "):
            in_inputs = False
        if not in_inputs:
            continue
        stripped = line.strip()
        if stripped.startswith("protein_sequence:"):
            protein_sequence = stripped.split(":", 1)[1].strip()
        elif stripped.startswith("ligand_smiles:"):
            ligand_smiles = stripped.split(":", 1)[1].strip()
    if not protein_sequence or not ligand_smiles:
        raise RuntimeError(f"Could not parse protein_sequence and ligand_smiles from {config_path}")
    return protein_sequence, ligand_smiles


def run_cli(
    desktop_bin: Path,
    data_dir: Path | None,
    command: str,
    payload: dict[str, Any],
) -> Any:
    cmd = [str(desktop_bin), "__biosimulant_cli__", "raw", command, "--json", "--input", json.dumps(payload)]
    if data_dir is not None:
        cmd.extend(["--data-dir", str(data_dir)])
    completed = subprocess.run(
        cmd,
        cwd=DESKTOP_SRC_TAURI,
        check=True,
        text=True,
        capture_output=True,
        env={**os.environ, **({"BIOSIMULANT_HUB_ACCESS_TOKEN": os.environ["BIOSIMULANT_HUB_ACCESS_TOKEN"]} if os.environ.get("BIOSIMULANT_HUB_ACCESS_TOKEN") else {})},
    )
    envelope = json.loads(completed.stdout)
    if not envelope.get("ok", False):
        raise RuntimeError(f"{command} failed: {envelope}")
    return envelope["data"]


def load_hub_access_token() -> str:
    completed = subprocess.run(
        ["security", "find-generic-password", "-s", "com.biosimulant.desktop", "-a", "hub_access_token", "-w"],
        check=True,
        text=True,
        capture_output=True,
    )
    token = completed.stdout.strip()
    if not token:
        raise RuntimeError("Could not load a Hub access token from the desktop keychain.")
    return token


def hub_request(
    method: str,
    path: str,
    token: str,
    *,
    json_body: dict[str, Any] | None = None,
) -> Any:
    url = f"{HUB_API_BASE}{path}"
    command = [
        "curl",
        "-sS",
        "-f",
        "-X",
        method,
        url,
        "-H",
        f"Authorization: Bearer {token}",
    ]
    if json_body is not None:
        command.extend(["-H", "Content-Type: application/json", "-d", json.dumps(json_body)])
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            f"Hub request returned non-JSON for {method} {path}: {completed.stdout[:500]}"
        ) from error


def hub_upload_package(package_path: Path, token: str, visibility: str = "private") -> dict[str, Any]:
    command = [
        "curl",
        "-sS",
        "-X",
        "POST",
        f"{HUB_API_BASE}/packages/upload",
        "-H",
        f"Authorization: Bearer {token}",
        "-F",
        f"package_file=@{package_path}",
        "-F",
        f"visibility={visibility}",
    ]
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    return json.loads(completed.stdout)


def hub_download_artifact(run_id: str, artifact_id: str, token: str, output_path: Path) -> None:
    command = [
        "curl",
        "-sS",
        "-L",
        f"{HUB_API_BASE}/runs/{run_id}/artifacts/{artifact_id}",
        "-H",
        f"Authorization: Bearer {token}",
        "-o",
        str(output_path),
    ]
    subprocess.run(command, check=True, text=True, capture_output=True)


def choose_remote_size(catalog: dict[str, Any], *, require_gpu: bool) -> dict[str, Any]:
    sizes = catalog.get("sizes") or []
    credit_balance = float(catalog.get("credit_balance") or 0)
    affordable_sizes = [
        size for size in sizes if size.get("is_active") and float(size.get("credit_cost") or 0) <= credit_balance
    ]
    active_gpu_sizes = [
        size
        for size in affordable_sizes
        if size.get("is_active") and size.get("gpu_type") and (size.get("gpu_count") or 1) > 0
    ]
    if active_gpu_sizes:
        active_gpu_sizes.sort(
            key=lambda size: (
                0 if str(size.get("gpu_type", "")).upper() == "A10" else 1,
                size.get("credit_cost_per_minute") or 0,
                size.get("credit_cost") or 0,
            )
        )
        return active_gpu_sizes[0]

    if require_gpu:
        raise RuntimeError("No affordable GPU remote sizes are available for this account.")

    affordable_cpu_sizes = [size for size in affordable_sizes if not size.get("gpu_type")]
    if not affordable_cpu_sizes:
        raise RuntimeError("No affordable remote sizes are available for this account.")
    affordable_cpu_sizes.sort(
        key=lambda size: (
            -(size.get("memory_mb") or 0),
            -(size.get("timeout_seconds") or 0),
            -(size.get("cpu_cores") or 0),
            size.get("credit_cost") or 0,
        )
    )
    return affordable_cpu_sizes[0]


def extract_structure_artifact(results_payload: dict[str, Any]) -> tuple[str, str | None]:
    visuals = results_payload.get("visuals") or []
    for module in visuals:
        if not isinstance(module, dict):
            continue
        for visual in module.get("visuals") or []:
            if not isinstance(visual, dict) or visual.get("render") != "structure3d":
                continue
            data = visual.get("data") or {}
            source = data.get("source") or {}
            artifact_id = source.get("artifact_id")
            if isinstance(artifact_id, str) and artifact_id.strip():
                file_name_hint = None
                source_path = source.get("path")
                if isinstance(source_path, str) and source_path.strip():
                    file_name_hint = Path(source_path).name
                return artifact_id, file_name_hint
    raise RuntimeError("Remote run results did not include a structure3d artifact.")


def assert_results_shape(results_payload: dict[str, Any]) -> None:
    outputs = results_payload.get("outputs")
    visuals = results_payload.get("visuals")
    if not isinstance(outputs, dict):
        raise RuntimeError("Run results are missing outputs.")
    for key in ("affinity_summary", "confidence_summary", "structure_artifacts"):
        if key not in outputs:
            raise RuntimeError(f"Run outputs are missing `{key}`.")
    if not isinstance(visuals, list) or not visuals:
        raise RuntimeError("Run results are missing visuals.")
    extract_structure_artifact(results_payload)


def build_manifest(space_id: str, imported_model: dict[str, Any], protein_sequence: str, ligand_smiles: str) -> dict[str, Any]:
    return {
        "schema_version": "2.0",
        "title": "Boltz-2 Remote Example",
        "description": "Portable single-model Boltz-2 lab with a real protein-ligand affinity example intended for remote GPU execution and inline 3D structure visualization.",
        "models": [
            {
                "alias": "boltz2_affinity",
                "path": imported_model["owned_path"],
                "provenance": {
                    "owner_space_id": space_id,
                    "owned_path": imported_model["owned_path"],
                    "imported_at": imported_model.get("imported_at"),
                    "local_revision": imported_model.get("local_revision"),
                    "dirty": False,
                },
                "parameters": {
                    "runtime_mode": "managed",
                    "use_msa_server": True,
                    "accelerator": "gpu",
                    "output_format": "mmcif",
                    "override": True,
                    "recycling_steps": 3,
                    "sampling_steps": 200,
                    "diffusion_samples": 1,
                    "default_protein_sequence": protein_sequence,
                    "default_ligand_smiles": ligand_smiles,
                },
            }
        ],
        "wiring": [],
        "runtime": {
            "duration": 0.01,
            "tick_dt": 0.01,
            "initial_inputs": {},
        },
        "scientific_context": {
            "question": "Can Boltz-2 predict a protein-ligand complex remotely and return a renderable 3D structure inside the desktop app?",
            "mode": "native-boltz",
            "assumptions": [
                "Remote runs should use GPU-backed execution.",
                "MSAs are generated through the configured MSA server for this example.",
            ],
            "expected_observables": [
                "Affinity summary values from Boltz.",
                "Confidence summary values from Boltz.",
                "A renderable structure3d visual backed by a persisted artifact.",
            ],
        },
    }


def build_layout() -> dict[str, Any]:
    return {
        "nodes": [
            {"id": "boltz2_affinity", "type": "model", "x": 180, "y": 120},
        ],
        "viewport": {"x": 0, "y": 0, "zoom": 1},
    }


def patch_owned_model_runtime_dependency(
    data_dir: Path,
    space_id: str,
    owned_path: str,
    package_spec: str,
) -> None:
    model_manifest_path = data_dir / "spaces" / space_id / Path(owned_path) / "model.yaml"
    if not model_manifest_path.is_file():
        raise RuntimeError(f"Owned model manifest is missing: {model_manifest_path}")
    contents = model_manifest_path.read_text()
    original = "- boltz[cuda]==2.0.2"
    replacement = f"- {package_spec}"
    if original in contents:
        contents = contents.replace(original, replacement)
    elif "- boltz==2.0.2" in contents:
        contents = contents.replace("- boltz==2.0.2", replacement)
    else:
        raise RuntimeError(f"Could not locate Boltz runtime dependency in {model_manifest_path}")
    model_manifest_path.write_text(contents)


def mirror_and_wait_for_remote_results(
    desktop_bin: Path,
    data_dir: Path,
    imported_space_id: str,
    stage_result: dict[str, Any],
    remote_size: dict[str, Any],
    timeout_seconds: int,
    poll_seconds: int,
) -> tuple[str, str, dict[str, Any], dict[str, Any]]:
    remote_run = run_cli(
        desktop_bin,
        data_dir,
        "hub_create_remote_run",
        {
            "payload": {
                "space_id": stage_result["hub_space_id"],
                "space_commit": stage_result["space_commit"],
                "simulation_config": {"duration": 0.01, "tick_dt": 0.01, "initial_inputs": {}},
                "remote_size_id": remote_size["id"],
            }
        },
    )
    remote_run_id = str(remote_run.get("id") or "").strip()
    if not remote_run_id:
        raise RuntimeError(f"Remote run creation did not return an id: {remote_run}")
    remote_status = str(remote_run.get("status") or "queued")

    local_run = run_cli(
        desktop_bin,
        data_dir,
        "create_run",
        {
            "space_id": imported_space_id,
            "status": remote_status,
            "execution_target": "remote",
            "hub_run_id": remote_run_id,
            "simulation_config": {"duration": 0.01, "tick_dt": 0.01, "initial_inputs": {}},
        },
    )
    local_run_id = local_run["id"]

    deadline = time.time() + timeout_seconds
    last_remote_state: dict[str, Any] | None = None
    while time.time() < deadline:
        current = run_cli(
            desktop_bin,
            data_dir,
            "hub_get_remote_run",
            {"run_id": remote_run_id},
        )
        last_remote_state = current
        status = str(current.get("status") or "")
        if status in {"completed", "failed", "cancelled"}:
            results_payload = run_cli(
                desktop_bin,
                data_dir,
                "hub_get_remote_run_results",
                {"run_id": remote_run_id},
            )
            run_cli(
                desktop_bin,
                data_dir,
                "sync_remote_run",
                {
                    "id": local_run_id,
                    "hub_run_id": remote_run_id,
                    "status": status,
                    "error_message": current.get("error_message"),
                    "duration_seconds": current.get("duration_seconds"),
                    "results_payload": results_payload,
                },
            )
            return local_run_id, remote_run_id, current, results_payload
        time.sleep(poll_seconds)
    raise TimeoutError(f"Remote run {remote_run_id} did not finish within {timeout_seconds} seconds. Last state: {last_remote_state}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build and validate a portable Boltz-2 .bsispace package.")
    parser.add_argument("--desktop-bin", type=Path, default=DEFAULT_DESKTOP_BIN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-remote-validation", action="store_true")
    parser.add_argument("--timeout-seconds", type=int, default=3600)
    parser.add_argument("--poll-seconds", type=int, default=15)
    parser.add_argument("--allow-cpu-fallback", action="store_true")
    args = parser.parse_args()

    protein_sequence, ligand_smiles = parse_example_inputs(EXAMPLE_CONFIG)
    build_root = Path(tempfile.mkdtemp(prefix="boltz2-bsispace-build-"))
    export_dir = build_root / "export"
    export_dir.mkdir(parents=True, exist_ok=True)
    package_tmp_path = export_dir / "Boltz2_Remote_GPU_Example.bsispace"

    data_dir = build_root / "data"
    import_data_dir = build_root / "imported-data"
    data_dir.mkdir(parents=True, exist_ok=True)
    import_data_dir.mkdir(parents=True, exist_ok=True)

    try:
        created_space = run_cli(
            args.desktop_bin,
            data_dir,
            "create_space",
            {
                "title": "Boltz-2 Remote Example",
                "description": "Portable remote GPU Boltz-2 example lab",
            },
        )
        space_id = created_space["id"]
        imported_model = run_cli(
            args.desktop_bin,
            data_dir,
            "import_model_into_space_from_path",
            {
                "spaceId": space_id,
                "path": str(MODEL_DIR),
                "alias": "boltz2_affinity",
            },
        )
        patch_owned_model_runtime_dependency(
            data_dir,
            space_id,
            imported_model["owned_path"],
            "boltz[cuda]==2.0.2",
        )
        run_cli(
            args.desktop_bin,
            data_dir,
            "save_space",
            {
                "id": space_id,
                "manifest": build_manifest(space_id, imported_model, protein_sequence, ligand_smiles),
                "wiringLayout": build_layout(),
            },
        )
        export_result = run_cli(
            args.desktop_bin,
            data_dir,
            "export_space_package",
            {"space_id": space_id, "output_path": str(package_tmp_path)},
        )
        if Path(export_result["path"]) != package_tmp_path:
            package_tmp_path = Path(export_result["path"])
        run_cli(
            args.desktop_bin,
            data_dir,
            "preview_package",
            {"package_path": str(package_tmp_path)},
        )
        imported_package = run_cli(
            args.desktop_bin,
            import_data_dir,
            "import_package",
            {"package_path": str(package_tmp_path)},
        )
        imported_space_id = (
            imported_package.get("local_space_id")
            or imported_package.get("space_id")
            or imported_package.get("local_id")
        )
        if not imported_space_id:
            raise RuntimeError(f"Could not determine imported space id from {imported_package}")
        imported_space = run_cli(
            args.desktop_bin,
            import_data_dir,
            "get_space",
            {"id": imported_space_id},
        )
        if len(imported_space.get("manifest", {}).get("models", [])) != 1:
            raise RuntimeError("Imported portable space does not contain exactly one model node.")

        validation_summary: dict[str, Any] | None = None
        if not args.skip_remote_validation:
            catalog = run_cli(
                args.desktop_bin,
                import_data_dir,
                "hub_get_remote_catalog",
                {},
            )
            gpu_size = choose_remote_size(catalog, require_gpu=not args.allow_cpu_fallback)
            stage_result = run_cli(
                args.desktop_bin,
                import_data_dir,
                "hub_stage_remote_space",
                {"space_id": imported_space_id},
            )
            local_run_id, remote_run_id, remote_state, results_payload = mirror_and_wait_for_remote_results(
                args.desktop_bin,
                import_data_dir,
                imported_space_id,
                stage_result,
                gpu_size,
                args.timeout_seconds,
                args.poll_seconds,
            )
            if remote_state.get("status") != "completed":
                raise RuntimeError(f"Remote run did not complete successfully: {remote_state}")
            assert_results_shape(results_payload)
            artifact_id, file_name_hint = extract_structure_artifact(results_payload)
            cached_artifact = run_cli(
                args.desktop_bin,
                import_data_dir,
                "hub_cache_remote_run_artifact",
                {
                    "run_id": remote_run_id,
                    "artifact_id": artifact_id,
                    "file_name_hint": file_name_hint,
                },
            )
            cached_path = Path(cached_artifact["local_path"])
            if not cached_path.is_file():
                raise RuntimeError(f"Desktop artifact cache did not produce a file at {cached_path}")
            local_results = run_cli(
                args.desktop_bin,
                import_data_dir,
                "get_run_results",
                {"run_id": local_run_id},
            )
            assert_results_shape(local_results)
            validation_summary = {
                "remote_size_id": gpu_size["id"],
                "remote_size_name": gpu_size["display_name"],
                "remote_run_id": remote_run_id,
                "local_run_id": local_run_id,
                "cached_artifact_path": str(cached_path),
                "status": remote_state["status"],
            }

        args.output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(package_tmp_path, args.output)
        summary = {
            "package_path": str(args.output),
            "temp_package_path": str(package_tmp_path),
            "validated_remote_run": validation_summary,
            "imported_space_id": imported_space_id,
        }
        print(json.dumps(summary, indent=2))
    finally:
        # Keep build_root for post-run inspection only when the caller exports `KEEP_BOLTZ2_BSISPACE_BUILD=1`.
        if not bool(int(__import__("os").environ.get("KEEP_BOLTZ2_BSISPACE_BUILD", "0"))):
            shutil.rmtree(build_root, ignore_errors=True)


if __name__ == "__main__":
    main()
