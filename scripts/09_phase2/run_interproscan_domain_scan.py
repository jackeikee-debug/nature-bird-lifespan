"""Run InterProScan for a local FASTA, using WSL automatically on Windows."""

from __future__ import annotations

import argparse
import os
import pathlib
import platform
import subprocess


def to_wsl_path(path: pathlib.Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    rest = resolved.as_posix().split(":", 1)[1]
    return f"/mnt/{drive}{rest}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--applications", required=True)
    parser.add_argument("--wsl-distro", default="Ubuntu-24.04")
    parser.add_argument("--cpu", type=int, default=2)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("", ".xml", ".gff3", ".json", ".tsv"):
        candidate = pathlib.Path(str(args.output) + suffix)
        if candidate.exists() and candidate != args.output:
            candidate.unlink()
    if args.output.exists():
        args.output.unlink()

    if platform.system().lower().startswith("win"):
        input_path = to_wsl_path(args.input)
        output_path = to_wsl_path(args.output)
        cwd = to_wsl_path(pathlib.Path.cwd())
        command = [
            "wsl",
            "-d",
            args.wsl_distro,
            "--",
            "bash",
            "-lc",
            (
                f"cd {cwd!r} && interproscan.sh "
                f"-i {input_path!r} -f tsv -dp -appl {args.applications!r} "
                f"-cpu {args.cpu} -o {output_path!r}"
            ),
        ]
    else:
        command = [
            "interproscan.sh",
            "-i",
            str(args.input),
            "-f",
            "tsv",
            "-dp",
            "-appl",
            args.applications,
            "-cpu",
            str(args.cpu),
            "-o",
            str(args.output),
        ]

    env = os.environ.copy()
    subprocess.run(command, check=True, env=env)


if __name__ == "__main__":
    main()
