#!/usr/bin/env python3

import os
import shutil
import subprocess as sp
import sys


SIM_ID = "smoke"
TEMPERATURE = 0.8
DURATION = 1000
FRAME_INTERVAL = 100
SEED = 1


def repo_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def upside_home():
    return os.environ.get("UPSIDE_HOME", repo_root())


def example_root():
    return os.path.abspath(os.path.dirname(__file__))


def input_dir():
    return os.path.join(example_root(), "inputs")


def output_dir():
    return os.path.join(example_root(), "outputs", SIM_ID)


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def main():
    ensure_dir(output_dir())

    config_path = os.path.join(input_dir(), "4w52_ligand.up")
    if not os.path.exists(config_path):
        sp.check_call([sys.executable, os.path.join(example_root(), "0.prepare.py")])

    run_path = os.path.join(output_dir(), "4w52_ligand.run.up")
    log_path = os.path.join(output_dir(), "4w52_ligand.run.log")
    shutil.copyfile(config_path, run_path)

    cmd = [
        os.path.join(upside_home(), "obj", "upside"),
        "--duration", str(DURATION),
        "--frame-interval", str(FRAME_INTERVAL),
        "--temperature", str(TEMPERATURE),
        "--seed", str(SEED),
        "--disable-recentering",
        run_path,
    ]
    print("Running:", " ".join(cmd))
    with open(log_path, "w") as log_handle:
        sp.check_call(cmd, stdout=log_handle, stderr=sp.STDOUT)
    print("Run written to:", run_path)


if __name__ == "__main__":
    main()
