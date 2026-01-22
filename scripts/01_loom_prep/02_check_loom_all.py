#!/usr/bin/env python3
# Requires: Python 3.x, loompy (pip install loompy). If NumPy ≥ 2, this script auto-handles np.string_.
"""
Tiny loom sanity checker (version / global attributes) -> to check loom was correctly produced

Usage:
  python check_loom.py --version /path/to/your.loom
  python check_loom.py --attrs   /path/to/your.loom
"""

import argparse
import sys
import loompy
import numpy as np

# Handle NumPy API changes (np.string_ → np.bytes_ in NumPy 2)
if not hasattr(np, "string_"):
    np.string_ = np.bytes_

def check_version(path: str):
    try:
        with loompy.connect(path, mode="r") as ds:
            ver = ds.attrs.get("LOOM_SPEC_VERSION", "unknown")
            print(f"LOOM_SPEC_VERSION: {ver}")
            # A common expectation for pySCENIC-compatible looms is "3.0.0"
            if ver != "3.0.0":
                print("NOTE: expected 3.0.0; please verify downstream tools.")
    except Exception as e:
        print(f"[ERROR] reading {path}: {e}", file=sys.stderr)
        sys.exit(1)

def list_attrs(path: str):
    try:
        with loompy.connect(path, mode="r") as ds:
            print("Global attributes:")
            for k in ds.attrs.keys():
                print(f"  {k} = {ds.attrs[k]}")
    except Exception as e:
        print(f"[ERROR] reading {path}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--version", action="store_true", help="print LOOM_SPEC_VERSION")
    g.add_argument("--attrs",   action="store_true", help="list global attributes")
    ap.add_argument("loom", help="path to loom file")
    args = ap.parse_args()

    if args.version:
        check_version(args.loom)
    elif args.attrs:
        list_attrs(args.loom)
