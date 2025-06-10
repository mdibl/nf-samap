#!/root/miniconda/bin/python
"""
Author : Ryan Sonderman
Date   : 2025-06-06
Purpose: Run SAMAP with command-line input
"""

import argparse
import json
from pathlib import Path
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple


class Args(NamedTuple):
    """Command-line arguments"""

    config: Path


# --------------------------------------------------
def get_args() -> Args:
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Run SAMAP with command-line input",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.json"),
        help="Path to JSON with species names and h5ad paths",
    )

    args = parser.parse_args()

    return Args(args.config)


# --------------------------------------------------
def main() -> None:
    """Run SAMap with command-line input"""
    args = get_args()

    with open(args.config, "r") as f:
        config = json.load(f)

    maps_dir = config.get("maps", "data/maps/")
    data_dict = config.get("species", {})

    sm = SAMAP(
        sams=data_dict,
        f_maps=maps_dir,
        save_processed=True,  # If False, do not save the processed results to `*_pr.h5ad`
    )

    sm.run()
    save_samap(sm, "samap_obj")


# --------------------------------------------------
if __name__ == "__main__":
    main()
