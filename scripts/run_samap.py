#!/root/miniconda/bin/python
"""
Author : Ryan Sonderman
Date   : 2025-06-06
Purpose: Run SAMAP with command-line input
"""

import argparse
import json
import datetime
from pathlib import Path
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple


class Args(NamedTuple):
    """Command-line arguments"""

    config: Path
    output_dir: Path


# --------------------------------------------------
def get_args() -> Args:
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Run SAMAP with command-line input",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        default=Path("config.json"),
        help="Path to JSON with species names and h5ad paths",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        type=Path,
        default=Path("results/samap_objects"),
        help="Directory to save the SAMAP pickle output",
    )

    args = parser.parse_args()

    return Args(args.config, args.output_dir)


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
        save_processed=False,  # If False, do not save the processed results to `*_pr.h5ad`
    )

    sm.run()

    save_results(sm, args.output_dir)


# --------------------------------------------------
def save_results(samap: SAMAP, output_dir: Path) -> None:
    """Save SAMAP results to the specified directory"""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a unique filename with timestamp
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = output_dir / f"samap_results_{timestamp}.pkl"
    save_samap(samap, str(output_file))
    print(f"SAMAP results saved to {output_file}")


# --------------------------------------------------
if __name__ == "__main__":
    main()
