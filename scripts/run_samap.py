#!/root/miniconda/bin/python
"""
Author : Ryan Sonderman
Date   : 2025-06-06
Purpose: Run SAMAP with command-line input
"""

import argparse
import datetime
import csv
from pathlib import Path
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple


class Args(NamedTuple):
    """Command-line arguments"""

    sample_sheet: Path
    maps: Path
    output_dir: Path


# --------------------------------------------------
def get_args() -> Args:
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Run SAMAP with command-line input",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--sample-sheet",
        dest="sample_sheet",
        type=Path,
        required=True,
        help="Path to the sample sheet CSV file",
    )

    parser.add_argument(
        "-m",
        "--maps",
        type=Path,
        default=Path("data/maps/"),
        help="Path to the maps directory",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        type=Path,
        default=Path("."),
        help="Directory to save the SAMAP pickle output",
    )

    args = parser.parse_args()

    return Args(args.sample_sheet, args.maps, args.output_dir)


# --------------------------------------------------
def load_species_dict(sample_sheet_path: Path) -> dict:
    """Load species dictionary from sample sheet CSV, sorted by key"""
    species = {}
    with open(sample_sheet_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            species[row["id2"]] = row["h5ad"]
    # Sort the dictionary by key (id2)
    return dict(sorted(species.items()))


# --------------------------------------------------
def main() -> None:
    """Run SAMap with command-line input"""
    args = get_args()

    species_dict = load_species_dict(args.sample_sheet)
    print(species_dict)
    print(args.maps)

    sm = SAMAP(
        sams=species_dict,
        f_maps=str(args.maps),
        save_processed=False,
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
