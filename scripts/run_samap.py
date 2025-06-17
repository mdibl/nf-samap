#!/root/miniconda/bin/python
"""
Author : Ryan Sonderman
Date   : 2025-06-06
Purpose: Run SAMAP post-processing with a pickled SAMAP object
"""

import argparse
import pickle
from pathlib import Path
from samap.mapping import SAMAP
from samap.utils import save_samap


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Run SAMAP post-processing with a pickled SAMAP object",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the pickled SAMAP object (.pkl)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        type=Path,
        default=Path("."),
        help="Directory to save the SAMAP pickle output",
    )
    return parser.parse_args()


# --------------------------------------------------
def main():
    """Run SAMap with command-line input"""
    args = get_args()
    print(f"[INFO] Loading SAMAP object from: {args.input}")
    with open(args.input, "rb") as f:
        samap = pickle.load(f)

    print("[INFO] Running SAMAP...")
    samap.run()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_file = args.output_dir / "samap_results.pkl"
    save_samap(samap, str(output_file))
    print(f"[INFO] SAMAP results saved to {output_file}")


# --------------------------------------------------
if __name__ == "__main__":
    main()
