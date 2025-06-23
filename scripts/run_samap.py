#!/root/miniconda/bin/python
"""
Author : Ryan Sonderman
Date   : 2025-06-06
Version: 1.0.0
Purpose: Run SAMAP post-processing with a pickled SAMAP object
"""

import argparse
import pickle
from log_utils import log, capture_output
from pathlib import Path
from samap.mapping import SAMAP  # noqa: F401
from samap.utils import save_samap


# --------------------------------------------------
def get_args():
    """
    Parse and return command-line arguments.

    Returns:
        argparse.Namespace: A Namespace object containing the parsed arguments.
    """
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
    parser.add_argument(
        "-n", 
        "--name",
        type=str,
        default="samap_results.pkl",
        help="Name of the saved SAMAP pickle"
    )
    return parser.parse_args()


# --------------------------------------------------
def main():
    """
    Main function to run SAMAP post-processing.

    1. Loads the pickled SAMAP object from the provided input path.
    2. Runs the SAMap post-processing.
    3. Saves the processed SAMAP object to the specified output directory.
    """
    # Get command-line arguments
    args = get_args()

    # Load the pickled SAMAP object
    log(f"Loading SAMAP object from {args.input}", "INFO")
    with open(args.input, "rb") as f:
        samap = pickle.load(f)

    # Run SAMap post-processing
    log("Attempting to run SAMap", "INFO")
    capture_output(samap.run)
    log("Successfully ran SAMap", "INFO")

    # Save the processed SAMAP object
    log("Attempting to save SAMAP object", "INFO")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_file = args.output_dir / args.name
    save_samap(samap, str(output_file))
    log(f"Successfully saved SAMAP results to {output_file}", "INFO")


# --------------------------------------------------
if __name__ == "__main__":
    main()
