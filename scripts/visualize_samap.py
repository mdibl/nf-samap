#!/usr/bin/env python3
"""
Author : Ryan Sonderman
Date   : 2025-06-09
Purpose: Visualize SAMAP results from pickle
"""

import argparse
import os
import pickle
import datetime
import json
from typing import NamedTuple, Optional
from samap.mapping import SAMAP
from samap.analysis import get_mapping_scores, sankey_plot
import matplotlib.pyplot as plt
import holoviews as hv

hv.extension("bokeh")


class Args(NamedTuple):
    input: str
    output_dir: Optional[str]


# --------------------------------------------------
def get_args() -> Args:
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Visualize SAMAP results from pickle",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i", "--input", metavar="PKL", required=True, help="Input SAMAP pickle file"
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        metavar="DIR",
        default=".",
        help="Optional output directory for visualizations",
    )

    args = parser.parse_args()

    return Args(args.input, args.output_dir)


# --------------------------------------------------
def create_output_dir(output_dir: str) -> None:
    """Create output directory if it does not exist"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Output directory created: {output_dir}")
    else:
        print(f"Output directory already exists: {output_dir}")


# --------------------------------------------------
def load_samap_pickle(pickle_file: str) -> SAMAP:
    """Load SAMAP object from a pickle file"""
    with open(pickle_file, "rb") as f:
        samap_obj = pickle.load(f)
    print(f"Loaded SAMAP object of type: {type(samap_obj)}")
    return samap_obj


# --------------------------------------------------
def save_sankey_plot(pairwise_mapping_scores, output_dir: str, timestamp: str) -> str:
    """Generate and save Sankey plot as HTML using holoviews"""
    sankey_obj = sankey_plot(
        pairwise_mapping_scores, align_thr=0.05, species_order=["pl", "hy", "sc"]
    )
    sankey_html_outfile = os.path.join(output_dir, f"sankey_{timestamp}.html")
    try:
        hv.save(sankey_obj, sankey_html_outfile, backend="bokeh")
        print(f"Saved Sankey plot as interactive HTML to: {sankey_html_outfile}")
    except Exception as e:
        print(f"Warning: Could not save Sankey plot as HTML. Error: {e}")
    return sankey_html_outfile


# --------------------------------------------------
def save_mapping_scores(samap: SAMAP, keys: dict, output_dir: str, timestamp: str):
    """Save mapping scores to CSV files and Sankey plot"""
    highest_mapping_scores, pairwise_mapping_scores = get_mapping_scores(
        samap, keys, n_top=0
    )

    # Save mapping dataframes to CSV
    hms_outfile = os.path.join(output_dir, f"hms_{timestamp}.csv")
    highest_mapping_scores.to_csv(hms_outfile)
    print(f"Saved highest mapping scores to: {hms_outfile}")

    pms_outfile = os.path.join(output_dir, f"pms_{timestamp}.csv")
    pairwise_mapping_scores.to_csv(pms_outfile)
    print(f"Saved pairwise mapping scores to: {pms_outfile}")

    # Save Sankey plot using the new function
    sankey_html_outfile = save_sankey_plot(pairwise_mapping_scores, output_dir, timestamp)

    return hms_outfile, pms_outfile, sankey_html_outfile


# --------------------------------------------------
def save_scatter_plot(samap: SAMAP, output_dir: str, timestamp: str):
    """Save scatter plot of SAMAP results"""
    samap.scatter()
    sc_outfile = os.path.join(output_dir, f"scatter_{timestamp}.png")
    plt.savefig(sc_outfile, dpi=300)
    print(f"Saved scatter plot to {sc_outfile}")
    plt.close()

    return sc_outfile


# --------------------------------------------------
def save_log(input_file: str, outputs: dict, log_file: str, timestamp: str) -> None:
    """Save a JSON log with input, outputs, and timestamp in nf module format"""
    log = {
        "input": input_file,
        "outputs": outputs,
        "timestamp": timestamp
    }
    with open(log_file, "w") as f:
        json.dump(log, f, indent=2)
    print(f"Saved log to {log_file}")


# --------------------------------------------------
def main() -> None:
    """Visualize SAMAP results from pickle"""

    args = get_args()

    # Create output directory if it does not exist
    create_output_dir(args.output_dir)
    # Load SAMAP object from pickle
    sm = load_samap_pickle(args.input)

    # Create a global timestamp for all outputs
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    # save mapping scores (currently hardcoded keys)
    keys = {"pl": "cluster", "hy": "Cluster", "sc": "tissue"}
    
    outputs = {}

    # Save mapping scores and Sankey
    hms_file, pms_file, sankey_file = save_mapping_scores(sm, keys, args.output_dir, timestamp)
    outputs["highest_mapping_scores"] = hms_file
    outputs["pairwise_mapping_scores"] = pms_file
    outputs["sankey_html"] = sankey_file

    # Save scatter plot
    scatter_file = save_scatter_plot(sm, args.output_dir, timestamp)
    outputs["scatter_plot"] = scatter_file

    # Save log 
    log_file = os.path.join(args.output_dir, f"visualization_{timestamp}.log")
    save_log(args.input, outputs, log_file, timestamp)


# --------------------------------------------------
if __name__ == "__main__":
    main()
