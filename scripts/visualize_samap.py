#!/usr/bin/env python3
"""
Author : Ryan Sonderman
Date   : 2025-06-09
Version: 1.0.0
Purpose: Visualize SAMAP results from pickle
"""

import argparse
import os
import pickle
import json
import csv
from typing import NamedTuple, Optional
from pathlib import Path
from samap.mapping import SAMAP
from samap.analysis import get_mapping_scores, sankey_plot, chord_plot
import matplotlib.pyplot as plt
import holoviews as hv

hv.extension("bokeh")


class Args(NamedTuple):
    """
    Command-line arguments for the script.

    Attributes:
        input: Path to the input SAMAP pickle file.
        output_dir: Optional directory where visualizations will be saved.
        sample_sheet: Path to the sample sheet CSV file for annotations.
    """
    input: str
    output_dir: Optional[str]
    sample_sheet: Path


# --------------------------------------------------
def get_args() -> Args:
    """
    Parse and return command-line arguments.

    Returns:
        Args: A NamedTuple containing the input file, output directory, and sample sheet path.
    """
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

    parser.add_argument(
        "--sample-sheet",
        type=Path,
        required=True,
        help="Path to the sample sheet CSV",
    )

    args = parser.parse_args()

    return Args(args.input, args.output_dir, args.sample_sheet)


# --------------------------------------------------
def create_output_dir(output_dir: str) -> None:
    """
    Create the specified output directory if it does not already exist.

    Args:
        output_dir (str): Path to the output directory.
    """
    """Create output directory if it does not exist"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Output directory created: {output_dir}")
    else:
        print(f"Output directory already exists: {output_dir}")


# --------------------------------------------------
def load_samap_pickle(pickle_file: str) -> SAMAP:
    """
    Load the SAMAP object from a pickle file.

    Args:
        pickle_file (str): Path to the SAMAP pickle file.

    Returns:
        SAMAP: Loaded SAMAP object.
    """
    with open(pickle_file, "rb") as f:
        samap_obj = pickle.load(f)
    print(f"Loaded SAMAP object of type: {type(samap_obj)}")
    return samap_obj


# --------------------------------------------------
def save_sankey_plot(mapping_table, 
                    output_dir: str,
                    align_thr=0.05,
                    file_ext='html',
                    file_name='sankey',
                    toolbar=True,
                    title='Sankey Plot',):
    """
    Generate and save a Sankey plot using Holoviews.

    Args:
        mapping_table (pandas.dataFrame): Pairwise mapping scores.
        output_dir (str): Directory where the plot will be saved. 
        align_thr (float, default=0.05): The alignment score threshold below which to remove cell type mappings.
        file_ext (str, default='html'): File extension applied to the plot.
        file_name (str, default='sankey'): Name of the saved file.
        toolbar (bool, default=True): Whether to include toolbar in the plot.
        title (str, default='Sankey Plot'): Custom title for exported HTML file.
    """
    file_fmt = file_ext.lstrip('.')
    sankey_obj = sankey_plot(mapping_table, align_thr=align_thr)
    sankey_outfile = os.path.join(output_dir, f"{file_name}.{file_fmt}")
    try:
        print(f"[INFO] Saving sankey plot to '{sankey_outfile}'")
        hv.save(sankey_obj, sankey_outfile, backend="bokeh", fmt=file_fmt, toolbar=toolbar, title=title)
        print(f"[INFO] Successfully saved sankey plot to '{sankey_outfile}'")
    except Exception as e:
        print(f"[ERROR] Could not save sankey plot. Error: {e}")


# --------------------------------------------------
def save_chord_plot(mapping_table, 
                    output_dir: str, 
                    align_thr=0.05, 
                    file_ext='html',
                    file_name='chord',
                    toolbar=True,
                    title='Chord Plot',):
    """
    Generate and save a Chord plot using Holoviews.

    Args:
        mapping_table (pandas.dataFrame): Pairwise mapping scores.
        output_dir (str): Directory where the plot will be saved.
        align_thr (float, default=0.05): The alignment score threshold below which to remove cell type mappings.
        file_ext (str, default='html'): File extension applied to the plot.
        file_name (str, default='chord'): Name of the saved file.
        toolbar (bool, default=True): Whether to include toolbar in the plot.
    """
    file_fmt = file_ext.lstrip('.')
    chord_obj = chord_plot(mapping_table, align_thr)
    chord_outfile = os.path.join(output_dir, f"{file_name}.{file_fmt}")
    try:
        print(f"[INFO] Saving chord plot to '{chord_outfile}'")
        hv.save(chord_obj, chord_outfile, backend="bokeh", fmt=file_fmt, toolbar=toolbar, title=title)
        print(f"[INFO] Successfully saved chord plot to '{chord_outfile}'")
    except Exception as e:
        print(f"[ERROR] Could not save chord plot. Error: {e}")


# --------------------------------------------------
def save_mapping_scores(samap: SAMAP, keys: dict, output_dir: str):
    """
    Save the highest mapping scores and pairwise mapping scores to CSV files,
    and generate and save a Sankey plot.

    Args:
        samap (SAMAP): The SAMAP object containing the results.
        keys (dict): Dictionary of annotations for each sample.
        output_dir (str): Directory where the results will be saved.

    Returns:
        tuple: Paths to the saved CSV files and Sankey plot HTML.
    """
    highest_mapping_scores, pairwise_mapping_scores = get_mapping_scores(
        samap, keys, n_top=0
    )

    # Save mapping dataframes to CSV
    hms_outfile = os.path.join(output_dir, "hms.csv")
    highest_mapping_scores.to_csv(hms_outfile)
    print(f"Saved highest mapping scores to: {hms_outfile}")

    pms_outfile = os.path.join(output_dir, "pms.csv")
    pairwise_mapping_scores.to_csv(pms_outfile)
    print(f"Saved pairwise mapping scores to: {pms_outfile}")

    # Save Sankey plot
    sankey_html_outfile = save_sankey_plot(
        pairwise_mapping_scores, output_dir
    )
    
    # Save chord plot
    save_chord_plot(pairwise_mapping_scores, output_dir=output_dir)

    return hms_outfile, pms_outfile, sankey_html_outfile


# --------------------------------------------------
def save_scatter_plot(samap: SAMAP, output_dir: str):
    """
    Save a scatter plot of the SAMAP results to the output directory.

    Args:
        samap (SAMAP): The SAMAP object containing the results.
        output_dir (str): Directory where the scatter plot will be saved.

    Returns:
        str: Path to the saved scatter plot image.
    """
    samap.scatter()
    sc_outfile = os.path.join(output_dir, "scatter.png")
    plt.savefig(sc_outfile, dpi=300)
    print(f"Saved scatter plot to {sc_outfile}")
    plt.close()

    return sc_outfile


# --------------------------------------------------
def save_log(
    input_file: str, outputs: dict, log_file: str, keys_json: dict
) -> None:
    """
    Save a JSON log with input, keys, and outputs in Nextflow module format.

    Args:
        input_file (str): Path to the input SAMAP pickle file.
        outputs (dict): Dictionary of paths to the generated output files.
        log_file (str): Path to the log file.
        keys_json (dict): Dictionary of sample annotations.

    Returns:
        None
    """
    log = {
        "input": input_file,
        "keys": keys_json,
        "outputs": outputs,
    }
    with open(log_file, "w") as f:
        json.dump(log, f, indent=2)
    print(f"Saved log to {log_file}")


# --------------------------------------------------
def load_keys_from_sample_sheet(sample_sheet_path: Path) -> dict:
    """
    Load a dictionary of keys from the sample sheet CSV. The keys are 
    mapped from id2 to their corresponding annotations.

    Args:
        sample_sheet_path (Path): Path to the sample sheet CSV.

    Returns:
        dict: Dictionary where the key is id2 and the value is the annotation.
    """
    keys = {}
    with open(sample_sheet_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            keys[row["id2"]] = row["annotation"]
    return keys


# --------------------------------------------------
def main() -> None:
    """
    Visualize SAMAP results from a pickle file.

    This function:
    1. Loads the SAMAP object from the provided pickle file.
    2. Loads sample annotations from the sample sheet.
    3. Saves mapping scores, Sankey plot, chord plot, and scatter plot to the output directory.
    4. Saves a log of the process in JSON format.
    """

    args = get_args()

    # Create output directory if it does not exist
    create_output_dir(args.output_dir)
    # Load SAMAP object from pickle
    sm = load_samap_pickle(args.input)

    # Load keys from sample sheet
    keys = load_keys_from_sample_sheet(args.sample_sheet)
    print(keys)  # or use as needed in your visualization logic

    outputs = {}

    # Save mapping scores and Sankey
    hms_file, pms_file, sankey_file = save_mapping_scores(
        sm, keys, args.output_dir
    )
    outputs["highest_mapping_scores"] = hms_file
    outputs["pairwise_mapping_scores"] = pms_file
    outputs["sankey_html"] = sankey_file

    # Save scatter plot
    scatter_file = save_scatter_plot(sm, args.output_dir)
    outputs["scatter_plot"] = scatter_file

    # Save log, now including keys json
    log_file = os.path.join(args.output_dir, "vis.log")
    save_log(args.input, outputs, log_file, keys)


# --------------------------------------------------
if __name__ == "__main__":
    main()
