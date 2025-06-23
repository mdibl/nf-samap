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
import csv
from log_utils import log
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
        log(f"Successfully created output directory at '{output_dir}'", "INFO")
    else:
        log(f"Output directory already exists at '{output_dir}, skipping step", "INFO")


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
    log(f"Successfully loaded pickle of type '{type(samap_obj)}", "INFO")
    return samap_obj


# --------------------------------------------------
def save_mapping_scores(samap: SAMAP, 
                        keys: dict, 
                        output_dir: str,
                        n_top=0,
                        hms_name='hms',
                        pms_name='pms'):
    """
    Save the highest mapping scores and pairwise mapping scores to CSV files.

    Args:
        samap (SAMAP): The SAMAP object containing the results.
        keys (dict): Dictionary of annotations for each sample.
        output_dir (str): Directory where the results will be saved.
        n_top (int, default=0): Average the alignment scores for the n top cells (0 averages all cells)
        hms_name (str, default='hms'): Name which the highest mapping scores table will be saved to.
        pms_name (str, default='pms'): Name which the pairwise mapping scores table will be saved to.

    Returns:
        tuple (pandas.dataFrame): Highest mapping scores and pairwise mapping scores.
    """
    hms, pms = get_mapping_scores(sm=samap, keys=keys, n_top=n_top)
    hms_outfile = os.path.join(output_dir, f"{hms_name}.csv")
    try: # Save the highest mapping scores to csv
        log(f"  Attempting to save highest mapping scores to '{hms_outfile}'", "INFO")
        hms.to_csv(hms_outfile)
        log(f"  Successfully saved highest mapping scores to '{hms_outfile}'", "INFO")
    except Exception as e:
        log(f"  Failed to save highest mapping scores to '{hms_outfile}'. Error: {e}", "ERROR")
    pms_outfile = os.path.join(output_dir, f"{pms_name}.csv")
    try: # Save the pairwise mapping scores to csv
        log(f"  Attempting to save pairwise mapping scores to '{pms_outfile}'", "INFO")
        pms.to_csv(pms_outfile)
        log(f"  Scussessfully save pairwise mapping scores to '{pms_outfile}'", "INFO")
    except Exception as e:
        print(f"[ERROR] Could not save pairwise mapping scores. Error: {e}")
        log(f"  Failed to save pairwise mapping scores to '{pms_outfile}'. Error: {e}", "ERROR")
    return hms, pms # Return the data frames


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
        log(f"  Attempting to save sankey plot to '{sankey_outfile}'", "INFO")
        hv.save(sankey_obj, sankey_outfile, backend="bokeh", fmt=file_fmt, toolbar=toolbar, title=title)
        log(f"  Successfully saved sankey plot to '{sankey_outfile}'", "INFO")
    except Exception as e:
        log(f"  Failed to save sankey plot to '{sankey_outfile}'. Error: {e}", "ERROR")


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
        log(f"  Attempting to save chord plot to '{chord_outfile}'", "INFO")
        hv.save(chord_obj, chord_outfile, backend="bokeh", fmt=file_fmt, toolbar=toolbar, title=title)
        log(f"  Successfully saved chord plot to '{chord_outfile}'", "INFO")
    except Exception as e:
        log(f"  Failed to save chord plot to '{chord_outfile}'. Error: {e}")


# --------------------------------------------------
def save_scatter_plot(samap: SAMAP, output_dir: str, out_name='scatter', dpi=300):
    """
    Save a scatter plot of the SAMAP results to the output directory.

    Args:
        samap (SAMAP): The SAMAP object containing the results.
        output_dir (str): Directory where the scatter plot will be saved.
        out_name (str, default='scatter'): Name of the saved file.
        dpi (int, default=300): DPI of the saved image.

    Returns:
        str: Path to the saved scatter plot image.
    """
    samap.scatter()
    scatter_outfile = os.path.join(output_dir, f"{out_name}.png")
    log(f"  Attempting to save scatter plot to '{scatter_outfile}'", "INFO")
    plt.savefig(scatter_outfile, dpi=dpi)
    log(f"  Successfully saved scatter plot to '{scatter_outfile}'", "INFO")
    plt.close()


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
            log(f"  Loading annotation for {row['id2']}: {row['annotation']}", "INFO")
    log("Successfully loaded all annotation keys")
    return keys


# --------------------------------------------------
def main() -> None:
    """
    Visualize SAMAP results from a pickle file.

    This function:
    1. Loads the SAMAP object from the provided pickle file.
    2. Loads sample annotations from the sample sheet.
    3. Saves mapping scores, Sankey plot, chord plot, and scatter plot to the output directory.
    """

    log("Beginning execution of script", "INFO")
    args = get_args()
    
    log(f"Attempting to create output directory at '{args.output_dir}'", "INFO")
    create_output_dir(args.output_dir)
    log(f"Attempting to load SAMAP pickle from '{args.input}'", "INFO")
    sm = load_samap_pickle(args.input)
    log(f"Attempting to load annotation keys from '{args.sample_sheet}'", "INFO")
    keys = load_keys_from_sample_sheet(args.sample_sheet)
    
    # Save and get mapping scores and use them to generate plots
    log("Attempting to save mapping scores")
    _, pms = save_mapping_scores(sm, keys, args.output_dir)
    log("Attempting to create chord plot", "INFO")
    save_chord_plot(mapping_table=pms, output_dir=args.output_dir)
    log("Attempting to create sankey plot", "INFO")
    save_sankey_plot(mapping_table=pms, output_dir=args.output_dir)
    log("Attempting to create scatter plot", "INFO")
    save_scatter_plot(samap=sm, output_dir=args.output_dir)
    log("Script complete", "INFO")


# --------------------------------------------------
if __name__ == "__main__":
    main()
