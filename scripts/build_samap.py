#!/usr/bin/env python3
"""
Author : Ryan Sonderman, Markus Sujansky
Date   : 2025-06-16
Version: 1.0.0
Purpose: Build a SAMAP object from sams and maps
"""

import argparse
import csv
import pickle
import os
import ast
from log_utils import log
from samap.mapping import SAMAP
from samap.utils import save_samap
from typing import NamedTuple
from pathlib import Path
from itertools import permutations


class Args(NamedTuple):
    """ Command-line arguments for the script"""
    
    sams_dir: Path   # List of SAM pickle files
    id2: str     # List of id2 strings
    mappings: Path #List of Mapping Files
    maps: Path          # Path to the maps directory
    name: str           # Name of the output pickle
    output_dir: Path    # Path to the output directory
    bit_threshold: float # Optional BLAST bitscore floor for map filtering (None = no filtering)



# --------------------------------------------------
def get_args() -> Args:
    """
    Parse and return command-line arguments.

    Returns:
        Args: A named tuple containing parsed command-line arguments for sams_dir, sample_sheet, and f_maps.
    """
    parser = argparse.ArgumentParser(
        description='Build a SAMAP object from a directory of SAMs and a sample sheet',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-d', '--sams-dir',
        required=True,
        type=Path,
        nargs='+',
        help='Directory containing SAM pickle files'
    )

    parser.add_argument(
        '-i', '--id2',
        required=True,
        type=str,
        nargs='+',
        help='list of id2 from the Sample Sheet'
    )

    parser.add_argument(
        '-p', '--mappings',
        required=False,
        type=str,
        nargs='+',
        help='list of mappings from the Sample Sheet'
    )

    parser.add_argument(
        '-m', '--maps',
        required=True,
        type=Path,
        help='Path to the maps directory'
    )
    
    parser.add_argument(
        '-n', '--name',
        required=False,
        type=str,
        help='Name of the output pickle',
        default='samap.pkl'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=False,
        type=Path,
        help='Path to the output directory',
        default=Path('.')
    )

    parser.add_argument(
        '-b', '--bit-threshold',
        required=False,
        type=float,
        default=None,
        help='If provided, filter all BLAST maps to hits with bitscore >= this before building SAMAP'
    )

    args = parser.parse_args()
    return Args(args.sams_dir, args.id2, args.mappings, args.maps, args.name, args.output_dir, args.bit_threshold)


# --------------------------------------------------
def load_species_dict(id2: str, sams_dir: Path) -> dict:
    """
    Load a dictionary of species, mapping id2 to corresponding SAM objects from the sams_dir directory.

    Args:
        id2 (list): List of ids in the sample sheet (FUTURE: change from id2 -> id)
        sams_dir (Path): Path to the directory containing the SAM pickle files.

    Returns:
        dict: A dictionary with id2 as the key and the corresponding SAM object as the value.
    """
    species = {}
    for val in id2:
        matching_files = [f for f in sams_dir if f.name.startswith(val) and f.suffix == ".pkl"]
        if not matching_files:
            log(f"  No SAM pickle found for '{val}' in provided files", "ERROR")
            continue
        sam_path = matching_files[0]
        with open(sam_path, "rb") as f:
            species[val] = pickle.load(f)
        log(f"  Loaded SAM for '{val}' from '{sam_path}'", "INFO")
    return species


    # --------------------------------------------------

def load_mapping_dict(id2: str, mapping_dir: list) -> dict:
    """
    Load a dictionary of mappings connecting protein/transcript ids in the BLAST maps to the format in the inputted h5ad files, creating a dictionary with the species id as the key

    Args:
        id2 (list): List of ids in the sample sheet (FUTURE: change from id2 -> id)
        mapping_dir (Path): Path to the directory containing the mapping files.

    Returns:
        dict: A dictionary with id2 as the key and the corresponding SAM object as the value.
    """
    mapping_dict = {}
    mapping_dir = [Path(p) for p in mapping_dir]
    
    for val in id2:
        # Find the mapping file whose name contains the species ID
        matching = [p for p in mapping_dir if val in p.name]
        if not matching:
            log(f"  No mapping file found for species '{val}'", "ERROR")
            continue
        map_path = matching[0]
        
        try:
            mapping_dict[val] = []
            with open(map_path, "r") as f:
                content = f.read().strip().strip("[]")
                pairs = content.split("), (")
                for p in pairs:
                    p = p.replace("(", "").replace(")", "")
                    fasta, gene = p.split(",")
                    mapping_dict[val].append((fasta.strip(), gene.strip()))
        except Exception as e:
            log(f"  Failed to parse mapping file '{map_path}' for species '{val}': {e}", "ERROR")
    
    return mapping_dict


    # --------------------------------------------------


def filter_map_file(map_file: Path, bit_threshold: float,
                src_root: Path, dst_root: Path) -> Path:

    """
    Filter one BLAST map (outfmt 6, tab-separated, no header) to rows whose
    bitscore (column 11) is >= bit_threshold. Writes to dst_root preserving the
    file's path relative to src_root. Returns the output path.
    """

    rel = map_file.relative_to(src_root)
    out_path = dst_root / rel
    out_path.parent.mkdir(parents=True, exist_ok=True)

    kept = total = 0
    with open(map_file) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if not line.strip():
                continue
            total += 1
            fields = line.rstrip('\n').split('\t')
            try:
                bitscore = float(fields[11])
            except (IndexError, ValueError):
                log(f"    Malformed line in {map_file.name}, skipped: {line!r}", "WARN")
                continue
            if bitscore >= bit_threshold:
                fout.write(line)
                kept += 1
    log(f"    {map_file.name}: kept {kept}/{total} hits (bitscore >= {bit_threshold})", "INFO")
    return out_path

# --------------------------------------------------
def main() -> None:
    """
    Main entry point for the script.

    This function:
    1. Parses command-line arguments.
    2. Loads the species dictionary from the sample sheet and SAM files.
    3. Validates the maps directory.
    4. Optionally can filter the BLAST maps inside the maps directory on a bit score threshold
    5. Creates a SAMAP object using the loaded species and maps data.
    6. Saves the SAMAP object to a pickle file.
    """

    # Parse command-line arguments
    log("Loading arguments", "INFO")
    args = get_args()
    sams_dir = args.sams_dir
    log(f"  Using SAMs directory '{sams_dir}'", "DEBUG")
    
    if args.mappings:
        log(f"  Mappings argument provided: {args.mappings}", "DEBUG")
        mapping_dir = args.mappings
    else:
        log("  No mappings argument provided", "DEBUG")
        mapping_dir = None

    maps = str(Path(args.maps).resolve())
    if not maps.endswith('/'): # SAMap *will* crash if passed a dir without a '/'
        maps += '/'
    log(f"  Using maps directory '{maps}'", "DEBUG")
    id2 = args.id2
    log(f"  Using id2 list '{id2}'", "DEBUG")
    name = args.name
    log(f"  SAMAP object will be saved with name '{name}'", "DEBUG")
    output_dir = args.output_dir
    log(f"  SAMAP object will be saved to '{output_dir}'", "DEBUG")
    
    # Load species dictionary from sample sheet
    log("Loading species dictionary from sample sheet", "INFO")
    species_dict = load_species_dict(id2, sams_dir)
    log(f"Loaded species dictionary with {len(species_dict)} entries", "INFO")

    # Ensure each SAM object knows its species
    for key, sam in species_dict.items():
        sam.species = key.strip()
        sam.species_id = key.strip()
        if hasattr(sam, 'adata') and hasattr(sam.adata, 'uns'):
            sam.adata.uns['species'] = key.strip()
    log("Set internal species identifiers for all SAM objects", "INFO")

    # Ensure maps is valid and formatted correctly
    log(f"Ensuring validity of '{maps}'", "INFO")
    if not maps.endswith('/'):
        maps += '/'
        log(f"Provided maps directory does not end with '/', changing to '{maps}'", "WARN")
    if not Path(maps).exists():
        error_message = f"Maps directory '{maps}' does not exist"
        log(error_message, "ERROR")
        raise FileNotFoundError(error_message)
    log(f"Maps directory found at '{maps}'", "INFO")

    bit_threshold = args.bit_threshold
    if bit_threshold is not None:
        src_root = Path(maps)
        dst_root = src_root.parent / f"{src_root.name}_bitfiltered"
        log(f"Filtering BLAST maps to bitscore >= {bit_threshold} -> '{dst_root}/'", "INFO")
        map_files = sorted(src_root.rglob('*.txt'))
        log(f"  Found {len(map_files)} map file(s) to filter", "INFO")
        if not map_files:
            log(f"  No .txt map files found under '{maps}'", "WARN")
        for map_file in map_files:
            filter_map_file(map_file, bit_threshold, src_root, dst_root)
        maps = str(dst_root) + '/'   # <-- SAMAP must read the FILTERED dir
        log(f"Maps directory for SAMAP set to filtered copy '{maps}'", "INFO")
    else:
        log("No bit-score threshold provided; using unfiltered maps", "INFO")




    if mapping_dir is not None:
        valid_mappings = [m for m in mapping_dir if m and m.lower() not in ('null', 'none', 'na', '')]

        #Load mapping dict from sample sheet
        if valid_mappings:
            log("Loading mapping dictionary from sample sheet to line up BLAST protein/transcript headers with SAM feature type", "INFO")
            mapping_dict = load_mapping_dict(id2, valid_mappings)
            log(f"Loaded mapping dictionary with {len(mapping_dict)} entries", "INFO")

            samap = SAMAP(sams=species_dict, f_maps=maps, save_processed=False, names=mapping_dict)
        else:
            log("No valid mapping dictionaries provided, skipping", "INFO")
            samap = SAMAP(sams=species_dict, f_maps=maps, save_processed=False)
    else:
        log("No mappings argument provided, skipping", "INFO")
        samap = SAMAP(sams=species_dict, f_maps=maps, save_processed=False)
    
    log(f"Successfully created SAMAP object with {len(samap.sams)} SAMs", "INFO")
        
    # Save SAMAP object
    log("Attempting to pickle SAMAP object", "INFO")
    save_samap(samap, os.path.join(output_dir, name))
    log(f"Successfully pickled SAMAP object '{name}' to '{output_dir}'")

# --------------------------------------------------
if __name__ == '__main__':
    main()
