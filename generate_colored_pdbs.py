#!/usr/bin/env python3
"""
Generate color-coded PDB files from epistatic antigenicity scores.

Usage:
    python3 generate_colored_pdbs.py <csv_file> <pdb_file> [pdb_file2 ...] [--attribute ATTR] [--threshold THR]

Example:
    python3 generate_colored_pdbs.py epistatic_scores.csv protein.pdb structure2.pdb --attribute significant_antigenicity_median --threshold 0.01
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import sys


def apply_scores_to_pdb(csv_file, pdb_file, attribute='significant_antigenicity_median', threshold=0.01,
                        position_cols=('pos', 'position'), make_positive=False):
    """
    Read csv_file (scores), apply the `attribute` values as B-factors into pdb_file.

    Parameters
    - csv_file: path to csv with epistatic antigenicity scores
    - pdb_file: path to input pdb
    - attribute: column name in csv to write into B-factor
    - threshold: absolute cutoff to include residues in PyMOL selection
    - position_cols: candidate column names for residue index (tries in order)
    - make_positive: if True use abs(attribute) when writing B-factors

    Returns
    - dict with keys: output_pdb, pymol_command, scores (DataFrame)
    """

    # load scores
    scores = pd.read_csv(csv_file, index_col=0).reset_index()

    # find position column
    position_col = next((c for c in position_cols if c in scores.columns), None)
    if position_col is None:
        raise ValueError("None of {} found in csv columns: {}".format(position_cols, scores.columns.tolist()))

    if attribute not in scores.columns:
        raise ValueError("Attribute '{}' not found in csv columns: {}".format(attribute, scores.columns.tolist()))

    # ensure numeric positions and numeric attribute
    scores = scores.copy()
    scores[position_col] = pd.to_numeric(scores[position_col], errors='coerce').astype('Int64')
    scores[attribute] = pd.to_numeric(scores[attribute], errors='coerce')

    # map pos -> value
    score_map = {int(r): float(v) for r, v in zip(scores[position_col], scores[attribute]) if pd.notna(r) and pd.notna(v)}

    # prepare output filename
    csv_base = Path(csv_file).stem
    pdb_base = Path(pdb_file).stem
    out_pdb = str(Path(pdb_file).parent / "{}_{}.pdb".format(pdb_base, attribute))

    # read and modify PDB
    out_lines = []
    with open(pdb_file, 'r') as fh:
        for line in fh:
            # treat ATOM and HETATM lines that contain residue number in cols 23-26 (1-based)
            if line.startswith(('ATOM', 'HETATM')):
                # safe parse residue number from columns 22:26 (0-based slice)
                resnum_str = line[22:26].strip()
                try:
                    resnum = int(resnum_str)
                except Exception:
                    resnum = None

                if resnum is not None and resnum in score_map:
                    change = score_map[resnum]
                    if make_positive:
                        change = abs(change)
                    # format B-factor field (columns 61-66, 0-based slice [60:66]) as 6.2f
                    bf = "{:6.2f}".format(change)
                    # ensure line long enough
                    if len(line) < 66:
                        line = line.rstrip('\n').ljust(66) + line[len(line):]
                    line = line[:60] + bf + line[66:]
            out_lines.append(line.rstrip('\n'))

    # write output pdb
    with open(out_pdb, 'w') as fh:
        fh.write("\n".join(out_lines) + "\n")

    # build pymol command: select residues with abs(attribute) >= threshold
    selected_pos = scores.loc[scores[attribute].abs() >= threshold, position_col].dropna().astype(int).astype(str).unique().tolist()
    
    # sort numerically for cleaner output
    selected_pos = sorted([int(p) for p in selected_pos])
    selected_pos_str = [str(p) for p in selected_pos]
    
    hashtag = '#1'
    if selected_pos_str:
        select_str = ",".join(selected_pos_str)
        command = (
            "hide atoms;hide /B cartoon;hide /C cartoon;"
            "select {}/A:{};show sel surfaces;"
            "select subtract sel;set bgColor white;lighting flat;lighting shadows true intensity 0.5;"
            "color by bfactor palette -1.04,blue:0,white:1,red;tile;".format(hashtag, select_str)
        )
    else:
        command = (
            "hide atoms;hide /B cartoon;hide /C cartoon;"
            "set bgColor white;lighting flat;lighting shadows true intensity 0.5;"
            "color by bfactor palette -1.04,blue:0,white:1,red;tile;"
        )

    return {
        "output_pdb": out_pdb,
        "pymol_command": command,
        "scores": scores,
        "pdb_file": pdb_file,
        "attribute": attribute
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate color-coded PDB files from epistatic antigenicity scores.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 generate_colored_pdbs.py scores.csv protein1.pdb protein2.pdb
  python3 generate_colored_pdbs.py scores.csv protein.pdb --attribute significant_antigenicity_median --threshold 0.05
        """
    )
    
    parser.add_argument('csv_file', help='CSV file with epistatic antigenicity scores (must have position and attribute columns)')
    parser.add_argument('pdb_files', nargs='+', help='One or more PDB files to process')
    parser.add_argument('--attribute', default='significant_antigenicity_median',
                        help='Column name in CSV to use for B-factor coloring (default: significant_antigenicity_median)')
    parser.add_argument('--threshold', type=float, default=0.01,
                        help='Absolute value threshold for residue selection in PyMOL (default: 0.01)')
    parser.add_argument('--position-col', default='position',
                        help='Column name for residue position (default: position, falls back to pos)')
    parser.add_argument('--output-commands', default='pymol_commands.txt',
                        help='Output file for PyMOL commands (default: pymol_commands.txt)')
    
    args = parser.parse_args()
    
    # Validate inputs
    csv_path = Path(args.csv_file)
    if not csv_path.exists():
        print("Error: CSV file not found: {}".format(args.csv_file), file=sys.stderr)
        sys.exit(1)
    
    pdb_paths = [Path(p) for p in args.pdb_files]
    for pdb_path in pdb_paths:
        if not pdb_path.exists():
            print("Error: PDB file not found: {}".format(pdb_path), file=sys.stderr)
            sys.exit(1)
    
    # Process each PDB
    results = []
    commands_lines = []
    
    print("Processing {} PDB file(s) with CSV: {}".format(len(pdb_paths), args.csv_file))
    print("Using attribute: {}, threshold: {}\n".format(args.attribute, args.threshold))
    
    for pdb_path in pdb_paths:
        try:
            result = apply_scores_to_pdb(
                str(csv_path),
                str(pdb_path),
                attribute=args.attribute,
                threshold=args.threshold,
                position_cols=(args.position_col, 'pos', 'position')
            )
            results.append(result)
            
            print("✓ Generated: {}".format(result['output_pdb']))
            
            # Store command with PDB name for clarity
            commands_lines.append("# PyMOL commands for {}".format(Path(result['output_pdb']).name))
            commands_lines.append(result['pymol_command'])
            commands_lines.append("")
            
        except Exception as e:
            print("✗ Error processing {}: {}".format(pdb_path, e), file=sys.stderr)
            sys.exit(1)
    
    # Write all commands to single output file
    output_commands_path = Path(args.output_commands)
    with open(output_commands_path, 'w') as fh:
        fh.write("\n".join(commands_lines))
    
    print("\n✓ PyMOL commands written to: {}".format(output_commands_path))
    print("\nSummary: Generated {} PDB file(s)".format(len(results)))
    for result in results:
        print("  • {}".format(Path(result['output_pdb']).name))


if __name__ == "__main__":
    main()
