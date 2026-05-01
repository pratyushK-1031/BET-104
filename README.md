# BET-104 Assignment

This workflow analyzes every Arginine residue located within an
alpha-helix triplet (HHH). For each case, it calculates the signed 3D
angle formed between neighboring C-α → centroid vectors, and then
categorizes the angle distribution based on the size group of the
preceding amino acid.

## Visualization

![density_plot](final/plot_for_arg_with_valid_runs.png)

List of PDB files contributing to the computed angles:
[`final/valid_pdbs.txt`](final/valid_pdbs.txt)

## How to Execute

Ensure all `.pdb.gz` structure files are placed inside a `pdbs/`
directory at the root of the project. Then run:

```bash
snakemake --cores 8
