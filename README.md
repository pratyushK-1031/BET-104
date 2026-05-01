# BET-104 Assignment

This workflow analyzes every Arginine residue located within an
alpha-helix triplet (HHH). For each case, it calculates the signed 3D
angle formed between neighboring C-α → centroid vectors, and then
categorizes the angle distribution based on the size group of the
preceding amino acid.

## Visualization

<img width="3000" height="1800" alt="plot_for_arg_with_valid_runs" src="https://github.com/user-attachments/assets/2963e40c-c67a-4868-a110-0a15fefaefaf" />

## How to Execute

Ensure all `.pdb.gz` structure files are placed inside a `pdbs/`
directory at the root of the project. Then run:

```bash
snakemake --cores 8
