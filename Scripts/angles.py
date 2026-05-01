import gzip
import multiprocessing as mp
import os
import sys
from functools import partial

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from tqdm import tqdm


AMINO_SIZE = {
    "K": "Large", "R": "Bulky", "D": "Intermediate", "E": "Large",
    "G": "Tiny", "A": "Tiny", "V": "Small", "L": "Intermediate", "I": "Intermediate",
    "M": "Large", "P": "Small", "S": "Small", "T": "Small", "N": "Intermediate",
    "Q": "Large", "C": "Small", "F": "Bulky", "Y": "Bulky", "W": "Bulky", "H": "Large",
}


def calculate_signed_angle(vec_a, vec_b, ref_axis, deg=True):
    vec_a = np.array(vec_a, dtype=float)
    vec_b = np.array(vec_b, dtype=float)
    ref_axis = np.array(ref_axis, dtype=float)

    vec_a /= np.linalg.norm(vec_a)
    vec_b /= np.linalg.norm(vec_b)
    ref_axis /= np.linalg.norm(ref_axis)

    cross_prod = np.cross(vec_a, vec_b)
    dot_prod = np.dot(vec_a, vec_b)

    theta = np.arctan2(np.dot(cross_prod, ref_axis), dot_prod)

    return np.degrees(theta) if deg else theta


def read_structure(pdb_code, pdb_folder):
    file_loc = os.path.join(pdb_folder, f"{pdb_code}.pdb.gz")
    if not os.path.exists(file_loc):
        return None

    parser = PDBParser(QUIET=True)
    try:
        with gzip.open(file_loc, "rt") as f:
            return parser.get_structure(pdb_code, f)
    except Exception:
        return None


def fetch_residue(structure_obj, chain_label, residue_num):
    try:
        for mdl in structure_obj:
            ch = mdl[chain_label]
            for residue in ch:
                if residue.id[1] == residue_num:
                    return residue
    except Exception:
        return None
    return None


def extract_ca(residue):
    if residue and "CA" in residue:
        return residue["CA"].get_coord()
    return None


def compute_centroid(residue):
    atom_coords = []
    for atom in residue:
        if atom.get_name() not in ["N", "CA", "C", "O"]:
            atom_coords.append(atom.get_coord())
    if not atom_coords:
        return None
    return np.mean(atom_coords, axis=0)


def handle_single_file(file_name, context_path, pdb_path, target_residue):
    file_full_path = os.path.join(context_path, file_name)

    if os.path.getsize(file_full_path) == 0:
        return ([], None)

    try:
        data = pd.read_csv(file_full_path, sep="\t", header=None)
    except Exception:
        return ([], None)

    if data.empty:
        return ([], None)

    pdb_code = data.iloc[0, 12]
    structure_obj = read_structure(pdb_code, pdb_path)
    if structure_obj is None:
        return ([], None)

    result_rows = []

    for idx in range(0, len(data), 3):
        try:
            res_prev = data.iloc[idx]
            res_mid = data.iloc[idx + 1]
            res_next = data.iloc[idx + 2]
        except Exception:
            continue

        if res_mid[0] != target_residue:
            continue
        if res_mid[11] != "HHH":
            continue

        chain_id = res_mid[1]
        prev_num = int(res_prev[2])
        mid_num = int(res_mid[2])
        next_num = int(res_next[2])

        r1 = fetch_residue(structure_obj, chain_id, prev_num)
        r2 = fetch_residue(structure_obj, chain_id, mid_num)
        r3 = fetch_residue(structure_obj, chain_id, next_num)

        if not r1 or not r2 or not r3:
            continue

        ca1 = extract_ca(r1)
        ca2 = extract_ca(r2)
        ca3 = extract_ca(r3)

        cen1 = compute_centroid(r1)
        cen2 = compute_centroid(r2)

        if any(x is None for x in [ca1, ca2, ca3, cen1, cen2]):
            continue

        vec1 = cen1 - ca1
        vec2 = cen2 - ca2
        axis_vec = ca2 - ca1

        if np.linalg.norm(axis_vec) == 0:
            continue

        try:
            ang = calculate_signed_angle(vec1, vec2, axis_vec)
        except Exception:
            continue

        prev_aa = res_prev[9]
        size_type = AMINO_SIZE.get(prev_aa, "Unknown")

        if size_type == "Unknown":
            continue

        result_rows.append([pdb_code, prev_aa, size_type, ang])

    return (result_rows, pdb_code if result_rows else None)


def run_pipeline():
    context_folder = sys.argv[1]
    output_path = sys.argv[2]
    target_res = sys.argv[3]
    pdb_folder = sys.argv[4] if len(sys.argv) > 4 else "pdbs"
    num_workers = int(sys.argv[5]) if len(sys.argv) > 5 else 12

    file_list = sorted(f for f in os.listdir(context_folder) if f.endswith(".tsv"))

    worker_func = partial(
        handle_single_file,
        context_path=context_folder,
        pdb_path=pdb_folder,
        target_residue=target_res,
    )

    collected_data = []
    valid_ids = set()

    with mp.Pool(num_workers) as pool:
        for rows, pdb_code in tqdm(
            pool.imap_unordered(worker_func, file_list, chunksize=16),
            total=len(file_list),
            desc="Processing files",
        ):
            if rows:
                collected_data.extend(rows)
            if pdb_code is not None:
                valid_ids.add(str(pdb_code))

    df_out = pd.DataFrame(
        collected_data,
        columns=["pdb_id", "previous_aa", "size_group", "angle_value"]
    )

    save_dir = os.path.dirname(output_path) or "."
    os.makedirs(save_dir, exist_ok=True)

    df_out.to_csv(output_path, sep="\t", index=False)

    valid_file = os.path.join(save_dir, "valid_pdb_ids.txt")
    with open(valid_file, "w") as f:
        for pid in sorted(valid_ids):
            f.write(f"{pid}\n")

    print(f"Output saved: {output_path}")
    print(f"Valid PDB list saved: {valid_file}")
    print(f"Total computed angles: {len(df_out)}")
    print(f"Total valid PDBs: {len(valid_ids)}")


if __name__ == "__main__":
    run_pipeline()
