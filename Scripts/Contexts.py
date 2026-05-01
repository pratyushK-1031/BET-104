import sys
import pandas as pd
import re
import os

input_ss = sys.argv[1]
out_path = sys.argv[2]
target_residue = sys.argv[3]

pdb_code = os.path.basename(input_ss).split(".")[0]

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

parsed_data = []

with open(input_ss) as file:
    for line in file:
        if not line.startswith("ASG"):
            continue

        tokens = line.split()
        if len(tokens) < 10:
            continue

        match_num = re.match(r"\d+", tokens[3])
        if not match_num:
            continue

        try:
            parsed_data.append({
                "aa3": tokens[1],
                "chain_id": tokens[2],
                "res_id": int(match_num.group()),
                "index": int(tokens[4]),
                "ss": tokens[5],
                "ss_label": tokens[6],
                "phi_angle": float(tokens[7]),
                "psi_angle": float(tokens[8]),
                "surface_area": float(tokens[9]),
            })
        except Exception:
            continue

df_data = pd.DataFrame(parsed_data)

final_rows = []

if not df_data.empty:
    for j in range(1, len(df_data) - 1):
        mid_res = df_data.iloc[j]

        if mid_res["aa3"] != target_residue:
            continue

        left_res = df_data.iloc[j - 1]
        right_res = df_data.iloc[j + 1]

        aa_triplet = (
            AA3_TO_AA1.get(left_res["aa3"], "X") +
            AA3_TO_AA1.get(mid_res["aa3"], "X") +
            AA3_TO_AA1.get(right_res["aa3"], "X")
        )

        ss_triplet = left_res["ss"] + mid_res["ss"] + right_res["ss"]

        position_info = (
            f"{mid_res['chain_id']}:{left_res['res_id']},"
            f"{mid_res['chain_id']}:{mid_res['res_id']},"
            f"{mid_res['chain_id']}:{right_res['res_id']}"
        )

        for residue in (left_res, mid_res, right_res):
            final_rows.append([
                residue["aa3"],
                residue["chain_id"],
                residue["res_id"],
                residue["index"],
                residue["ss"],
                residue["ss_label"],
                residue["phi_angle"],
                residue["psi_angle"],
                residue["surface_area"],
                AA3_TO_AA1.get(residue["aa3"], "X"),
                aa_triplet,
                ss_triplet,
                pdb_code,
                position_info,
            ])

os.makedirs(os.path.dirname(out_path), exist_ok=True)

output_df = pd.DataFrame(final_rows)
output_df.to_csv(out_path, sep="\t", index=False, header=False)
