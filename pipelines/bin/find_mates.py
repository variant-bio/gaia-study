#!/usr/bin/env python3

import pandas as pd
import argparse
from scipy.optimize import linear_sum_assignment


def idealize_csv(input_list, min_obs_threshold=500):
    sample_list = []
    for x in input_list:
        sample_df = pd.read_csv(x)
        sample_df["SAMPLE"] = x.split("-")[0].split("_")[0]
        sample_df["SEQ_COV"] = x.split("-")[0].split("_")[1]
        sample_df["HAP"] = x.split("-")[2].split(".")[0]
        sample_df["REGION"] = x.split("-")[1]
        sample_df["NAME_Q"] = (
            sample_df["SAMPLE"]
            + "_"
            + sample_df["SEQ_COV"].astype(str)
            + "_"
            + sample_df["REGION"]
            + "_"
            + sample_df["Q_COUNT"].astype(str)
        )
        sample_df["NAME_T"] = (
            sample_df["NAME_T"]
            + "_"
            + sample_df["HAP"]
            + "_"
            + sample_df["T_COUNT"].astype(str)
        )
        sample_list.append(sample_df)
    df = pd.concat(sample_list)
    df = df[df["OBS"] >= min_obs_threshold]
    df["ERR"] = df["FP"] + df["FN"]

    check_df = df.set_index(["SAMPLE", "SEQ_COV", "REGION"]).copy()
    df = df.set_index(["NAME_T", "NAME_Q", "ERR"])
    results = []

    for index in check_df.index.unique():
        hung_df = check_df.loc[index].copy()
        cost_matrix = hung_df.pivot_table(
            index="NAME_T", columns="NAME_Q", values="ERR", aggfunc="min"
        )
        # Ensure no NaNs (you could use a very large number instead)
        cost_matrix = cost_matrix.fillna(1e9)
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        for r, c in zip(row_ind, col_ind):
            name_t = cost_matrix.index[r]
            name_q = cost_matrix.columns[c]
            err = cost_matrix.iat[r, c]
            results.append((name_t, name_q, err))

    df = df.loc[results].copy().reset_index()

    return df


# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Process regions and generate summary files."
)
parser.add_argument(
    "--output", "-o", type=str, required=True, help="output directory name"
)
parser.add_argument(
    "--input", "-i", type=str, required=True, nargs="+", help="Input CSV file(s)"
)
parser.add_argument(
    "--min_obs_threshold",
    "-m",
    type=int,
    default=500,
    required=False,
    help="Minimum number of observations to include in the output",
)
args = parser.parse_args()

# Update variables with command line arguments
output = args.output
input = args.input
min_obs_threshold = args.min_obs_threshold

# Call the function with the updated regions list
if __name__ == "__main__":
    combined_df = idealize_csv(input, min_obs_threshold)
    combined_df.to_csv(args.output)
