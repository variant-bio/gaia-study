#!/usr/bin/env python3

import pandas as pd
import argparse
import os




def write_csv(input_list, table, output):
    df = pd.concat([pd.read_csv(x) for x in input_list])
    df_dict = {}
    df_dict["sum"] = df[
        ["SAMPLE", "SEQ_COV", "TP", "TN", "FP", "FN", "OBS", "ERR", "REGION", "HAP"]
    ].sort_values(["SAMPLE", "SEQ_COV", "HAP", "REGION"])
    df_dict["per_hap"] = df[["SAMPLE", "SEQ_COV", "TP", "TN", "FP", "FN", "OBS", "ERR", "REGION"]].groupby(
        ["SAMPLE", "SEQ_COV", "REGION"]
    ).mean().reset_index().sort_values(["SAMPLE", "SEQ_COV", "REGION"])
    df_dict["cov_reg"] = df[["SEQ_COV", "TP", "TN", "FP", "FN", "OBS", "ERR", "REGION"]].groupby(
        ["REGION", "SEQ_COV"]
    ).mean().reset_index().sort_values(["SEQ_COV", "REGION"])
    df_dict["sample_cov"] = df[["SEQ_COV", "TP", "TN", "FP", "FN", "OBS", "ERR", "SAMPLE"]].groupby(
        ["SAMPLE", "SEQ_COV"]
    ).mean().reset_index().sort_values(["SAMPLE", "SEQ_COV"])
    df_dict["cov"] = df[["SEQ_COV", "TP", "TN", "FP", "FN", "OBS", "ERR"]].groupby(
        ["SEQ_COV"]
    ).mean().reset_index().sort_values(["SEQ_COV"])
    df_dict["sample"] = df[["TP", "TN", "FP", "FN", "OBS", "ERR", "SAMPLE"]].groupby(
        ["SAMPLE"]
    ).mean().reset_index().sort_values(["SAMPLE"])
    df_dict["reg"] = df[["TP", "TN", "FP", "FN", "OBS", "ERR", "REGION"]].groupby(
        ["REGION"]
    ).mean().reset_index().sort_values(["REGION"])
    # Ensure output directory exists
    if not os.path.exists(output):
        os.makedirs(output)
    if table == "all":
        for tb in df_dict:
            df_dict[tb].to_csv(f"{output}/{tb}.csv", index=False, header=True)
    else:
        df_dict[table].to_csv(f"{output}/{table}.csv", index=False, header=True)


# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Process regions and generate summary files."
)
parser.add_argument(
    "--output", "-o", type=str, required=True, help="output directory name"
)
parser.add_argument(
    "--table", "-t", type=str, required=False, help="Which table to produce", default="all", choices=["sum", "per_hap", "cov_reg", "sample_cov", "cov", "sample", "reg", "all"]
)
parser.add_argument(
    "--input", "-i", type=str, required=True, nargs="+", help="Input CSV file(s)"
)
args = parser.parse_args()

# Update variables with command line arguments
output = args.output
table = args.table
input = args.input

# Call the function with the updated regions list
if __name__ == "__main__":
    write_csv(input, table, output)
