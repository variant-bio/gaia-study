#!/usr/bin/env python3


import argparse
import logging
import sys
import pysam
import numpy as np
import pandas as pd


def _init_logger():  # Create a logger named 'app'
    logger = logging.getLogger("app")
    # Set the threshold logging level of the logger to INFO
    logger.setLevel(logging.DEBUG)
    # Create a stream-based handler that writes the log entries    #into the standard output stream
    handler = logging.StreamHandler(sys.stdout)
    # Create a formatter for the logs
    formatter = logging.Formatter("%(levelname)s:%(name)s:%(module)s:%(message)s")
    # Set the created formatter as the formatter of the handler
    handler.setFormatter(formatter)
    # Add the created handler to this logger
    logger.addHandler(handler)


_init_logger()
_logger = logging.getLogger("app")


def cigar_similarity(tread, qread, gff=None, name=None, gff_filter=False, n_window=5):
    """
    Compare two CIGAR strings and their corresponding sequences.
    """
    exon_coords = {}
    _logger.debug("Comparing CIGAR strings and sequences. GFF filter: %s", gff_filter)
    if gff_filter and gff is not None:
        _logger.debug("Loading GFF file: %s", gff)
        gff_df = pd.read_csv(
            gff,
            sep="\t",
            header=None,
            names=[
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
            comment="#",
        )
        _logger.debug("GFF file loaded.")
        gff_df["gene_name"] = gff_df["attributes"].str.extract(r"gene_name=([^;]+)")[0]
        gff_df = gff_df.loc[gff_df["gene_name"] == name].copy()
        exonic_types = ["exon", "CDS"]
        gff_df = gff_df.loc[gff_df["type"].isin(exonic_types)].copy()
        _logger.debug(
            "Filtered GFF file for gene name [%s] with %d entries.",
            name,
            len(gff_df),
        )
        if len(gff_df) == 0:
            _logger.warning(
                "No exonic regions found in GFF file for gene name [%s].",
                name,
            )
            gff_filter = False
        else:
            gff_df["start"] = gff_df["start"].astype(int)
            gff_df["end"] = gff_df["end"].astype(int)
            for _, row in gff_df.iterrows():
                for i in range(row["start"], row["end"] + 1):
                    exon_coords[i] = 1
            _logger.debug(
                "Found %d exonic regions in GFF file for gene name [%s].",
                len(exon_coords),
                name,
            )
    tread_cigar = tread.cigartuples
    tread_cigar_flat = flatten_cigar(tread_cigar)
    tread_seq = tread.query_sequence
    qread_cigar = qread.cigartuples
    qread_cigar_flat = flatten_cigar(qread_cigar)
    qread_seq = qread.query_sequence
    tread_ref_coord = cigar_consume(tread_cigar, tread.reference_start, "R")
    qread_ref_coord = cigar_consume(qread_cigar, qread.reference_start, "R")
    assert len(tread_ref_coord) == len(tread_cigar_flat), (
        "tread_ref_coord and tread_cigar_flat lengths do not match [%s, %s]"
        % (len(tread_ref_coord), len(tread_cigar_flat))
    )
    assert len(qread_ref_coord) == len(qread_cigar_flat), (
        "qread_ref_coord and qread_cigar_flat lengths do not match [%s, %s]"
        % (len(qread_ref_coord), len(qread_cigar_flat))
    )
    tread_q_coord = cigar_consume(tread_cigar, 0, "Q")
    qread_q_coord = cigar_consume(qread_cigar, 0, "Q")
    if tread_seq is None or qread_seq is None:
        return 0, 0, 0, 0, 0
    # Find all indexes where tread_seq has an 'N'
    tread_n_coord = np.array([i for i, base in enumerate(tread_seq) if base == "N"])
    # Add +/-5 integers on either side of each N position
    tread_n_coord = add_n_pad(tread_n_coord, n_window)
    qread_n_coord = np.array([i for i, base in enumerate(qread_seq) if base == "N"])
    qread_n_coord = add_n_pad(qread_n_coord, n_window)
    qread_n_indexes = np.array(
        [i for i, coord in enumerate(qread_q_coord) if coord in qread_n_coord]
    )
    tread_n_indexes = np.array(
        [i for i, coord in enumerate(tread_q_coord) if coord in tread_n_coord]
    )
    tread_skip_indexes = np.array(
        [i for i, base in enumerate(tread_cigar_flat) if base == 3]
    )
    qread_skip_indexes = np.array(
        [i for i, base in enumerate(qread_cigar_flat) if base == 3]
    )
    ref_coord_skip = np.array([])
    for array in [tread_n_indexes, tread_skip_indexes]:
        if len(array) > 0:
            ref_coord_skip = np.append(ref_coord_skip, tread_ref_coord[array])
    for array in [qread_n_indexes, qread_skip_indexes]:
        if len(array) > 0:
            ref_coord_skip = np.append(ref_coord_skip, qread_ref_coord[array])

    if gff_filter:
        # Filter the reference coordinates based on the GFF file
        _logger.debug("Filtering reference coordinates based on GFF file.")
        qread_gff_filt = np.array([coord in exon_coords for coord in qread_ref_coord])
        qread_ref_coord = qread_ref_coord[qread_gff_filt]
        qread_cigar_flat = qread_cigar_flat[qread_gff_filt]
        tread_gff_filt = np.array([coord in exon_coords for coord in tread_ref_coord])
        tread_ref_coord = tread_ref_coord[tread_gff_filt]
        tread_cigar_flat = tread_cigar_flat[tread_gff_filt]
        tread_q_coord = tread_q_coord[tread_gff_filt]
        qread_q_coord = qread_q_coord[qread_gff_filt]
    # Extract the CIGAR strings and sequences
    shared_coords = set(tread_ref_coord) & set(qread_ref_coord)
    shared_coords = shared_coords - set(ref_coord_skip)
    if len(shared_coords) == 0:
        _logger.debug(
            "No matching reference coordinates found between tread and qread."
        )
        return 0, 0, 0, 0, 0
    elif tread_ref_coord[0] == qread_ref_coord[0]:
        tread_ref_pointer, qread_ref_pointer = (
            0,
            0,
        )
    else:
        first_coord = list(shared_coords)[0]
        tread_ref_pointer = np.where(tread_ref_coord == first_coord)[0][0]
        qread_ref_pointer = np.where(qread_ref_coord == first_coord)[0][0]
    # Extract the CIGAR strings and sequences
    observations = 0
    false_neg = 0
    true_neg = 0
    false_pos = 0
    true_pos = 0
    reading = True
    prev_qcigar = None
    prev_tcigar = None
    _logger.debug(
        "tread_ref_coords: %s, qread_ref_coord: %s",
        len(tread_ref_coord),
        len(qread_ref_coord),
    )
    combined_matches = find_all_match_indexes(
        tread_ref_coord, qread_ref_coord, ref_coord_skip
    )
    _logger.debug("Found combined_matches")
    while (
        tread_ref_pointer < len(tread_cigar_flat)
        and qread_ref_pointer < len(qread_cigar_flat)
        and reading
    ):
        _logger.debug(
            "tread_ref_pointer: %s, qread_ref_pointer: %s, tread_ref_coord: %s, qread_ref_coord: %s, reading: %s",
            tread_ref_pointer,
            qread_ref_pointer,
            len(tread_ref_coord),
            len(qread_ref_coord),
            reading,
        )
        # Check if the current bases are 'N' in either sequence
        if tread_cigar_flat[tread_ref_pointer] == qread_cigar_flat[qread_ref_pointer]:
            if tread_cigar_flat[tread_ref_pointer] in [1, 8]:
                if (
                    tread_seq[tread_q_coord[tread_ref_pointer] - 1]
                    != qread_seq[qread_q_coord[qread_ref_pointer] - 1]
                ):
                    false_neg += 1
                    false_pos += 1
                else:
                    true_pos += 1
            elif tread_cigar_flat[tread_ref_pointer] in [0, 7]:
                true_neg += 1
            else:
                true_pos += 1
            observations += 1
        else:
            if (qread_cigar_flat[qread_ref_pointer] == 2 and prev_qcigar == 2) or (
                tread_cigar_flat[tread_ref_pointer] == 2 and prev_tcigar == 2
            ):
                observations += 0
            elif tread_cigar_flat[tread_ref_pointer] in [0, 7]:
                observations += 1
                false_pos += 1
            elif qread_cigar_flat[qread_ref_pointer] in [0, 7]:
                false_neg += 1
                observations += 1
            else:
                false_pos += 1
                false_neg += 1
                observations += 1
        _logger.debug(
            "Began search for next match index: tread_ref_pointer: %s, qread_ref_pointer: %s, tread_ref_coord: %s, qread_ref_coord: %s",
            tread_ref_pointer,
            qread_ref_pointer,
            tread_ref_coord[tread_ref_pointer],
            qread_ref_coord[qread_ref_pointer],
        )
        prev_qcigar = qread_cigar_flat[qread_ref_pointer]
        prev_tcigar = tread_cigar_flat[tread_ref_pointer]
        tread_ref_pointer, qread_ref_pointer, reading = find_next_match_index(
            combined_matches, tread_ref_coord[tread_ref_pointer]
        )
        _logger.debug(
            "Found next match index: tread_ref_pointer: %s, qread_ref_pointer: %s",
            tread_ref_pointer,
            qread_ref_pointer,
        )

    # Return the false positive and false negative rates
    return false_pos, true_pos, false_neg, true_neg, observations


def find_all_match_indexes(list1, list2, ref_coord_skip=None):
    """
    Find all match indexes between two lists.
    """
    shared = set(list1) & set(list2)
    if ref_coord_skip is not None:
        shared = shared - set(ref_coord_skip)
    first1 = {v: np.where(list1 == v)[0][0] for v in shared}
    first2 = {v: np.where(list2 == v)[0][0] for v in shared}
    combined = {v: (first1[v], first2[v]) for v in shared}
    return combined


def find_next_match_index(combined, w):
    """
    Find the next match index in the combined list.
    Args:
        combined (dict): Dictionary of combined indexes.
        w (int): Current index.
    Returns:
        tuple: Next match index and its corresponding values.
    """
    next_values = sorted(v for v in combined if v > w)
    if next_values:
        next_v = next_values[0]
        return combined[next_v][0], combined[next_v][1], True
    else:
        return -1, -1, False


def flatten_cigar(cigar):
    """
    Flatten the CIGAR string into a list of tuples.
    Args:
        cigar (list): CIGAR string as a list of tuples.
    Returns:
        list: Flattened CIGAR string.
    """
    flat_cigar = []
    for op, length in cigar:
        flat_cigar += [op] * length

    return np.array(flat_cigar)


def add_n_pad(narray, n_window=5):
    """
    Add a padding of +/-5 integers to each element in the array.
    Args:
        narray (np.ndarray): Array of integers.
    Returns:
        np.ndarray: Padded array.
    """
    padded_array = np.array([])
    for i in narray:
        for j in range(-n_window, n_window + 1):
            padded_array = np.append(padded_array, i + j)
    return np.unique(padded_array).astype(int)  # Remove duplicates


def cigar_consume(cigar, coord_start, qorr) -> np.ndarray:
    """
    Calculate the number of bases consumed by a CIGAR operation.
    Args:
        op (int): CIGAR operation code.
        length (int): Length of the operation.
    Returns:
        int: Number of bases consumed.
    """
    consume_dict = {
        0: {"CIGAR": "M", "Q": 1, "R": 1},  # Match
        1: {"CIGAR": "I", "Q": 1, "R": 0},  # Insertion
        2: {"CIGAR": "D", "Q": 0, "R": 1},  # Deletion
        3: {"CIGAR": "N", "Q": 0, "R": 1},  # Skipped region
        4: {"CIGAR": "S", "Q": 1, "R": 0},  # Soft clip
        5: {"CIGAR": "H", "Q": 0, "R": 0},  # Hard clip
        6: {"CIGAR": "P", "Q": 0, "R": 0},  # Padding
        7: {"CIGAR": "=", "Q": 1, "R": 1},  # Sequence match
        8: {"CIGAR": "X", "Q": 1, "R": 1},  # Sequence mismatch
    }

    coord_list = [coord_start]
    # Iterate through the CIGAR operations
    for op, length in cigar:
        if op not in consume_dict:
            raise ValueError(f"Invalid CIGAR operation: {op}")
        else:
            coord_list += [consume_dict[op][qorr]] * length

    if coord_list:
        coord_list.pop()

    return np.cumsum(coord_list)


def compare_bam_files(
    truthset_bam,
    comparison_bam,
    output_file,
    gff=None,
    name=None,
    gff_filter=False,
    n_filter=5,
):
    """
    Compare two BAM files: a truthset and a comparison file.
    """
    # Create an output file
    with open(output_file, "w", encoding="utf-8") as out:
        # out.write("NAME_T,NAME_Q,FP,FN,ER,PREC,RECALL,F1,OBS\n")
        out.write("NAME_T,T_COUNT,NAME_Q,Q_COUNT,FP,TP,FN,TN,OBS\n")
        # Open the truthset BAM file
        truthset = pysam.AlignmentFile(truthset_bam, "rb")
        # Iterate through the reads in the truthset BAM file
        for i, tread in enumerate(truthset.fetch()):
            comparison = pysam.AlignmentFile(comparison_bam, "rb")
            for j, qread in enumerate(comparison.fetch()):
                # Check if the read names match
                if tread.query_name != qread.query_name:
                    # Compare the CIGAR strings and sequences
                    fp, tp, fn, tn, obs = cigar_similarity(
                        tread, qread, gff, name, gff_filter, n_filter
                    )
                    out.write(
                        "%s,%s,%s,%s,%s,%s,%s,%s,%s\n"
                        % (
                            tread.query_name,
                            i,
                            qread.query_name,
                            j,
                            fp,
                            tp,
                            fn,
                            tn,
                            obs,
                        )
                    )
            comparison.close()
        truthset.close()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Compare two BAM files: a truthset and a comparison file."
    )
    parser.add_argument(
        "--truthset",
        "-t",
        required=True,
        type=str,
        help="Path to the truthset BAM file.",
    )
    parser.add_argument(
        "--comparison",
        "-c",
        required=True,
        type=str,
        help="Path to the comparison BAM file.",
    )
    parser.add_argument(
        "--gff",
        "-g",
        required=False,
        type=str,
        help="Path to the GFF file [optional].",
    )
    parser.add_argument(
        "--name",
        "-n",
        required=False,
        type=str,
        help="Gene name to filter GFF by",
    )
    parser.add_argument(
        "--n_filter",
        "-f",
        required=False,
        type=int,
        default=5,
        help="Number of bases to skip after an 'N' in the sequence string [default: 5].",
    )
    parser.add_argument(
        "--output", "-o", required=True, type=str, help="Path to the output file."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    GFF_FILTER = False
    # Check if both GFF and gene name are provided
    if args.gff and not args.name:
        _logger.warning(
            "GFF file provided but no gene name specified. Proceeding without filtering."
        )
    if args.name and not args.gff:
        _logger.warning(
            "Gene name provided but no GFF file specified. Proceeding without filtering."
        )
    if args.gff and args.name:
        GFF_FILTER = True
        _logger.info(
            "GFF file and gene name provided. Filtering GFF by gene name [%s].",
            args.name,
        )

    # Compare the BAM files
    compare_bam_files(
        args.truthset,
        args.comparison,
        args.output,
        args.gff,
        args.name,
        GFF_FILTER,
        args.n_filter,
    )
