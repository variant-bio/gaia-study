#!/usr/bin/env python3


import argparse
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed


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


def cigar_consume(op, length, qorr) -> int:
    """
    Calculate the number of bases consumed by a CIGAR operation.
    Args:
        op (int): CIGAR operation code.
        length (int): Length of the operation.
    Returns:
        int: Number of bases consumed.
    """

    return length * consume_dict[op][qorr]


def trim_cigar(cigar, ref_start, window_start, window_end):
    """
    Adjust the CIGAR string to reflect trimming.
    Args:
        cigar (list): The CIGAR string as a list of tuples (operation, length).
        trim_start (int): Number of bases to trim from the start.
        trim_end (int): Number of bases to trim from the end.
    Returns:
        list: Adjusted CIGAR string.
    """

    # Check forward
    new_cigar = []
    ref_pointer = ref_start
    query_end_index = 0
    for op, length in cigar:
        if op == 5:  # Hard clip
            continue
        ref_pointer_check = ref_pointer + cigar_consume(op, length, "R")
        if ref_pointer_check < window_end:
            new_cigar.append((op, length))
            ref_pointer += cigar_consume(op, length, "R")
            query_end_index += cigar_consume(op, length, "Q")
        elif ref_pointer_check >= window_end:
            new_cg_len = window_end - ref_pointer
            new_cigar.append((op, new_cg_len))
            query_end_index += cigar_consume(op, new_cg_len, "Q")
            break
    # Check reverse
    final_cigar = []
    ref_pointer = window_end
    query_start_index = query_end_index
    new_cigar = reversed(new_cigar)
    for op, length in new_cigar:
        ref_pointer_check = ref_pointer - cigar_consume(op, length, "R")
        if ref_pointer_check > window_start:
            final_cigar.append((op, length))
            ref_pointer -= cigar_consume(op, length, "R")
            query_start_index -= cigar_consume(op, length, "Q")
        elif ref_pointer_check <= window_start:
            new_cg_len = ref_pointer - window_start
            final_cigar.append((op, new_cg_len))
            query_start_index -= cigar_consume(op, new_cg_len, "Q")
            break

    # Remove any operations with zero length
    final_cigar = reversed([(op, length) for op, length in final_cigar if length > 0])

    # need to track query start and end positions
    ## As well as reference start and end positions
    return final_cigar, query_start_index, query_end_index


def trim_bam_by_bed(bam_file, bed_file, output_bam, threads=4):
    # Parse BED file to extract regions
    regions = []
    with open(bed_file, "r", encoding="utf-8") as bed:
        for line in bed:
            if line.strip():
                chrom, start, end = line.strip().split()[:3]
                regions.append((chrom, int(start), int(end)))

    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        header = bam_in.header

    # Prepare region data for parallel processing
    region_data_list = [
        (chrom, start, end, bam_file, consume_dict) for chrom, start, end in regions
    ]

    # Process regions in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(process_region, region_data): region_data
            for region_data in region_data_list
        }

        # Write results to output BAM as they complete
        with pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:
            for future in as_completed(futures):
                processed_reads = future.result()
                for read_dict in processed_reads:
                    read = pysam.AlignedSegment.from_dict(read_dict, header)
                    bam_out.write(read)
            bam_out.close()


def process_region(region_data):
    """Process a single region and return processed reads"""
    chrom, start, end, bam_file, consume_dict = region_data
    processed_reads = []
    q_name_seen = {}

    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        for read in bam_in.fetch(contig=chrom, start=start, stop=end):
            # Calculate trimming positions
            ref_start = read.reference_start
            # Adjust the CIGAR string
            if read.cigartuples is None:
                continue
            new_cigar_tuples, query_start_index, query_end_index = trim_cigar(
                read.cigartuples, ref_start, start, end
            )
            # Convert CIGAR tuples to CIGAR string
            read.cigarstring = "".join(
                f"{length}{consume_dict[op]['CIGAR']}"
                for op, length in new_cigar_tuples
            )
            # Update the read's reference start position
            read.reference_start = max(read.reference_start, start)
            if read.query_name not in q_name_seen:
                q_name_seen[read.query_name] = 1
            else:
                q_name_seen[read.query_name] += 1
                read.query_name = f"{read.query_name}_{q_name_seen[read.query_name]}"
            if read.is_supplementary:
                read.flag = read.flag - 2048
            if read.query_sequence is not None:
                read.query_sequence = read.query_sequence[
                    query_start_index:query_end_index
                ]
            else:
                with pysam.AlignmentFile(bam_file, "rb") as seq_search:
                    for seq_search_read in seq_search.fetch():
                        if (
                            seq_search_read.query_name == read.query_name
                            and seq_search_read.cigarstring is not None
                        ):
                            if (
                                seq_search_read.query_sequence is not None
                                and "H" not in seq_search_read.cigarstring
                            ):
                                read.query_sequence = seq_search_read.query_sequence[
                                    query_start_index:query_end_index
                                ]
                                break
                seq_search.close()
            processed_reads.append(read.to_dict())

    return processed_reads

    # Open the input BAM file and get header for output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Trim BAM file reads based on BED file regions."
    )
    parser.add_argument("-b", "--bam_file", required=False, help="Input BAM file")
    parser.add_argument(
        "-d", "--bed_file", required=False, help="Input BED file with regions"
    )
    parser.add_argument(
        "-o", "--output_bam", required=False, help="Output BAM file with trimmed reads"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of threads for parallel processing",
    )
    args = parser.parse_args()

    trim_bam_by_bed(args.bam_file, args.bed_file, args.output_bam, args.threads)
