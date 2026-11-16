#!/usr/bin/env python3

import pysam
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Add CR, CB, UB and UR tags to a BAM file by parsing CB and UMI from read names"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input BAM file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output BAM file with CR/CB/UR/UB tags added"
    )

    args = parser.parse_args()

    bam_in = pysam.AlignmentFile(args.input, "rb")
    bam_out = pysam.AlignmentFile(args.output, "wb", template=bam_in)

    for read in bam_in:
        qname = read.query_name

        # CB_UMI#original_read_id_strand
        try:
            cb_umi, _ = qname.split("#", 1)
            cb, umi = cb_umi.split("_", 1)
        except ValueError:
            # If read name doesn't match expected pattern, keep read unchanged
            bam_out.write(read)
            continue

        # Add raw tags
        read.set_tag("CR", cb, value_type="Z")
        read.set_tag("CB", cb, value_type="Z")
        read.set_tag("UR", umi, value_type="Z")
        read.set_tag("UB", umi, value_type="Z")

        bam_out.write(read)

    bam_in.close()
    bam_out.close()

    print(f"[DONE] CR/CB/UR/UB tagging finished successfully.")

if __name__ == "__main__":
    main()
