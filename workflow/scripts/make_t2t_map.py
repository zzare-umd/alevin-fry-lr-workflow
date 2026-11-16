#!/usr/bin/env python3

import argparse


def generate_t2t_transcriptome(fasta_file, output_file):
    """
    Generate transcript-to-transcript (t2t) map from a transcriptome FASTA.
    Each line: <transcript_id>\t<transcript_id>
    """
    with open(fasta_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                # transcript ID is the first token after '>'
                transcript_id = line[1:].strip().split()[0]
                fout.write(f"{transcript_id}\t{transcript_id}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate transcript-to-transcript (t2t) map from transcriptome FASTA"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Transcriptome FASTA file (e.g. *.fa or *.fasta)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output t2t file"
    )

    args = parser.parse_args()
    generate_t2t_transcriptome(args.input, args.output)


if __name__ == "__main__":
    main()
