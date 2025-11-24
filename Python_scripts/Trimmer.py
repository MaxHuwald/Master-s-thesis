import pysam
import argparse

def trim_polyA(seq, qual, max_wiggle=2):
    wiggle_count = 0
    trim_pos = len(seq)
    for i in range(len(seq) - 1, -1, -1):
        if seq[i] == 'A':
            wiggle_count = 0
            trim_pos = i
        else:
            wiggle_count += 1
            if wiggle_count > max_wiggle:
                break
    return seq[:trim_pos], qual[:trim_pos]

def main(input_bam, output_bam, max_wiggle):
    in_bam = pysam.AlignmentFile(input_bam, "rb", check_sq=False)
    out_bam = pysam.AlignmentFile(output_bam, "wb", template=in_bam)

    for read in in_bam.fetch(until_eof=True):
        seq = read.query_sequence
        qual = read.query_qualities

        if seq:
            trimmed_seq, trimmed_qual = trim_polyA(seq, qual, max_wiggle=max_wiggle)
            trimmed_bases = len(seq) - len(trimmed_seq)

            read.query_sequence = trimmed_seq
            read.query_qualities = trimmed_qual

            # Add tag for trimmed bases (XB)
            read.set_tag("XB", trimmed_bases, value_type='i')  # integer tag

        out_bam.write(read)

    in_bam.close()
    out_bam.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim polyA tails from BAM")
    parser.add_argument("--input", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output trimmed BAM file")
    parser.add_argument("--max_wiggle", type=int, default=2, help="Max number of non-A bases allowed in polyA tail trimming (default: 2)")
    args = parser.parse_args()

    main(args.input, args.output, args.max_wiggle)
