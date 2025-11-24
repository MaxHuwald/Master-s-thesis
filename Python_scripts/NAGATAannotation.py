import pandas as pd
import argparse
import os

def load_and_parse_gff3(filepath, seqid_override=None):
    header_lines = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                line = line.rstrip('\n')
                if seqid_override and line.startswith('##sequence-region'):
                    parts = line.split()
                    if len(parts) >= 2:
                        parts[1] = seqid_override.strip()
                        line = ' '.join(parts)
                header_lines.append(line)
            else:
                break

    df = pd.read_csv(filepath, sep='\t', comment='#', header=None)
    df.columns = [
        "seqid", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes"
    ]

    def parse_attributes(attr_string):
        pairs = attr_string.split(';')
        return {k: v for k, v in (field.split('=') for field in pairs if '=' in field)}
    parsed = df['attributes'].apply(parse_attributes)
    parsed_df = pd.DataFrame(parsed.tolist())
    df = pd.concat([df, parsed_df], axis=1)
    df = df.drop(columns=['attributes'])
    return df, header_lines


def main(nagata_path, reference_path, output_path, seqid_override=None):
    NAGATA, _ = load_and_parse_gff3(nagata_path)  # Input file, no seqid override needed

    reference, reference_header = load_and_parse_gff3(reference_path, seqid_override=seqid_override)

    cds_rows = reference[reference['type'] == 'CDS']
    NAGATA = pd.concat([NAGATA, cds_rows], ignore_index=True)

    mrna_df = NAGATA[NAGATA['type'] == 'mRNA'].copy()
    cds_df = NAGATA[NAGATA['type'] == 'CDS'].copy()

    # Group CDS and mRNA info by ID for start, end, strand
    cds_groups = cds_df.groupby('ID').agg({
        'start': 'min',
        'end': 'max',
        'strand': 'first'
    }).reset_index()

    mrna_groups = mrna_df.groupby('ID').agg({
        'start': 'min',
        'end': 'max',
        'strand': 'first'
    }).reset_index()

    # Map CDS ID to list of parent mRNA IDs
    cds_to_parents = {cds_id: [] for cds_id in cds_groups['ID']}

    ambiguous_mrnas = []

    # For each mRNA, find fully contained CDS groups on same strand, assign parent mRNA IDs
    for _, mrna_row in mrna_groups.iterrows():
        mrna_id = mrna_row['ID']
        mrna_start = mrna_row['start']
        mrna_end = mrna_row['end']
        mrna_strand = mrna_row['strand']

        candidates = cds_groups[
            (cds_groups['strand'] == mrna_strand) &
            (cds_groups['start'] >= mrna_start) &
            (cds_groups['end'] <= mrna_end)
        ]

        if not candidates.empty:
            if mrna_strand == '+':
                min_dist = (candidates['start'] - mrna_start).abs().min()
                best_cds = candidates[(candidates['start'] - mrna_start).abs() == min_dist]
            else:  # reverse strand
                min_dist = (candidates['end'] - mrna_end).abs().min()
                best_cds = candidates[(candidates['end'] - mrna_end).abs() == min_dist]

            # Track ambiguous cases where multiple CDS tie for min distance
            if len(best_cds) > 1:
                ambiguous_mrnas.append((mrna_id, best_cds['ID'].tolist()))

            for cds_id in best_cds['ID']:
                cds_to_parents[cds_id].append(mrna_id)

    # After processing, print ambiguous mRNA ties if any
    if ambiguous_mrnas:
        print("\nWarning: Found mRNAs with multiple equally close CDS candidates:")
        for mrna_id, tied_cds_ids in ambiguous_mrnas:
            print(f"  mRNA ID '{mrna_id}' assigned to multiple CDS IDs: {', '.join(tied_cds_ids)}")

    # Track fallback parents from reference for CDS with no assigned mRNAs
    fallback_parent_ids = set()
    for cds_id, parents in cds_to_parents.items():
        if not parents:
            original_parent = reference.loc[
                (reference['type'] == 'CDS') & (reference['ID'] == cds_id), 'Parent'
            ]
            if not original_parent.empty:
                fallback_parent_ids.add(original_parent.iloc[0])

    # Update 'Parent' attribute for CDS rows in NAGATA
    def parents_for_cds_row(cds_id):
        parents = cds_to_parents.get(cds_id, [])
        if parents:
            return ','.join(parents)
        else:
            # fallback to original Parent in NAGATA (copied from reference)
            original = cds_df.loc[cds_df['ID'] == cds_id, 'Parent']
            if not original.empty:
                return original.iloc[0]
            return ''

    NAGATA.loc[cds_df.index, 'Parent'] = cds_df['ID'].apply(parents_for_cds_row)

    # Append missing Genes from reference
    missing_parents = reference[
        (reference['type'] == 'gene') &
        (reference['ID'].isin(fallback_parent_ids))
    ]
    NAGATA = pd.concat([NAGATA, missing_parents], ignore_index=True)

    # Sort by start, end, type (gene < mRNA < exon < others)
    type_order = {'gene': 0, 'mRNA': 1, 'exon': 2}
    NAGATA['type_order'] = NAGATA['type'].map(type_order).fillna(99)

    NAGATA.sort_values(by=['start', 'end', 'type_order'], ascending=[True, True, True], inplace=True)
    NAGATA.drop(columns=['type_order'], inplace=True)

    cols_after_phase = NAGATA.columns.tolist()
    phase_idx = cols_after_phase.index('phase')
    cols_to_encode = cols_after_phase[phase_idx + 1:]

    # Rebuild attributes column from extra columns
    def rebuild_attributes(row):
        pairs = []
        for col in cols_to_encode:
            val = row[col]
            if pd.notna(val) and val != '':
                pairs.append(f"{col}={val}")
        return ';'.join(pairs)

    NAGATA['attributes'] = NAGATA.apply(rebuild_attributes, axis=1)

    # Keep only original columns
    cols_to_keep = [
        "seqid", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes"
    ]
    NAGATA = NAGATA[cols_to_keep]

    # Override seqid column if specified
    if seqid_override:
        NAGATA['seqid'] = seqid_override

    # Handle output file extension
    def ensure_gff3_extension(path):
        if not path.endswith('.gff3'):
            return path + '.gff3'
        return path

    base_output = output_path
    output_path = ensure_gff3_extension(output_path)

    # Write three output files:
    # 1. full combined output
    # 2. forward strand only (_forward.gff3)
    # 3. reverse strand only (_reverse.gff3)

    # Prepare output filenames
    base_name = os.path.splitext(base_output)[0]
    forward_output = base_name + '_forward.gff3'
    reverse_output = base_name + '_reverse.gff3'

    # Function to write GFF3 file with header and no dataframe header line
    def write_gff3(df, path, header_lines):
        with open(path, 'w') as f:
            for hline in header_lines:
                f.write(hline + '\n')
            df.to_csv(f, sep='\t', header=False, index=False)

    # Write full output
    write_gff3(NAGATA, output_path, reference_header)

    # Filter forward and reverse strands
    forward_df = NAGATA[NAGATA['strand'] == '+']
    reverse_df = NAGATA[NAGATA['strand'] == '-']

    write_gff3(forward_df, forward_output, reference_header)
    write_gff3(reverse_df, reverse_output, reference_header)

    print(f"Formatted files saved to:\n- {output_path}\n- {forward_output}\n- {reverse_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process annotated NAGATA files and merge with reference files.")
    parser.add_argument("-i", "--input", required=True,
                        help="Input NAGATA file path")
    parser.add_argument("-r", "--reference", required=True,
                        help="Reference file path")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file path (without extension)")
    parser.add_argument("--seqid", default=None,
                        help="Optional seqid to override in output (also updates ##sequence-region header)")

    args = parser.parse_args()
    main(args.input, args.reference, args.output, seqid_override=args.seqid)
