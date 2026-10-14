"""Join feature counts, SINTAX ranks, reference hits and clustered DB metadata."""
import csv
from pathlib import Path
from .taxonomy import TAX_COLUMNS, taxonomy_fields
from .utils import PMCSMRTError
from .log import log

EXTRA_COLUMNS = TAX_COLUMNS + [
    'blast_query', 'blast_percent_identity', 'blast_alignment_length',
    'ncbi_avglinkage_otu_id', 'num_of_sp_in_cluster',
    'num_of_sp_in_db', 'num_of_otus_with_same_sp']


def final_parsing(reference_hits: Path, annotated_table: Path, db_clust_file: str,
                  output_dir: Path, verbose: bool, feature_type: str = 'zOTU') -> None:
    """Write fixed-width annotations and collapse counts by species."""
    hits = {}
    with reference_hits.open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 4:
                hits.setdefault(fields[0].split(';')[0], fields[1:4])
    clustered = {}
    with open(db_clust_file) as handle:
        next(handle, None)
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 15:
                raise PMCSMRTError('Clustered database rows must have at least 15 columns')
            species = fields[12]
            if species in clustered:
                clustered[species][2] += 1
            else:
                clustered[species] = [fields[1], fields[-3], 1, fields[-1]]
    collapsed, members = {}, {}
    with annotated_table.open(newline='') as source, (output_dir / 'post_final_results.tsv').open('w', newline='') as dest:
        reader = csv.reader(source, delimiter='\t')
        header = next(reader)
        count_columns = header[1:-1]  # size and sample columns, positional to allow arbitrary sample IDs
        writer = csv.writer(dest, delimiter='\t')
        writer.writerow(header[:-1] + EXTRA_COLUMNS)
        for row in reader:
            if len(row) != len(header):
                raise PMCSMRTError('Annotated feature table has inconsistent columns')
            tax = taxonomy_fields(row[-1])
            species = tax[12]
            writer.writerow(row[:-1] + tax + hits.get(row[0], ['NA'] * 3)
                            + clustered.get(species, ['NA'] * 4))
            values = [int(value) for value in row[1:-1]]
            if species not in collapsed:
                collapsed[species] = [0] * len(values)
                members[species] = []
            collapsed[species] = [a + b for a, b in zip(collapsed[species], values)]
            members[species].append(row[0])
    totals = [sum(values[i] for values in collapsed.values()) for i in range(len(count_columns))]
    with (output_dir / 'post_taxa_collapse.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['Species'] + count_columns + [c + '_rel' for c in count_columns]
                        + [feature_type + '_member'])
        for species, values in collapsed.items():
            relative = [value / total * 100 if total else 0.0 for value, total in zip(values, totals)]
            writer.writerow([species] + values + relative + [';'.join(members[species])])
    log('Final annotations and species collapse written.', verbose=verbose)
