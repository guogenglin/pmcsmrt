"""Method-neutral feature tables and SINTAX summary reports."""
import csv
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from .utils import PMCSMRTError, read_sequences


def read_table(path: Path) -> tuple[list[str], dict]:
    """Read integer feature counts without interpreting sample names as metadata."""
    with path.open(newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, [])
        if len(header) < 2 or len(set(header[1:])) != len(header[1:]):
            raise PMCSMRTError(f'Invalid sample columns in {path}')
        counts = {}
        for row in reader:
            if not row:
                continue
            if len(row) != len(header):
                raise PMCSMRTError(f'Inconsistent table columns in {path}')
            key = row[0].split(';')[0]
            if key in counts:
                raise PMCSMRTError(f'Duplicate feature {key} in {path}')
            values = [int(value) for value in row[1:]]
            if any(value < 0 for value in values):
                raise PMCSMRTError(f'Negative feature counts in {path}')
            counts[key] = dict(zip(header[1:], values))
    return header[1:], counts


def write_table(path: Path, feature_type: str, samples: list[str], counts: dict) -> None:
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([f'#{feature_type} ID'] + samples)
        for feature, row in counts.items():
            writer.writerow([feature] + [row.get(s, 0) for s in samples])


def get_report(result, sintax_path: Path, unmapped_path: Path,
               temp_dir: Path, output_dir: Path) -> Path:
    samples, counts = read_table(result.table)
    representatives = {r.id.split(';')[0] for r in read_sequences(result.representatives, 'fasta')}
    if not counts or set(counts) - representatives:
        raise PMCSMRTError('Feature table is empty or contains IDs absent from representative FASTA')
    assignments = {}
    with sintax_path.open() as handle:
        for line in handle:
            row = line.rstrip('\n').split('\t')
            if len(row) >= 2:
                assignments[row[0].split(';')[0]] = row[1] or 'NA'
    annotated = temp_dir / f'post_{result.feature_type}_table_sintax_map.tsv'
    with annotated.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([f'#{result.feature_type} ID', 'size'] + samples + ['Assignment'])
        for feature, row in counts.items():
            writer.writerow([feature, sum(row.values())] + [row[s] for s in samples]
                            + [assignments.get(feature, 'NA')])
    unmapped = Counter()
    for record in read_sequences(unmapped_path, 'fasta'):
        attrs = dict(part.split('=', 1) for part in record.id.split(';')[1:] if '=' in part)
        unmapped[attrs.get('barcodelabel', 'Unknown')] += 1
    with (output_dir / 'post_general_report.txt').open('w') as handle:
        handle.write(f'Feature type = {result.feature_type}\n')
        handle.write(f'Number of {result.feature_type}s overall = {len(representatives)}\n')
        handle.write(f'Assigned reads in final abundance table = {sum(sum(r.values()) for r in counts.values())}\n')
        if result.dereplicated:
            unique_count = sum(1 for _ in read_sequences(result.dereplicated, 'fasta'))
            handle.write(f'Unique sequences after dereplication = {unique_count}\n')
        if result.diagnostics and result.diagnostics.exists():
            with result.diagnostics.open() as diagnostics:
                statuses = Counter(line.rstrip('\n').split('\t')[1]
                                   for line in diagnostics if len(line.rstrip('\n').split('\t')) > 1)
            handle.write('Inference diagnostic records by status (unique sequences):\n')
            for status, count in sorted(statuses.items()):
                handle.write(f'\t{status}: {count}\n')
        handle.write('Unmapped eligible reads with USEARCH:\n')
        for sample, count in sorted(unmapped.items()):
            handle.write(f'\t{sample}: {count}\n')
    return annotated
