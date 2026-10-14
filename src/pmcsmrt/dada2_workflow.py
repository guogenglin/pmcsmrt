"""Invoke native R DADA2 on cutadapt-trimmed, oriented PacBio FASTQs."""
import csv
from pathlib import Path
import subprocess
from Bio import SeqIO
from .usearch_workflow import InferenceResult, require_sequences
from .utils import PMCSMRTError, read_sequences

R_SCRIPT = Path(__file__).with_name('run_dada2.R')


def check_dada2() -> None:
    """Check R, dada2 and its PacBio estimator, and print installed versions."""
    subprocess.run(['Rscript', str(R_SCRIPT), '--check'], check=True)


def dada2_clustering(reads: Path, samples: list[str], all_reads: dict,
                     args, temp_dir: Path, output_dir: Path) -> InferenceResult:
    """Write safe sample filenames and a manifest; preserve original sample IDs."""
    require_sequences(reads, 'fastq')
    demux = temp_dir / 'dada2_input'
    demux.mkdir()
    sample_paths = {sample: demux / f'sample_{i:06d}.fastq'
                    for i, sample in enumerate(samples)}
    handles = {}
    try:
        for record in read_sequences(reads, 'fastq'):
            sample = all_reads[record.id.split(';')[0]].basename
            if sample not in handles:
                handles[sample] = sample_paths[sample].open('w')
            SeqIO.write(record, handles[sample], 'fastq')
    finally:
        for handle in handles.values():
            handle.close()
    manifest = temp_dir / 'dada2_manifest.tsv'
    with manifest.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample', 'fastq'])
        for sample, path in sample_paths.items():
            writer.writerow([sample, str(path.resolve()) if path.exists() else ''])
    representatives = output_dir / 'post_ASV.fa'
    table = output_dir / 'post_ASV_table.tsv'
    filtered = temp_dir / 'dada2_filtered.fq'
    stats_path = output_dir / 'post_dada2_stats.tsv'
    command = ['Rscript', str(R_SCRIPT), str(manifest), str(temp_dir.resolve()),
               str(representatives.resolve()), str(table.resolve()),
               str(stats_path.resolve()), str(filtered.resolve()),
               str((output_dir / 'post_base_transition_stats.tsv').resolve()),
               str(args.lengthmin), str(args.lengthmax), str(args.eevalue),
               str(args.abundance), str(args.threads)]
    log_path = temp_dir / 'dada2.log'
    with log_path.open('w') as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        tail = log_path.read_text(errors='replace')[-6000:]
        raise PMCSMRTError(f'Native R DADA2 failed (exit {result.returncode}):\n{tail}')
    require_sequences(representatives)
    with stats_path.open(newline='') as handle:
        stats = {row['sample']: {key: int(value) for key, value in row.items()
                                if key != 'sample'}
                 for row in csv.DictReader(handle, delimiter='\t')}
    if set(stats) != set(samples):
        raise PMCSMRTError('DADA2 statistics did not retain all input samples')
    return InferenceResult('ASV', representatives, table, filtered, stats=stats)
