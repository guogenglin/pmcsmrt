"""Common filtering, dispatch, SINTAX annotation, reports and read extraction."""
from pathlib import Path
import csv
import subprocess

from Bio import SeqIO

from .dada2_workflow import dada2_clustering
from .usearch_workflow import map_reads, require_sequences, usearch_clustering
from .usearch_get_report import get_report, read_table, write_table
from .usearch_final_parsing import final_parsing
from .merge_all_info_and_sintax import merge_all_info
from .utils import write_to_fastq, PMCSMRTError, read_sequences
from .log import green_bold


def sample_names(args, report_hash: dict) -> list[str]:
    """Include samples whose input FASTQ contains no reads."""
    if args.samplelist:
        with open(args.samplelist) as handle:
            names = [line.rstrip('\n').split('\t')[0] for line in handle
                     if len(line.rstrip('\n').split('\t')) >= 2]
    else:
        names = [path.stem for path in Path(args.foldername).iterdir() if path.is_file()]
    for name in names:
        if not name or any(c in name for c in '\t\r\n'):
            raise PMCSMRTError('Sample IDs must be nonempty and contain no tabs or newlines')
        report_hash[name]  # defaultdict registers empty samples
    return sorted(report_hash)


def prepare_reads(all_reads: dict, trimmed: dict, samples: list[str], args,
                  temp_dir: Path) -> tuple[Path, Path, dict]:
    """Write inference and mapping inputs; DADA2 performs its own EE/length filter."""
    inference = temp_dir / 'pre_sequences_for_clustering.fq'
    mapping = temp_dir / 'pre_sequences_for_usearch_global.fq'
    stats = {s: dict(total=0, ccs_filt=0, primer_filt=0, host_filt=0,
                     size_filt=0, ee_filt=0, unique_filtered=0) for s in samples}
    unique = {s: set() for s in samples}
    for rs in all_reads.values():
        stats[rs.basename]['total'] += 1
        if rs.ccs >= args.ccsvalue:
            stats[rs.basename]['ccs_filt'] += 1
    with inference.open('w') as infer_handle, mapping.open('w') as map_handle:
        for name, tr in trimmed.items():
            rs = all_reads[name]
            st = stats[rs.basename]
            if rs.ccs < args.ccsvalue:
                continue
            st['primer_filt'] += 1
            if rs.host_map == 'True':
                continue
            st['host_filt'] += 1
            size_ok = args.lengthmin <= len(tr.seq) <= args.lengthmax
            if size_ok:
                st['size_filt'] += 1
                write_to_fastq(map_handle, tr.full_name, tr.seq, tr.qual)
            passes = size_ok and tr.ee <= args.eevalue
            if passes:
                st['ee_filt'] += 1
                unique[rs.basename].add(str(tr.seq))
            if args.cluster_method == 'dada2' or passes:
                write_to_fastq(infer_handle, tr.full_name, tr.seq, tr.qual)
    for sample in samples:
        stats[sample]['unique_filtered'] = len(unique[sample])
    return inference, mapping, stats


def sintax(query: Path, database: str, target: Path, threads: int) -> None:
    subprocess.run(['usearch', '-sintax', str(query), '-db', database,
                    '-tabbedout', str(target), '-sintax_cutoff', '0.8',
                    '-strand', 'both', '-threads', str(threads)], check=True)


def split_reads(representatives: Path, mapping: Path, fastq: Path,
                output_dir: Path, feature_type: str, stage: str) -> None:
    """Extract original reads using exact feature IDs, never substring matching."""
    targets = {record.id.split(';')[0] for record in read_sequences(representatives, 'fasta')}
    assignments = {}
    with mapping.open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 2 and fields[1].split(';')[0] in targets:
                assignments[fields[0].split(';')[0]] = fields[1].split(';')[0]
    directory = output_dir / f'split_{feature_type.lower()}s'
    directory.mkdir(exist_ok=True)
    for target in targets:
        (directory / f'{target}_{stage}_subset.fasta').write_text('')
        (directory / f'{target}_{stage}_headers.tsv').write_text('')
    # Open on demand to avoid exhausting file descriptors with many ASVs.
    for record in read_sequences(fastq, 'fastq'):
        target = assignments.get(record.id.split(';')[0])
        if target is None:
            continue
        with (directory / f'{target}_{stage}_subset.fasta').open('a') as fasta:
            SeqIO.write(record, fasta, 'fasta')
        with (directory / f'{target}_{stage}_headers.tsv').open('a') as headers:
            headers.write(record.id + '\n')


def write_step_report(stats: dict, result, samples: list[str], output_dir: Path) -> None:
    """Use read counts for filtering and label unique sequence counts separately."""
    _, counts = read_table(result.table)
    totals = {s: sum(row.get(s, 0) for row in counts.values()) for s in samples}
    if result.feature_type == 'ASV':
        unique = {s: set() for s in samples}
        for record in read_sequences(result.filtered_fastq, 'fastq'):
            attrs = dict(part.split('=', 1) for part in record.id.split(';')[1:] if '=' in part)
            unique[attrs['barcodelabel']].add(str(record.seq))
        for s in samples:
            stats[s]['ee_filt'] = result.stats[s]['filtered']
            stats[s]['unique_filtered'] = len(unique[s])
    keys = ['total', 'ccs_filt', 'primer_filt', 'host_filt', 'size_filt', 'ee_filt', 'unique_filtered']
    with (output_dir / 'post_each_step_report.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample', 'method'] + keys + ['denoised_reads', 'non_chimeric_reads', 'assigned_reads'])
        for s in samples:
            native = result.stats.get(s, {})
            writer.writerow([s, {'ASV': 'dada2', 'OTU': 'uparse', 'zOTU': 'unoise'}[result.feature_type]]
                            + [stats[s][key] for key in keys]
                            + [native.get('denoised', 'NA'), native.get('non_chimeric', 'NA'), totals[s]])


def process_reads(all_fastq: Path, info: Path, all_reads: dict, trimmed: dict,
                  report_hash: dict, args, temp_dir: str, output_dir: Path) -> None:
    """Run one inference method and the shared SINTAX/report workflow."""
    temp = Path(temp_dir)
    samples = sample_names(args, report_hash)
    inference, mapping_reads, stats = prepare_reads(all_reads, trimmed, samples, args, temp)
    print(green_bold(f'Inferring features with {args.cluster_method}...'))
    if args.cluster_method == 'dada2':
        result = dada2_clustering(inference, samples, all_reads, args, temp, output_dir)
        # Native R output is authoritative. Mapping cannot overwrite its table.
        map_table = temp / 'post_ASV_mapping_table.tsv'
    else:
        result = usearch_clustering(inference, args.cluster_method, temp, output_dir, args.threads)
        map_table = result.table
    # DADA2-filtered sequences can differ after quality truncation; use those for
    # the after-filter extraction map, independently of the before-filter map.
    require_sequences(mapping_reads, 'fastq')
    mapping = map_reads(mapping_reads, result.representatives, args.min_identity,
                        temp, output_dir, args.threads, map_table)
    table_samples, counts = read_table(result.table)
    if set(table_samples) - set(samples):
        raise PMCSMRTError('Abundance table contains unexpected sample IDs')
    write_table(result.table, result.feature_type, samples, counts)
    taxa = temp / 'post_reads.sintax'
    sintax(result.representatives, args.taxadbfile, taxa, args.threads)
    annotated = get_report(result, taxa, temp / 'post_unmapped_userach_global.fa', temp, output_dir)
    reference_hits = temp / f'post_{result.feature_type}_blast.txt'
    subprocess.run(['usearch', '-usearch_global', str(result.representatives), '-db', args.taxadbfile,
                    '-top_hit_only', '-id', '0.97', '-blast6out', str(reference_hits),
                    '-strand', 'both', '-evalue', '0.01', '-threads', str(args.threads)], check=True)
    final_parsing(reference_hits, annotated, args.dbclusteredfile, output_dir,
                  args.verbose, result.feature_type)
    sintax(all_fastq, args.taxadbfile, output_dir / 'all_bc_reads.sintax', args.threads)
    merge_all_info(info, output_dir / 'all_bc_reads.sintax', output_dir, args.verbose)
    if args.split_switch == 'yes':
        source = mapping_reads if args.splitmethod == 'before' else result.filtered_fastq
        if args.splitmethod == 'after':
            after_dir = temp / 'after_mapping'
            after_dir.mkdir()
            mapping = map_reads(source, result.representatives, args.min_identity,
                                after_dir, after_dir, args.threads, after_dir / 'table.tsv')
        split_reads(result.representatives, mapping, source, output_dir,
                    result.feature_type, args.splitmethod)
    write_step_report(stats, result, samples, output_dir)
