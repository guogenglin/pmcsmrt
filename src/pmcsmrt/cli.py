# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:00:53 2026

@author: Genglin Guo
@e-mail: gg599@drexel.edu
"""

"""Process PacBio reads with UNOISE, UPARSE or native R DADA2, then SINTAX.

Run from the package parent:
    PYTHONPATH=/mnt/f/pmcsmrt/src python -m pmcsmrt_codex_final.cli \
    -c unoise \
    -f reads \
    -p /mnt/f/pmcsmrt/primers16s.fa \
    -d /mnt/f/pmcsmrt/ncbi_16s_sintax_uniques.fasta \
    -b /mnt/f/pmcsmrt/ncbi_clustered_table.tsv \
    -g /mnt/f/pmcsmrt/all_required_mcsmrt_files/hg38/hg38.fasta \
    -t 24 \
    -V
"""

import argparse
import sys
import shutil
import subprocess
import re
import pysam
import math
from pathlib import Path
import tempfile
from tqdm import tqdm
import pyfastx
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict


from io import TextIOWrapper
from collections.abc import Iterator


from .log import log, bold, formatted_description, green_bold, magenta_bold, yellow
from .workflow import process_reads
from .dada2_workflow import check_dada2
from .utils import check_cpus, check_file, PMCSMRTError, check_programs_shutil


__version__ = '3.0.0'
# ─────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────
def parse_args(a) -> argparse.Namespace:
    """Parse the unified CLI and resolve defaults for the selected method."""
    parser = argparse.ArgumentParser(
        prog='pmcsmrt', description=formatted_description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    io_group = parser.add_argument_group('Input and Output')
    in_out_opts(io_group)
    params = parser.add_argument_group('Parameters')
    params.add_argument('-c', '--cluster-method', choices=('unoise', 'uparse', 'dada2'),
                        default='unoise', help='Feature inference method (default: unoise)')
    params.add_argument('-e', '--eevalue', type=float,
                        help='Maximum expected errors (default: 1.0; dada2: 5.0)')
    params.add_argument('-m', '--min_identity', type=float,
                        help='Read mapping identity (default: 0.99; uparse: 0.97)')
    params.add_argument('-a', '--abundance', type=float, default=2,
                        help='DADA2 chimera minimum parental fold abundance (default: 2)')
    parameter_opts(params)
    params.add_argument('-w', '--split_switch', choices=('yes', 'no'), default='no',
                        help='Split mapped reads into per-feature FASTAs (default: no)')
    params.add_argument('-j', '--splitmethod', choices=('before', 'after'), default='after',
                        help='Split reads before or after EE filtering (default: after)')
    parser.add_argument('-V', '--verbose', action='store_true',
                        help='Keep intermediate files and print debug messages.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    args = parser.parse_args(a)
    if args.eevalue is None:
        args.eevalue = 5.0 if args.cluster_method == 'dada2' else 1.0
    if args.min_identity is None:
        args.min_identity = 0.97 if args.cluster_method == 'uparse' else 0.99
    if not math.isfinite(args.eevalue) or not math.isfinite(args.abundance) or args.eevalue <= 0 or args.abundance < 1:
        parser.error('eevalue must be positive and abundance must be at least 1')
    if not math.isfinite(args.min_identity) or not 0 < args.min_identity <= 1:
        parser.error('min_identity must be in (0, 1]')
    if args.lengthmin < 1 or args.lengthmax < args.lengthmin:
        parser.error('Require 1 <= lengthmin <= lengthmax')
    if args.threads < 1 or args.ccsvalue < 0:
        parser.error('threads must be positive and ccsvalue must be non-negative')
    return args

def in_out_opts(opts: argparse.ArgumentParser):
    '''
    Add input and output options.
    '''
    opts.add_argument('-f', '--foldername', type = str, 
                        help = 'Folder with fastq files, one per sample (-f or -i required)')
    opts.add_argument('-i', '--samplelist', type = str,
                        help = 'List of fastq files (tab-separated sample_name/filepath, -f or -i required)')
    opts.add_argument('-p', '--primerfile', required = True, type = str, 
                    help = 'Primer sequence fasta; headers must be \'forward\' or \'reverse\'')
    opts.add_argument('-d', '--taxadbfile', required = True, type = str,
                        help = 'Fasta version of taxonomy database')
    opts.add_argument('-g', '--host_db', 
                        help = 'Fasta of host genome (can be tar.gz)')
    opts.add_argument('-b', '--dbclusteredfile', required = True, 
                        help = 'Clustered database tsv')
    opts.add_argument('-o', '--output', type = str, default = 'pmcsmrt_output',
                        help = 'The name of the output dict (default: pmcsmrt_output)')


def parameter_opts(opts: argparse.ArgumentParser):
    '''
    Add parameter options.
    '''
    opts.add_argument('-t', '--threads', type = int, default = check_cpus(),
                        help = 'Number of threads (default: all available cpus)')
    opts.add_argument('-s', '--ccsvalue', type = int, default = 5, 
                        help = 'Minimum number of CCS passes (default: 5)')
    opts.add_argument('-x', '--lengthmax', type = int, default = 2000,
                        help = 'Maximum read length (default: 2000)')
    opts.add_argument('-n', '--lengthmin', type = int, default = 500,
                        help = 'Minimum read length (default: 500)')


# ─────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────
class ReadSequence:
    '''
    Holds per-read attributes (mirrors Ruby Read_sequence class).
    '''

    def __init__(self, read_name: str = 'NA', basename: str = 'NA', ccs: int = -1, ee_pretrim: float = -1, 
                 ee_posttrim: float = -1, length_pretrim: int = -1, length_posttrim: int = -1, 
                 host_map: str = 'False', read_orientation: str = 'NA', primer_note: str = 'no_primer_hits'):
        self.read_name          = read_name
        self.basename           = basename
        self.ccs                = ccs
        self.ee_pretrim         = ee_pretrim
        self.ee_posttrim        = ee_posttrim
        self.length_pretrim     = length_pretrim
        self.length_posttrim    = length_posttrim
        self.host_map           = host_map
        self.read_orientation   = read_orientation
        self.primer_note        = primer_note


class TrimmedRead:
    '''
    Holds attributes of a read after trimming.
    '''

    def __init__(self, read_name: str = 'NA', full_name: str = 'NA', seq: Seq = Seq(''), qual: str = '', ee: float = 0.00):
        self.read_name = read_name
        self.full_name = full_name
        self.seq       = seq
        self.qual      = qual
        self.ee        = ee


class Report:
    '''
    Counts reads passing each filter per sample (mirrors Ruby Report class).
    '''

    def __init__(self, total: int = 0, ee_filt: int = 0, size_filt: int = 0, host_filt: int = 0, 
                 primer_filt: int = 0, derep_filt: int = 0, dada2_filt: int = 0, noise_filt: int = 0):
        self.total        = total
        self.ee_filt      = ee_filt
        self.size_filt    = size_filt
        self.host_filt    = host_filt
        self.primer_filt  = primer_filt
        self.derep_filt   = derep_filt
        self.dada2_filt   = dada2_filt
        self.noise_filt   = noise_filt


# ─────────────────────────────────────────────
# FASTQ helpers
# ─────────────────────────────────────────────
def parse_fastq(reads_file: Path) -> Iterator[tuple[str, str]]:
    '''
    Generator that yields (header, sequence) tuples from a FASTQ file.
    '''
    for name, seq, qual, comment in pyfastx.Fastx(str(reads_file), comment = True):
        if comment:
            yield f'{name}_{comment}', seq.upper()
        else:
            yield name, seq.upper()


def process_fastq(fastq_path: str, barcode: str, out_handle: TextIOWrapper):
    '''
    Read a FASTQ file, add/update barcodelabel in the header, and write to out_handle.
    '''
    with open(fastq_path) as handle:
        for record in SeqIO.parse(handle, 'fastq'):

            desc = record.description.replace(' ', '_')  # Replace spaces with underscores in the description

            if 'barcodelabel=' in desc:
                desc = re.sub(
                    r'barcodelabel=([^;]+)',
                    f'barcodelabel={barcode}',
                    desc
                )
            else:
                if desc.endswith(';'):
                    desc += f'barcodelabel={barcode};'
                else:
                    desc += f';barcodelabel={barcode};'

            record.description = desc
            record.id = desc.split()[0]
            record.name = record.id

            SeqIO.write(record, out_handle, 'fastq')
# ─────────────────────────────────────────────
# Core methods
# ─────────────────────────────────────────────
def concat_files(folder_name: str, sample_list: str, temp_dir: str, verbose: bool = False) -> Path:
    '''
    Concatenate all per-sample FASTQ files into one file (pre_demultiplexed_ccsfilt.fq), adding barcodelabel tags.
    '''
    out_path = Path(temp_dir) / 'pre_demultiplexed_ccsfilt.fq'
    with open(out_path, 'w') as out_handle:
        if sample_list is not None:
            with open(sample_list) as f:
                for line in f:
                    arr = line.rstrip('\n').split('\t')
                    if len(arr) < 2:
                        continue
                    barcode = arr[0]
                    fastq_path = arr[1]

                    process_fastq(fastq_path, barcode, out_handle)

        else:
            folder = Path(folder_name.rstrip('/'))
            for fastq_file in folder.iterdir():
                if not fastq_file.is_file():
                    continue
                barcode = fastq_file.stem
                process_fastq(
                    str(fastq_file),
                    barcode,
                    out_handle
                )

    if not Path(out_path).exists() or Path(out_path).stat().st_size == 0:
        raise PMCSMRTError('Error: Concatenating fastq files failed')
    log(f'Concatenated fastq files into {out_path.name}', verbose = verbose)
    return out_path


def get_ee_from_fq_file(all_bc_reads_file: Path, verbose: bool = False) -> dict[str, float]:
    '''
    Extract error rates (EE) from a FASTQ file.
    '''

    ee_hash = {}
    for name, seq, qual, comment in tqdm(pyfastx.Fastx(str(all_bc_reads_file), comment = True), 
                                         desc = 'Calculating EE', 
                                         unit = 'reads'):
        if comment:
            read_name = f'{name}_{comment}'.split(';')[0]
        else:
            read_name = name.split(';')[0]

        # convert Phred quality scores to error probabilities
        ee_hash[read_name] = sum(10 ** (-(ord(c) - 33) / 10.0) for c in qual)
    
    log(f'Extracted EE values for {len(ee_hash)} reads', verbose = verbose)
    return ee_hash


def map_to_human_genome(all_bc_reads_file: Path, human_db: str, temp_dir: str, thread: int, 
                        verbose: bool = False) -> tuple[str, list[str] | None]:
    '''
    Aligns reads to a host genome with bwa mem + sambamba.
    Returns (mapped_count_str, mapped_id_list).
    '''

    # Step 1: bwa mem
    bam_file = Path(temp_dir) / "pre_map_to_host.bam"
    output_fq = Path(temp_dir) / "pre_filt_non_host.fq"

    # Step 1: minimap2 HiFi alignment
    p1 = subprocess.Popen(["minimap2", "-ax", "map-hifi", "-t", str(thread), human_db, 
                           str(all_bc_reads_file)], stdout=subprocess.PIPE)

    try:
        subprocess.run(["samtools", "view", "-b", "-o", str(bam_file), "-"], stdin=p1.stdout, 
                       check=True)
    finally:
        p1.stdout.close()

    if p1.wait() != 0:
        raise RuntimeError("minimap2 alignment failed")

    # Step 2: identify high-confidence human reads
    mapped_set = set()
    host_names = set()

    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary:
                continue

            if read.query_length is None or read.query_length == 0:
                continue

            aligned_length = read.query_alignment_length
            coverage = aligned_length / read.query_length

            if not read.has_tag("NM"):
                continue

            nm = read.get_tag("NM")
            aligned_bases = (read.query_alignment_length + read.get_cigar_stats()[0][2])

            if aligned_bases == 0:
                continue

            identity = 1 - nm / aligned_bases

            if coverage >= 0.90 and identity >= 0.99:
                host_names.add(read.query_name)
                mapped_set.add(read.query_name.split(";")[0])

    # Step 3: recover non-host reads from original FASTQ
    with pysam.FastxFile(str(all_bc_reads_file)) as fq:
        with open(output_fq, "w") as out:
            for record in fq:
                if record.name not in host_names:
                    out.write(str(record) + "\n")

    mapped_count = str(len(mapped_set))
    mapped_id = list(mapped_set)

    log(f"High-confidence human reads: {len(host_names)}", verbose = verbose)

    return mapped_count, mapped_id


def primer_match(all_bc_reads_file: Path, primer_file: str, output_dir: Path, thread: int, 
                 verbose: bool = False) -> tuple[Path, dict[str, TrimmedRead]]:
    '''
    Runs usearch local alignment against primers, then calls primer_parse.rb.
    Returns tuple (pre_trimmed_and_oriented_path, trimmed_length, trimmed_fullhead).
    '''
    forward_primer, reverse_primer = '', ''
    with open(primer_file) as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id.lower() == 'forward':
                forward_primer = str(record.seq)
            elif record.id.lower() == 'reverse':
                reverse_primer = str(record.seq.reverse_complement())
    pre_trimmed_and_oriented_path = Path(output_dir) / 'pre_trimmed_and_oriented.fq'

    cmd = ['cutadapt', '-j', str(thread), '-g', f'{forward_primer}...{reverse_primer}', '--discard-untrimmed', '-o', 
           str(pre_trimmed_and_oriented_path), str(all_bc_reads_file), '--revcomp']

    r = subprocess.run(cmd, check = True)
    
    if r.returncode != 0:
        raise PMCSMRTError("primer cut failed")

    if not pre_trimmed_and_oriented_path.exists() or pre_trimmed_and_oriented_path.stat().st_size == 0:
        raise PMCSMRTError("Error: None of the reads mapped to the primers")
    
    trimmed_hash = {}

    for name, seq, qual, comment in tqdm(pyfastx.Fastx(str(pre_trimmed_and_oriented_path), comment = True),
                    desc = 'Parsing trimmed reads',
                    unit = 'reads'):
        description =  f'{name}_{comment}' if comment else name
        tr = TrimmedRead(
            read_name = description.split(';')[0],
            full_name = description,
            seq = Seq(seq),
            qual = qual,
            ee = sum(10 ** (-(ord(c) - 33) / 10.0) for c in qual)
        )
        trimmed_hash[description.split(';')[0] if comment else name.split(';')[0]] = tr

    log(f'Primer matching completed. Trimmed reads in {pre_trimmed_and_oriented_path.name}', verbose = verbose)
    return pre_trimmed_and_oriented_path, trimmed_hash


def parse_fastq_info(all_bc_reads_file: Path, mapped_id: list[str] | None, trimmed_hash: dict[str, TrimmedRead], 
                     verbose: bool = False) -> tuple[dict[str, ReadSequence], dict[str, ReadSequence], dict[str, Report]]:
    '''
    Parses FASTQ headers to extract read info and populate all_reads_hash, singletons_hash.
    '''
    all_reads_hash, singletons_hash, report_hash = {}, {}, defaultdict(Report)

    count_primer_match = 0
    for header, seq in tqdm(parse_fastq(all_bc_reads_file),
                            desc = 'Parsing FASTQ info',
                            unit = 'reads'):
        def_split = header.split(';')
        read_name = def_split[0]
        
        basename = ''
        ccs = -1

        for item in def_split[1:]:
            k, sep, v = item.partition('=')
            if k == 'barcodelabel':
                basename = v
            elif k == 'ccs':
                ccs = int(v)

            if basename != '' and ccs != -1:
                break
        '''
        if 'barcodelabel' not in attrs or 'ccs' not in attrs:
            raise PMCSMRTError(
                'Header does not have barcodelabel and ccs counts. '
                'Headers should have the format: '
                '@read_name;barcodelabel=sample_name;ccs=ccs_passes'
            )
        '''
        all_reads_hash[read_name] = ReadSequence(
            read_name = read_name, 
            basename = basename, 
            ccs = ccs, 
            length_pretrim = len(seq)
        )

        # Report hash
        report_hash[all_reads_hash[read_name].basename].total += 1

    if mapped_id:
        for map_id in mapped_id:
            all_reads_hash[map_id].host_map = 'True'

    for name in trimmed_hash:
        all_reads_hash[name].primer_note = 'good'
        count_primer_match += 1
        
    log(f'Parsed FASTQ info for {len(all_reads_hash)} reads', verbose = verbose)
    log(f'{count_primer_match} reads had primer matches', verbose = verbose)

    return all_reads_hash, singletons_hash, report_hash


def process_all_bc_reads_file(all_bc_reads_file: Path, output_dir: Path, temp_dir: str, human_db: str, 
                              primer_file: str, thread: int, verbose: bool = False) -> tuple[dict, dict, dict, Path]:
    '''
    Main per-read processing method.
    Returns (all_reads_hash, trimmed_hash, report_hash).
    '''
    all_info_out = open(output_dir / 'pre_all_reads_info.tsv', 'w')
    all_info_out.write(
        'read_name\tbasename\tccs\t'
        'ee_pretrim\tee_posttrim\tlength_pretrim\tlength_posttrim\t'
        'host_map\tread_orientation\tprimer_note\n'
    )

    # Host genome mapping
    print(green_bold('Mapping reads to reference genome...'))
    mapped_count, mapped_id = (map_to_human_genome(all_bc_reads_file, human_db, temp_dir, thread, verbose)
                               if human_db else ('0', []))
    print(green_bold('Done.'))

    # Primer matching
    print(green_bold('Aligning primers...'))
    pre_trimmed_and_oriented_path, trimmed_hash = primer_match(all_bc_reads_file, primer_file, output_dir, thread, verbose)
    print(green_bold('Done.'))

    all_reads_hash = {}

    print(green_bold('Parsing fastq headers and adding info to all_read table...'))
    all_reads_hash, singletons_hash, report_hash = parse_fastq_info(all_bc_reads_file, mapped_id, trimmed_hash, verbose)
    print(green_bold('Done.'))

    # EE pre-trim
    print(green_bold('Calculating EE for all fastq sequences...'))
    ee_pretrim_hash = get_ee_from_fq_file(all_bc_reads_file, verbose)
    print(green_bold('Done.'))

    for k, v in ee_pretrim_hash.items():
        if k in all_reads_hash:
            all_reads_hash[k].ee_pretrim = v

    # Trimming
    print(green_bold('Parsing ee and trimming info for every read...'))
    if not all_reads_hash:
        raise PMCSMRTError('There was a problem parsing, trimming, or orienting the reads')

    for k, rs in all_reads_hash.items():
        if k in trimmed_hash:
            rs.length_posttrim = len(trimmed_hash[k].seq)
            rs.ee_posttrim = trimmed_hash[k].ee
            rs.read_orientation = '-' if trimmed_hash[k].full_name.endswith('rc') else '+'
        else:
            rs.length_posttrim = -1
            rs.ee_posttrim     = -1
            rs.read_orientation = 'NA'

        all_info_out.write('\t'.join(str(x) for x in [
            k, rs.basename, rs.ccs, rs.ee_pretrim, rs.ee_posttrim, rs.length_pretrim, rs.length_posttrim,
            rs.host_map, rs.read_orientation, rs.primer_note
        ]) + '\n')

    print(green_bold('Done.'))
    all_info_out.close()

    if not Path(pre_trimmed_and_oriented_path).exists() or Path(pre_trimmed_and_oriented_path).stat().st_size == 0:
        raise PMCSMRTError('FAIL: No reads written to the pre_trimmed_and_oriented.fq file')
    if not Path(output_dir / 'pre_all_reads_info.tsv').exists() or Path(output_dir / 'pre_all_reads_info.tsv').stat().st_size == 0:
        raise PMCSMRTError('FAIL: No records were written to the \'pre_all_reads_info.tsv\' file.')

    log('All per-read processing completed.', verbose = verbose)
    return all_reads_hash, trimmed_hash, report_hash, Path(output_dir / 'pre_all_reads_info.tsv')


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main(argv: list[str] | None = None):
    '''
    Main method for pmcsmrt.
    '''
    argv = sys.argv[1:] if argv is None else argv
    args = parse_args(argv)

    output_dir = Path.cwd() / args.output
    output_dir.mkdir(exist_ok = True)

    check_file(args)
    programs = ['minimap2', 'samtools'] if args.host_db else []
    check_programs_shutil(programs + ['cutadapt', 'usearch'], verbose=args.verbose)
    if args.cluster_method == 'dada2':
        check_programs_shutil(['Rscript'], verbose=args.verbose)
        check_dada2()

    with tempfile.TemporaryDirectory() as temp_dir:
        # ── Step 1: Concatenate FASTQ files ──────────────────────────────────────
        print(green_bold('Concatenating fastq files..'))
        all_bc_reads_file = concat_files(args.foldername, args.samplelist, temp_dir, args.verbose)
        print(green_bold('Done.'))

        # ── Step 2: Main per-read processing ─────────────────────────────────────
        print(magenta_bold('Starting main method...'))
        all_reads_hash, trimmed_hash, report_hash, pre_all_reads_info_file = process_all_bc_reads_file(
            all_bc_reads_file, output_dir, temp_dir, args.host_db, args.primerfile, args.threads, args.verbose)
        print(magenta_bold('Finished main method.'))

        try:
            process_reads(all_bc_reads_file, pre_all_reads_info_file, all_reads_hash,
                          trimmed_hash, report_hash, args, temp_dir, output_dir)
        finally:
            # Preserve backend diagnostics on failure as well as success.
            if args.verbose:
                shutil.copytree(temp_dir, output_dir / 'mcsmrt_files', dirs_exist_ok=True)

        # ── Step 12: Clean up (if not verbose) ───────────────────────────────────

        print(yellow('mcsmrt successfully completed'))
        log('mcsmrt successfully completed.', verbose = args.verbose)


if __name__ == '__main__':
    raise SystemExit(main())
