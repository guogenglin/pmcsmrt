# -*- coding: utf-8 -*-
"""
Created on Sun May 10 15:23:46 2026

@author: Genglin Guo
@e-mail: gg599@drexel.edu
"""


"""
get_fastqs.py — Python translation of get_fastqs.rb

Reads a sample key file, unpacks PacBio barcoded FASTQ archives,
attaches CCS-pass counts and barcode labels to each read header,
and writes one FASTQ per sample into an output folder.

Usage:
    PYTHONPATH=/mnt/f/pmcsmrt/src python -m pmcsmrt.get_fastqs -s sample_key_file.tsv -o output_folder
"""

import os
import argparse
import re
import shutil
import subprocess
import sys
import h5py
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from typing import Generator
from io import TextIOWrapper


# ─────────────────────────────────────────────
# Data class
# ─────────────────────────────────────────────
class Barcode16sRecord:
    """Mirrors Ruby Barcode_16s_record class."""
    def __init__(self, fow_barcode: str | None = None, rev_barcode: str | None = None, sample: str | None = None, 
                 path: str | None = None):
        self.fow_barcode = fow_barcode
        self.rev_barcode = rev_barcode
        self.sample      = sample
        self.path        = path


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────
SPECIAL_CHARS = r"?!<>',?\[\]}{=\-\)\(*&\^%\$#`~{}"
SPECIAL_RE    = re.compile(f"[{SPECIAL_CHARS}]")


def parse_fastq(reads_file: Path | str) -> Generator[tuple[str, Seq, str]]:
    """Generator that yields (header, sequence, quality) tuples from a FASTQ file."""
    for read in SeqIO.parse(reads_file, 'fastq'):
        yield read.description, read.seq.upper(), read.letter_annotations['phred_quality']


def write_to_fastq(fh: TextIOWrapper, header: str, sequence: Seq, quality: str):
    """Write one record in FASTQ format."""
    record = SeqRecord(sequence, id = header, description = '')
    record.letter_annotations['phred_quality'] = quality
    SeqIO.write(record, fh, 'fastq')


# ─────────────────────────────────────────────
# Methods
# ─────────────────────────────────────────────
def prompt_to_change_sample_names(sample_filepath: str) -> dict[str, str]:
    """
    Checks each sample name for leading digits or special characters.
    If any are found, asks the user whether to auto-fix them.
    Returns a dict {sample_name: 'yes'|'no'} indicating which need fixing.
    Mirrors Ruby prompt_to_change_sample_names().
    """
    prompt_hash  = {}
    sample_ind   = None

    with open(sample_filepath) as sf:
        header_split = next(sf).strip('\n').split('\t')
        sample_ind = header_split.index('sample_name')
        for line in sf:
            entry_split = line.strip('\n').split('\t')
            sample = entry_split[sample_ind]
            needs_fix = (bool(re.match(r"^[0-9]", sample)) or bool(SPECIAL_RE.search(sample)))
            prompt_hash[sample] = "yes" if needs_fix else "no"

    if 'yes' in prompt_hash.values():
        print('Sample names need to be changed. '
              'Should I go ahead and change it myself? Answer in yes or no.')
        answer = input().strip()
        if answer == 'yes':
            print('Slightly modifying sample names and moving ahead....')
        elif answer == 'no':
            sys.exit(
                "Correct the sample names manually and re-run this script. "
                "Use sample names without special characters "
                "(underscores are allowed). "
                "Sample names should always begin with an alphabet."
            )
        else:
            sys.exit(
                "You did not answer in yes or no. "
                "Aborting. Run the script again."
            )

    return prompt_hash


def num_of_projects(sample_filepath: str, prompt_hash: dict[str, str]) -> dict[str, list]:
    """
    Parses the sample key file and groups samples by PacBio job ID.
    Returns dict {pb_jobid: [data_path, [Barcode16sRecord, ...]]}.
    Mirrors Ruby num_of_projects().
    """
    pb_projects      = {}
    fow_barcode_ind  = 0
    rev_barcode_ind  = 0
    sample_ind       = 0
    pb_jobid_ind     = 0
    data_path_ind    = 0

    with open(sample_filepath) as sf:
        first_line = next(sf).strip().split('\t')
        fow_barcode_ind = first_line.index('forward_barcode')
        rev_barcode_ind = first_line.index('reverse_barcode')
        sample_ind      = first_line.index('sample_name')
        pb_jobid_ind    = first_line.index('PB_jobid')
        data_path_ind   = first_line.index('data_path')

        for entry in sf:
            cols = entry.rstrip().split('\t')
            job_id = cols[pb_jobid_ind]
            if job_id not in pb_projects:
                pb_projects[job_id] = ['', []]

            rec = Barcode16sRecord()
            rec.fow_barcode = cols[fow_barcode_ind]
            rec.rev_barcode = cols[rev_barcode_ind]

            sample_raw = cols[sample_ind]
            if prompt_hash.get(sample_raw) == "yes":
                # Auto-fix: prepend 'Samp_' if starts with digit
                if re.match(r"^[0-9]", sample_raw):
                    sample_mod = f"Samp_{sample_raw}"
                else:
                    sample_mod = sample_raw
                # Replace any remaining special/non-word chars with '_'
                if SPECIAL_RE.search(sample_mod):
                    sample_mod = re.sub(r"[^\w]", "_", sample_mod)
                rec.sample = sample_mod
            else:
                rec.sample = sample_raw

            pb_projects[job_id][0] = cols[data_path_ind]
            pb_projects[job_id][1].append(rec)

    print(f'Sanity check, number of projects are: {len(pb_projects)}')
    return pb_projects


def get_tarred_folder(samps: list):
    """
    Finds and unpacks the barcoded-fastqs.tgz archive into raw_data/.
    Mirrors Ruby get_tarred_folder().
    """
    data_path = samps[0].rstrip("/")
    curr_file = f"{data_path}/barcoded-fastqs.tgz"
    if not Path(curr_file).exists():
        raise FileNotFoundError(f"The file {curr_file} does not exist!")
    subprocess.run(f"tar xvf {curr_file} -C raw_data/", shell=True, check=True)


def get_ccs_passes(filename):
    """
    Reads a PacBio .ccs.h5 file and yields (read_header, num_passes) pairs.
    Python 3 translation of ccs_passes.py.
    """
    with h5py.File(filename, "r") as hf:
        basename = os.path.basename(filename)[:-9]   # strip '.X.ccs.h5'
 
        suffix = filename[-9:]
        assert suffix in (".1.ccs.h5", ".2.ccs.h5", ".3.ccs.h5"), \
            f"Unexpected suffix: {suffix}"
 
        hole_numbers = hf["/PulseData/ConsensusBaseCalls/ZMW/HoleNumber"][()]
        num_passes_arr = hf["/PulseData/ConsensusBaseCalls/Passes/NumPasses"][()]
 
        for hole_num, num_passes in zip(hole_numbers, num_passes_arr):
            if num_passes > 0:
                yield f"{basename}/{hole_num}/ccs", num_passes


def get_ccs_counts(samps: list):
    """
    Reads CCS pass counts from all .ccs.h5 files under samps[0]
    and returns a dict {read_header: num_passes_str}.
    Mirrors Ruby get_ccs_counts(); calls get_ccs_passes() directly
    instead of shelling out to ccs_passes.py.
    """
    passes_path = Path("passes")
    if passes_path.exists():
        passes_path.unlink()

    data_dir = Path(samps[0].rstrip("/"))
    h5_files = sorted(data_dir.glob("*.ccs.h5"))

    if not h5_files:
        print(f"Warning: no .ccs.h5 files found in {data_dir}", file=sys.stderr)

    with passes_path.open("w") as pf:
        for h5_file in h5_files:
            for read_header, num_passes in get_ccs_passes(str(h5_file)):
                pf.write(f"{read_header}\t{num_passes}\n")

    ccs_hash = {}
    with passes_path.open() as pf:
        for line in pf:
            parts = line.split()
            if len(parts) >= 2:
                ccs_hash[parts[0]] = parts[1]

    return ccs_hash


def process_each_file(samps, output_folder, ccs_hash):
    """
    For each sample record, attaches barcode/CCS tags to the FASTQ headers
    and writes a per-sample FASTQ into output_folder.
    Mirrors Ruby process_each_file().
    """
    with open("log.txt", "w") as log:
        log.write("STYLE: Sample;Barcode\n")

        for rec in samps[1]:
            out_path = Path(output_folder) / f"{rec.sample}.fq"
            if out_path.exists():
                sys.exit(
                    f"Unique sample names must be provided. "
                    f"Sample name {rec.sample} exists multiple times."
                )

            with open(out_path, "w") as all_bc_reads:
                bc  = f"barcodelabel={rec.sample};barcode={rec.fow_barcode}--{rec.rev_barcode};"
                fq  = f"raw_data/{rec.fow_barcode}--{rec.rev_barcode}.fastq"

                if Path(fq).exists():
                    log.write(f"FOUND: {rec.sample};{rec.fow_barcode}--{rec.rev_barcode}\n")
                    print(f"{rec.sample};{rec.fow_barcode}--{rec.rev_barcode}")

                    for header, seq, qual in parse_fastq(fq):
                        raw_name = header.split()[0]
                        ccs_val  = ccs_hash.get(raw_name, "nil")
                        ccs_tag  = f"ccs={ccs_val};"
                        new_header = f"{raw_name};{ccs_tag}{bc}"
                        write_to_fastq(all_bc_reads, new_header, seq.upper(), qual)
                else:
                    log.write(
                        f"NOT FOUND: {rec.sample};"
                        f"{rec.fow_barcode}--{rec.rev_barcode}\n"
                    )
    

# ─────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────
def get_argument() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare per-sample FASTQ files from PacBio barcoded data."
    )
    parser.add_argument("-s", "--samplefile", required = True,
        help="Tab-delimited file with sample information")
    parser.add_argument("-o", "--outfolder", required = True,
        help="Output folder for per-sample FASTQ files")
    return parser


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main(argv: list[str] | None = None):

    argv = sys.argv[1:] if argv is None else argv
    args = get_argument().parse_args(argv)

    if not Path(args.samplefile).exists():
        raise FileNotFoundError('Can\'t open the sample pool file given with \'-s\'!')

    # Create output folder (abort if already exists)
    if Path(args.outfolder).exists():
        raise FileExistsError(
            'Output folder already exists. '
            'Either delete it if not needed, or enter a new folder name.'
        )
    Path(args.outfolder).mkdir(parents = True, exist_ok = False)

    # Check / fix sample names
    prompt_hash = prompt_to_change_sample_names(args.samplefile)

    # Parse sample key file
    pb_projects = num_of_projects(args.samplefile, prompt_hash)
    print(pb_projects)

    # Process each PacBio job
    for samps in pb_projects.values():
        # (Re-)create raw_data directory
        if Path('raw_data').exists():
            shutil.rmtree('raw_data')
        Path('raw_data').mkdir(parents = True, exist_ok = True)

        get_tarred_folder(samps)
        ccs_hash = get_ccs_counts(samps)
        process_each_file(samps, args.outfolder, ccs_hash)

    # Clean up
    shutil.rmtree('raw_data')


if __name__ == "__main__":
    main()