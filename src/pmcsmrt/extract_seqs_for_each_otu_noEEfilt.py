# -*- coding: utf-8 -*-
"""
Created on Wed May  6 23:23:34 2026

@author: Genglin Guo
@e-mail: gg599@drexel.edu
"""


"""
extract_seqs_for_each_zotu_noEEfilt.py
Python translation of extract_seqs_for_each_otu_noEEfilt.rb

For each zOTU, extracts the read headers from a usearch_global results file
(without EE filtering) and then retrieves those sequences from a FASTQ file.

"""

from pathlib import Path
from Bio import SeqIO

from collections.abc import Iterator

from .log import log

# ─────────────────────────────────────────────
# Minimal FASTA header reader
# ─────────────────────────────────────────────
def read_fasta_definitions(filepath) -> Iterator[str]:
    '''
    Yields the definition line (without '>') for each FASTA record.
    '''
    for record in SeqIO.parse(filepath, 'fasta'):
        yield record.description


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def extract_seqs_noee(post_zOTU_file: Path, post_usearch_glob_results_file: Path, 
                      usearch_global_path: Path, output_dir: Path, verbose: bool):
    '''
    Extracts sequences for each zOTU from the usearch_global file.
    '''
    zotu_list = list(read_fasta_definitions(post_zOTU_file))

    Path(output_dir / 'split_zotus').mkdir(exist_ok = True)

    for description in zotu_list:
        print(description)
        zotu_name = description.split(';')[0]

        headers_file = Path(output_dir) / 'split_zotus' / f'noEE_{zotu_name}_headers.tsv'
        zout_fasta    = Path(output_dir) / 'split_zotus' / f'noEE_{zotu_name}_eefilt_subset.fasta'

        # extract matching headers from uglob_file (replace grep)
        # exact zOTU definition line in the usearch_global file,
        # then keep only the first column (read name)
        matched_ids = []

        with open(post_usearch_glob_results_file) as f:
            for line in f:
                if description in line:
                    matched_ids.append(line.split('\t')[0].rsplit(' ', 1)[0])

        with open(headers_file, 'w') as f:
            f.write('\n'.join(matched_ids) + '\n')

        # retrieve those sequences from the FASTQ file
        with open(zout_fasta, 'w') as f:
            selected_seqs = [record for record in SeqIO.parse(usearch_global_path, 'fastq') if record.id.strip().rsplit(' ', 1)[0] in matched_ids]
            SeqIO.write(selected_seqs, f, 'fasta')
        log(f'Extracted {len(selected_seqs)} sequences for zOTU {zotu_name} (no EE filtering).', verbose = verbose)

