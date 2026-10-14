# -*- coding: utf-8 -*-
"""
Created on Wed May  6 23:25:37 2026

@author: Genglin Guo
@e-mail: gg599@drexel.edu
"""


"""
extract_seqs_for_each_zotu.py
Python translation of extract_seqs_for_each_otu.rb
 
For each zOTU, extracts matching read headers from a uparse .up file,
then retrieves those sequences from a dereplicated FASTA file.

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
    Yield the definition line (without '>') for each FASTA record.
    '''
    for record in SeqIO.parse(filepath, 'fasta'):
        yield record.description
 
 
# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def extract_seqs(post_zOTU: Path, post_unoise_file: Path, post_dereplicated_file: Path, 
                 post_usearch_glob_results_file: Path, output_dir: Path, verbose: bool):
    '''
    Extracts sequences for each zOTU from the usearch_global file.
    '''
    # Ruby: entry.definition.split(";")[0]  — strip size/other tags from header
    zotu_list = [d.split(';')[0] for d in read_fasta_definitions(post_zOTU)]

    Path(output_dir / 'split_zotus').mkdir(exist_ok = True)
 
    for description in zotu_list:
        print(description)

        headers_file = Path(output_dir) / 'split_zotus' / f'{description}_headers.tsv'
        out_fasta    = Path(output_dir) / 'split_zotus' / f'{description}_eefilt_subset.fasta'
 
        # extract matching headers from uglob_file (replace grep)
        # exact zOTU definition line in the usearch_global file,
        # then keep only the first column (read name)

        matched_ids = []

        with open(post_usearch_glob_results_file) as f:
            for line in f:
                if description in line:
                    matched_ids.append(line.split('\t')[0].rsplit(' ', 1)[0])

        filtered_out_ids = []
        with open(post_unoise_file) as f:
            for line in f:
                filtered_out_ids.append(line.strip().split('\t')[0].rsplit(' ', 1)[0])
 
        # retrieve those sequences from the dereplicated FASTA
        count = 0
        with open(out_fasta, 'w') as handle1, open(headers_file, 'w') as handle2:
            headers_file, selected_seqs = [], []
            for record in SeqIO.parse(post_dereplicated_file, 'fasta'):

                if record.id.strip().rsplit(' ', 1)[0] not in filtered_out_ids and record.id.strip().rsplit(' ', 1)[0] in matched_ids:
                    SeqIO.write(record, handle1, 'fasta')
                    handle2.write(f'{record.id}\n')
                    count += 1
        log(f'Extracted {count} sequences for zOTU {description}.', verbose = verbose)