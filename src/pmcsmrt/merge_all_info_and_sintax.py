# -*- coding: utf-8 -*-
"""
Created on Wed May  6 23:21:38 2026

@author: Genglin Guo
@e-mail: gg599@drexel.edu
"""


"""
merge_all_info_and_sintax.py — Python translation of merge_all_info_and_sintax.rb

Merges per-read info (pre_all_reads_info.tsv) with SINTAX taxonomic
classifications (all_bc_reads.sintax) into a single output file
(pre_all_bc_reads_info_with_taxa.tsv).

"""

from pathlib import Path

from .log import log
from .taxonomy import taxonomy_fields


# ─────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────
SINTAX_HEADER = (
    'Domain\tDomain_conf\tPhylum\tPhylum_conf\tClass\tClass_conf\t'
    'Order\tOrder_conf\tFamily\tFamily_conf\tGenus\tGenus_conf\t'
    'Species\tSpecies_conf'
)

# 14 NA columns — one per SINTAX_HEADER field
NA_TAXA = '\t'.join(['NA'] * 14)

# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def merge_all_info(pre_all_reads_info_file: Path, all_bc_reads_sintax_file: Path, output_dir: Path, verbose: bool = False) -> None:
    '''
    Merges per-read info with SINTAX taxonomic classifications.
    '''
    # ── Build SINTAX hash ───────────────────────────────────────────────────────
    # key  : read name (before first ';')
    # value: tab-joined taxonomy capture groups
    sintax_hash = {}
    with open(all_bc_reads_sintax_file) as sintax_fh:
        for line in sintax_fh:
            line_split = line.strip('\n').split('\t')
            if len(line_split) < 2:
                continue
            key = line_split[0].split(';')[0]
            sintax_hash[key] = '\t'.join(taxonomy_fields(line_split[1]))

    # ── Build all-reads hash ──────────────────────────────────────────────────
    # key  : read name (column 0)
    # value: remaining columns joined by tab
    all_reads_hash = {}
    info_header = ""
    with open(pre_all_reads_info_file) as info_fh:
        info_header = next(info_fh).strip('\n')
        for line in info_fh:
            line_split = line.split('\t')
            all_reads_hash[line_split[0]] = '\t'.join(line_split[1:]).strip('\n')

    # ── Write merged output ───────────────────────────────────────────────────
    with open(Path(output_dir) / 'pre_all_bc_reads_info_with_taxa.tsv', 'w') as out:
        out.write(f'{info_header}\t{SINTAX_HEADER}\n')
        for key, value in all_reads_hash.items():
            if key in sintax_hash:
                out.write(f'{key}\t{value}\t{sintax_hash[key]}\n')
            else:
                out.write(f'{key}\t{value}\t{NA_TAXA}\n')

    log('Merged pre_all_reads_info.tsv with all_bc_reads.sintax into pre_all_bc_reads_info_with_taxa.tsv', verbose = verbose)
