# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:10:06 2026

@author: Genglin Guo
@e-mail: gg599@drexel.edu
"""

import os
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from typing import Any
from io import TextIOWrapper

from .log import log

class PMCSMRTError(RuntimeError):
    pass


def read_sequences(path: Path, fmt: str):
    """Yield records while closing the input handle on exhaustion or failure."""
    with Path(path).open() as handle:
        yield from SeqIO.parse(handle, fmt)


def check_cpus(cpus: Any = None, max_cpus: int = 32, verbose: bool = False) -> int:
    '''
    Check the number of CPUs to use for parallel processing. If cpus is None, use all available.
    '''
    avail_cpus = os.cpu_count() or max_cpus
    if isinstance(cpus, str):
        cpus = int(cpus) if cpus.isdigit() else avail_cpus
    elif isinstance(cpus, float):
        cpus = int(cpus)
    else:
        cpus = avail_cpus
    cpus = min(cpus, avail_cpus, max_cpus)
    log(f'Using {cpus=}', verbose)
    return cpus


def check_file(args):
    '''
    Check that all required files exist, and that either a folder of fastqs or a sample list is provided (but not both).
    '''
    if not args.foldername and not args.samplelist:
        raise PMCSMRTError('Must supply a directory of fastqs, or a tab-separated file of sample/filepaths')
    if args.samplelist and args.foldername:
        raise PMCSMRTError('Please only supply either a list of fastqs, or a single directory of all fastqs (not both)')
    
    for path, flag in [
        (args.taxadbfile,       "-d (taxonomy db)"),
        (args.primerfile,       "-p (primer file)"),
        (args.dbclusteredfile,  "-b (databse clustered file)"),
    ]:
        if not path or not Path(path).is_file():
            raise PMCSMRTError(f"Must supply an existing file with {flag}: {path}")

    if args.host_db:
        if not Path(args.host_db).exists():
            raise PMCSMRTError(f"Cannot find host genome file: {args.host_db}")
    log(f'All required parameters were found', verbose = args.verbose)


def check_programs_shutil(progs: list[str], verbose: bool = False):
    '''
    Check if programs are installed and executable using shutil.which.
    '''
    missing = []
    for prog in progs:
        path = shutil.which(prog)
        if path:
            print(f'{prog}: {path}')
        else:
            missing.append(prog)
    if not missing:
        log(f'All required programs found: {", ".join(progs)}', verbose = verbose)
    if missing:
        raise PMCSMRTError(f'Missing programs: {", ".join(missing)}')
    

def dic_check(dictonary: Path) -> list:

    files_in_dic = [f for f in dictonary.iterdir() if f.is_file()]

    if len(files_in_dic) > 1:
        raise PMCSMRTError(f'Expected only one file in {str(dictonary)}, but found {len(files_in_dic)} files.')
    elif len(files_in_dic) == 0:
        raise PMCSMRTError(f'No file found in {str(dictonary)} directory.')
    else:
        return files_in_dic
    
    
def write_to_fastq(fh: TextIOWrapper, header: str, sequence: Seq, quality: str):
    """Write one record in FASTQ format."""
    record = SeqRecord(sequence, id = header, description = '')
    record.letter_annotations['phred_quality'] = [ord(c) - 33 for c in quality]
    SeqIO.write(record, fh, 'fastq')
