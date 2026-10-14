"""USEARCH feature inference and shared read mapping commands."""

from dataclasses import dataclass, field
from pathlib import Path
import subprocess
from Bio import SeqIO
from .utils import PMCSMRTError


@dataclass
class InferenceResult:
    """Backend outputs; counts are feature rows by sample columns."""
    feature_type: str
    representatives: Path
    table: Path
    filtered_fastq: Path
    stats: dict = field(default_factory=dict)
    dereplicated: Path | None = None
    diagnostics: Path | None = None


def require_sequences(path: Path, fmt: str = 'fasta') -> None:
    if not path.exists():
        raise PMCSMRTError(f'No sequences remain in {path.name}; check filtering/inference settings')
    with path.open() as handle:
        if next(SeqIO.parse(handle, fmt), None) is None:
            raise PMCSMRTError(f'No sequences remain in {path.name}; check filtering/inference settings')


def map_reads(reads: Path, representatives: Path, identity: float,
              temp_dir: Path, output_dir: Path, threads: int, table: Path) -> Path:
    """Map eligible reads; the supplied table may be an auxiliary DADA2 table."""
    mapping = output_dir / 'post_usearch_glob_results.tsv'
    subprocess.run([
        'usearch', '-usearch_global', str(reads), '-db', str(representatives),
        '-strand', 'plus', '-id', str(identity), '-maxaccepts', '1', '-top_hit_only',
        '-uc', str(temp_dir / 'post_readmap.uc'), '-otutabout', str(table),
        '-notmatched', str(temp_dir / 'post_unmapped_userach_global.fa'),
        '-userout', str(mapping), '-userfields', 'query+target',
        '-threads', str(threads)], check=True)
    return mapping


def usearch_clustering(filtered: Path, method: str, temp_dir: Path,
                       output_dir: Path, threads: int) -> InferenceResult:
    """Infer UNOISE zOTUs or UPARSE OTUs without altering UNOISE parameters."""
    require_sequences(filtered, 'fastq')
    dereplicated = temp_dir / 'post_dereplicated.fa'
    subprocess.run(['usearch', '-fastx_uniques', str(filtered), '-fastaout',
                    str(dereplicated), '-sizeout', '-threads', str(threads)], check=True)
    feature_type = 'zOTU' if method == 'unoise' else 'OTU'
    representatives = output_dir / f'post_{feature_type}.fa'
    if method == 'unoise':
        diagnostics = temp_dir / 'post_unoise3.txt'
        command = ['usearch', '-unoise3', str(dereplicated), '-zotus',
                   str(representatives), '-tabbedout', str(diagnostics)]
    elif method == 'uparse':
        diagnostics = temp_dir / 'post_uparse.txt'
        command = ['usearch', '-cluster_otus', str(dereplicated), '-otus',
                   str(representatives), '-relabel', 'OTU_', '-uparseout', str(diagnostics), '-minsize', '1']
    else:
        raise PMCSMRTError(f'Unsupported USEARCH method: {method}')
    subprocess.run(command, check=True)
    require_sequences(representatives)
    return InferenceResult(feature_type, representatives,
                           output_dir / f'post_{feature_type}_table.tsv', filtered,
                           dereplicated=dereplicated, diagnostics=diagnostics)
