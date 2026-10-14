# PMCSMRT

PMCSMRT is a Python pipeline for processing PacBio full-length 16S rRNA amplicon reads, inferring OTUs, zOTUs or ASVs, and assigning taxonomy with USEARCH SINTAX. It produces representative sequences, per-sample abundance tables, species summaries, and detailed read-level and filtering reports.

# Preface

This project is an upgraded Python version of [MCSMRT](https://github.com/jpearl01/mcsmrt), originally developed by Joshua P. Earl and colleagues. The original tool established a workflow for clustering PacBio full-length 16S reads and combining taxonomic assignments with read-specific information.

In this version, the Ruby workflow has been translated into Python, USEARCH has been updated from version 8 to version 12, and additional feature inference methods have been integrated. The aim is to retain the useful reporting features of MCSMRT while making the workflow easier to run and extend.

For questions about the original methodology, please refer to the original repository and publication. For questions or feedback about this modified version, you are welcome to contact me.

# Major Improvements

## Python implementation

The workflow and supporting scripts now use Python 3, replacing the original Ruby scripts. Shared preprocessing, taxonomy parsing, reporting and sequence extraction are used across all three inference methods.

## USEARCH 12 and SINTAX

The supplied Conda environment includes USEARCH `12.0_beta`. Taxonomic classification uses a SINTAX-formatted reference FASTA rather than the original USEARCH 8 UTAX UDB workflow. Both representative sequences and input reads receive taxonomic annotations.

## Three feature inference methods

| Method | Implementation | Feature type | Counting method |
|---|---|---|---|
| `unoise` (default) | USEARCH UNOISE3 | zOTU | Eligible reads mapped to representatives |
| `uparse` | USEARCH UPARSE | OTU | Eligible reads mapped to representatives |
| `dada2` | Native R DADA2 with `PacBioErrfun` | ASV | Native counts after denoising and chimera removal |

DADA2 runs directly through `Rscript`; QIIME2 and QZA/QZV artifacts are not required. Its native abundance table is preserved independently of the additional USEARCH mapping results.

## Shared preprocessing and reports

Cutadapt matches linked primers, trims reads and handles reverse-complement orientation. Optional host filtering uses minimap2 and samtools. Reports distinguish cumulative read filtering counts, unique sequence counts and final assigned reads. Reads can also be exported into separate FASTA files for each inferred feature.

# Dependencies and Installation

## Conda environment:

The supplied [pmcsmrt.yml](pmcsmrt.yml) contains Linux-specific packages and build identifiers. Use it in Linux or inside WSL on Windows.

Clone the repository into a directory named `pmcsmrt`:

```bash
git clone https://github.com/guogenglin/pmcsmrt.git pmcsmrt
cd pmcsmrt
conda env create -f pmcsmrt.yml
conda activate usearch
cd ..
PYTHONPATH=/path/to/package_parent python -m pmcsmrt.cli --help
```

The environment is named **`usearch`**, as specified in the YAML file. Its main dependencies include:

| Dependency | Version in `pmcsmrt.yml` | Purpose |
|---|---|---|
| Python | 3.11.15 | Pipeline and supporting scripts |
| USEARCH | 12.0_beta | UNOISE, UPARSE, mapping and SINTAX |
| cutadapt | 5.2 | Primer matching, trimming and orientation |
| minimap2 | 2.31 | Optional host alignment |
| samtools | 1.23.1 | Optional host alignment processing |
| R | 4.5.3 | DADA2 backend |
| DADA2 | 1.38.0 | PacBio error learning and ASV inference |

The YAML already includes DADA2 and the main Python dependencies. `h5py`, needed by the optional `get_fastqs` utility, is installed separately as described below. These versions describe the supplied environment specification; successful environment creation has not been independently verified in this checkout.

## Manual setup

Alternatively, use Python 3.10+ and install the Python libraries:

```bash
python -m pip install biopython pandas tqdm pyfastx pysam h5py
```

Install `usearch` and `cutadapt` and make them available on `PATH`. The USEARCH executable must be named `usearch`. Host filtering additionally requires `minimap2` and `samtools`.

For DADA2, install R with `Rscript` on `PATH`, then install DADA2 in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("dada2")
```

The pipeline checks that DADA2 provides `PacBioErrfun`. USEARCH is required even for DADA2 runs because mapping and taxonomy annotation use USEARCH. Neither `blastn` nor `biom` is required.

# Input Requirements

## Per-sample FASTQ files

Provide demultiplexed CCS reads, one plain-text FASTQ file per sample. Choose either a directory (`-f`) or a sample list (`-i`), but not both.

FASTQ headers must include a numeric CCS pass count:

```text
@movie_name/123/ccs;ccs=17;
```

The pipeline adds or updates `barcodelabel` using the filename stem or the sample name in the input list. Read identifiers before the first semicolon must be unique across all samples. Missing CCS annotations fail the default CCS filter; nonnumeric values cause a parsing error. Pass counts are not inferred from quality scores.

For directory input, include only sample FASTQs. For example, `sample_A.fq` becomes sample `sample_A`. The current directory reader does not select files by extension, so do not mix other files into this directory.

Alternatively, provide a **headerless**, tab-separated `samples.tsv`:

```tsv
sample_A	reads/sample_A.fq
sample_B	reads/sample_B.fq
```

Paths are resolved from the working directory. Use distinct sample names; letters, digits and underscores are recommended.

## Primers and reference databases

| Input | Option | Required format |
|---|---|---|
| Primers | `-p` | FASTA with records named `forward` and `reverse` |
| Taxonomy reference | `-d` | SINTAX-formatted reference FASTA |
| Clustered reference metadata | `-b` | TSV with the legacy clustered database layout below |
| Host reference | `-g` | Optional host genome FASTA |

Example primer FASTA (replace with your experimental primers):

```fasta
>forward
AGAGTTTGATCMTGGCTCAG
>reverse
TACGGYTACCTTGTTACGACTT
```

Supply both primers in their synthesized 5′–3′ orientation. The pipeline reverse-complements the reverse primer when constructing the cutadapt linked adapter.

An illustrative SINTAX reference header is:

```text
>reference_001;tax=d:Bacteria,p:Bacillota,c:Bacilli,o:Bacillales,f:Staphylococcaceae,g:Staphylococcus,s:Staphylococcus_epidermidis;
```

Use reference files that agree on species names. The clustered TSV parser skips the first line as a header and requires at least 15 columns per data row. It reads these positions; a generic two-column taxonomy table is insufficient:

| Position (1-based) | Information read by PMCSMRT |
|---|---|
| Column 2 | Reference cluster ID |
| Column 13 | Species name used for the join |
| Third column from the end | Number of species in the cluster |
| Last column | Number of OTUs with the same species |

`num_of_sp_in_db` is computed as the number of metadata rows sharing that species name. Primers and reference databases must be supplied separately.

## Preparing legacy PacBio archives with `get_fastqs`

Already prepared sample FASTQs with valid `ccs=` annotations can be used directly. For legacy SMRTportal/SMRTlink job outputs containing `barcoded-fastqs.tgz` and `.ccs.h5` files, this optional utility extracts sample FASTQs and adds CCS counts and barcode labels.

Install its additional dependency:

```bash
conda install -n usearch -c conda-forge h5py
```

The sample key has a header and five tab-separated columns:

```tsv
PB_jobid	data_path	forward_barcode	reverse_barcode	sample_name
17094	/data/pacbio/jobs/017094	0001_Forw_bc	0002_Rev_bc	sample_A
17094	/data/pacbio/jobs/017094	0003_Forw_bc	0004_Rev_bc	sample_B
```

Each `data_path` must contain the archive and matching `.ccs.h5` files. Archive FASTQs are expected to be named `<forward_barcode>--<reverse_barcode>.fastq`. The utility reads legacy HDF5 metadata; it does not extract pass counts from modern BAM files or perform demultiplexing.

Run it from a dedicated preparation working directory, with the package parent on `PYTHONPATH`:

```bash
PYTHONPATH=/path/to/package_parent python -m pmcsmrt.get_fastqs \
    -s sample_key.tsv -o prepared_reads
```

Use a new output directory. The utility uses and removes a working directory named `raw_data`, and writes `passes` and `log.txt`. Sample names should begin with a letter and use letters, digits and underscores. Missing CCS metadata produces `ccs=nil`; supply valid numeric pass counts before running the pipeline, which cannot parse `nil`.

# Usage

## Quick start

From the package parent directory, run UNOISE:

```bash
PYTHONPATH=/path/to/package_parent python -m pmcsmrt.cli \
    -c unoise \
    -f reads \
    -p primers.fa \
    -d references/sintax_db.fa \
    -b references/clustered_db.tsv \
    -o results_unoise \
    -t 8 -V
```

UPARSE:

```bash
PYTHONPATH=/path/to/package_parent python -m pmcsmrt.cli \
    -c uparse \
    -f reads \
    -p primers.fa \
    -d references/sintax_db.fa \
    -b references/clustered_db.tsv \
    -o results_uparse \
    -t 8 -V
```

PYTHONPATH=/path/to/package_parent DADA2:

```bash
python -m pmcsmrt.cli \
    -c dada2 \
    -f reads \
    -p primers.fa \
    -d references/sintax_db.fa \
    -b references/clustered_db.tsv \
    -o results_dada2 \
    -t 8 -V
```

Omitting `-c` selects UNOISE. To use a sample list, replace `-f reads` with `-i samples.tsv`. To filter host-mapping reads, add `-g references/host.fa`. Host filtering is skipped when `-g` is omitted; a BWA index is not required.

## Command-line options

| Option | Description | Default |
|---|---|---|
| `-f`, `--foldername` | Directory containing sample FASTQs | Supply `-f` or `-i` |
| `-i`, `--samplelist` | Headerless sample-name / FASTQ-path TSV | Supply `-f` or `-i` |
| `-p`, `--primerfile` | Forward and reverse primer FASTA | Required |
| `-d`, `--taxadbfile` | SINTAX taxonomy reference FASTA | Required |
| `-b`, `--dbclusteredfile` | Clustered reference metadata TSV | Required |
| `-g`, `--host_db` | Host genome FASTA | No host filtering |
| `-o`, `--output` | Output directory | `pmcsmrt_output` |
| `-c`, `--cluster-method` | `unoise`, `uparse` or `dada2` | `unoise` |
| `-t`, `--threads` | Requested thread count | Available CPUs, capped at 32 |
| `-s`, `--ccsvalue` | Minimum CCS passes | `5` |
| `-n`, `--lengthmin` | Minimum trimmed length (bp) | `500` |
| `-x`, `--lengthmax` | Maximum trimmed length (bp) | `2000` |
| `-e`, `--eevalue` | Maximum expected errors (EE) | `1.0`; DADA2: `5.0` |
| `-m`, `--min_identity` | Read mapping identity, as a fraction | `0.99`; UPARSE: `0.97` |
| `-a`, `--abundance` | DADA2 chimera minimum parental fold abundance | `2` |
| `-w`, `--split_switch` | Export per-feature reads: `yes` or `no` | `no` |
| `-j`, `--splitmethod` | Export reads `before` or `after` inference filtering | `after` |
| `-V`, `--verbose` | Keep intermediate files and print debug messages | Off |
| `-v`, `--version` | Print pipeline version | Current version: `3.0.0` |
| `-h`, `--help` | Show help | — |

Explicit `-e` and `-m` values override method-specific defaults. Length boundaries are inclusive. **`-m` changes mapping identity, not the UPARSE clustering threshold.** `-a` affects DADA2 chimera removal only, not UNOISE or UPARSE minimum abundance. The SINTAX command uses a fixed cutoff of `0.8`.

## Export reads for individual features

For example, export UNOISE reads after inference filtering:

```bash
python -m pmcsmrt.cli \
    -c unoise -f reads -p primers.fa \
    -d references/sintax_db.fa -b references/clustered_db.tsv \
    -o results_unoise_split -t 8 -w yes -j after -V
```

`-j before` exports mapped reads passing CCS, primer, host and length requirements before EE filtering. `-j after` exports mapped inference-filtered reads; DADA2 uses the actual R-filtered records. Output directories are `split_zotus`, `split_otus` or `split_asvs`, containing `<feature>_<stage>_subset.fasta` and `<feature>_<stage>_headers.tsv`.

These exports use USEARCH mapping assignments. For DADA2, they are mapping-based subsets and should not be used to reconstruct the native abundance table.

# Outputs

## Main files

| File | Description |
|---|---|
| `post_zOTU.fa` / `post_OTU.fa` / `post_ASV.fa` | Representative sequences for the selected method |
| `post_zOTU_table.tsv` / `post_OTU_table.tsv` / `post_ASV_table.tsv` | Feature-by-sample integer abundance table |
| `post_final_results.tsv` | Feature counts, total size, taxonomy and confidence, reference matches and cluster metadata |
| `post_taxa_collapse.tsv` | Species counts, relative abundances (%) and feature membership |
| `post_each_step_report.tsv` | Per-sample filtering and assignment statistics |
| `post_general_report.txt` | Feature and assigned-read totals and USEARCH mapping diagnostics |
| `post_usearch_glob_results.tsv` | Headerless query ID / target ID read-to-feature mapping |
| `pre_all_reads_info.tsv` | Read metrics before and after primer trimming |
| `pre_all_bc_reads_info_with_taxa.tsv` | Read metrics joined to read-level taxonomy |
| `all_bc_reads.sintax` | Raw SINTAX results for input reads |
| `pre_trimmed_and_oriented.fq` | Cutadapt output before inference filtering |
| `mcsmrt_files/` | Intermediate files retained with `-V` |

Only the representative FASTA and abundance table matching the selected method are produced. Zero-count sample columns are retained. DADA2 ASV IDs are assigned in lexical sequence order and are deterministic for the same set of sequences.

## Abundance and taxonomy tables

Example from `test_results/pmcsmrt_unoise/post_zOTU_table.tsv`:

```tsv
#zOTU ID	A_22_accu	A_22_gotaq	A_hum_22_accu	A_hum_22_gotaq
Zotu22	131	116	127	79
Zotu11	265	244	255	168
```

The final annotated table adds:

| Column or group | Meaning |
|---|---|
| Feature ID | Identifier matching the representative FASTA |
| `size` | Sum of the feature's counts across sample columns |
| Sample columns | Assigned read counts per sample |
| `domain` through `species` | SINTAX classifications at seven ranks |
| `domain_conf` through `species_conf` | Confidence for each rank |
| `blast_query` | Matched reference sequence header |
| `blast_percent_identity` | Representative-to-reference alignment identity (%) |
| `blast_alignment_length` | Alignment length (bp) |
| `ncbi_avglinkage_otu_id` | Cluster ID from reference metadata |
| `num_of_sp_in_cluster` | Species count in the reference cluster |
| `num_of_sp_in_db` | Metadata row count for the assigned species |
| `num_of_otus_with_same_sp` | Same-species OTU count from reference metadata |

The `blast_*` names are retained for compatibility, but the values come from **USEARCH `-usearch_global`**, not `blastn`. Reference matching uses identity `0.97`, independently of the read-mapping `-m` setting.

Species summaries contain counts, corresponding `*_rel` columns and `zOTU_member`, `OTU_member` or `ASV_member`. Relative abundance is `100 × species count / total count` for each corresponding count column. Empty samples receive `0`. Unassigned species are grouped under `NA` and included in the denominator.

For example, the supplied UNOISE output groups `Zotu22;Zotu19;Zotu12;Zotu29;Zotu25` under `Staphylococcus_epidermidis`, with total `size` of `2617` and `size_rel` of approximately `7.02%`.

## Read-level information

| Column | Meaning |
|---|---|
| `read_name`, `basename` | Read identifier and sample label |
| `ccs` | Pass count from the FASTQ header |
| `ee_pretrim`, `ee_posttrim` | Expected errors before and after trimming |
| `length_pretrim`, `length_posttrim` | Read lengths before and after trimming |
| `host_map` | Whether the read was marked as host-mapping |
| `read_orientation` | `+`, `-` or `NA` after primer matching |
| `primer_note` | Primer matching status |
| `Domain` through `Species`, plus `*_conf` | Read taxonomy and confidence in the joined table |

Read-level taxonomy is assigned to input reads, including reads later excluded by filtering. It is separate from representative-based abundance reporting.

## Filtering report

Example from the supplied UNOISE result:

```tsv
sample	method	total	ccs_filt	primer_filt	host_filt	size_filt	ee_filt	unique_filtered	denoised_reads	non_chimeric_reads	assigned_reads
A_22_accu	unoise	14865	14865	13346	13346	13342	5712	3328	NA	NA	11185
```

| Column | Meaning |
|---|---|
| `total` | Input reads |
| `ccs_filt` | Reads passing the CCS threshold |
| `primer_filt` | CCS-passing reads with primer matches |
| `host_filt` | Reads remaining after optional host exclusion |
| `size_filt` | Reads passing trimmed length limits |
| `ee_filt` | Inference-filtered reads; for DADA2, the complete R quality/length filter result |
| `unique_filtered` | Distinct filtered sequences per sample, not a read count |
| `denoised_reads` | Native DADA2 denoised reads; `NA` for USEARCH methods |
| `non_chimeric_reads` | Native DADA2 reads after chimera removal; `NA` for USEARCH methods |
| `assigned_reads` | Sum of the sample's final abundance counts |

The shared preprocessing columns through `size_filt` are cumulative read counts. DADA2's `size_filt` is the pre-R length check, while `ee_filt` reflects R filtering, which can also truncate low-quality tails.

## Additional DADA2 outputs

| File | Contents |
|---|---|
| `post_dada2_stats.tsv` | `sample`, `input`, `filtered`, `denoised`, `non_chimeric` |
| `post_dada2_versions.txt` | R and DADA2 versions used by the run |
| `post_base_transition_stats.tsv` | `transition`, `quality`, `count`, `estimated`, `observed`, `nominal` |

DADA2's `input` is the input to R after shared CCS, primer and host filtering, rather than the original FASTQ count. Transition statistics describe the learned error model by base transition and quality score. With `-V`, the R log, error-model RDS, manifest and auxiliary mapping table are retained in `mcsmrt_files/`.

## Supplied example results

The [test_results](test_results) directory contains outputs from all three methods:

| Result directory | Features | Assigned reads in the abundance table |
|---|---:|---:|
| `pmcsmrt_unoise` | 57 zOTUs | 37,279 |
| `pmcsmrt_uparse` | 19 OTUs | 41,802 |
| `pmcsmrt_dada2` | 52 ASVs | 25,826 |

These are output examples, not a benchmark or evidence of higher accuracy. Methods use different inference, filtering and counting procedures. Full execution parameters and external tool versions cannot be inferred from these totals alone.

# Tips and Notes

## Interpreting counts across methods

For UNOISE and UPARSE, EE-filtered reads are used to infer representatives. Abundance is then calculated by mapping reads passing CCS, primer, host and length requirements, including reads that failed EE filtering. Consequently, `assigned_reads` can exceed `ee_filt`, as in the example above.

DADA2 counts come from its native non-chimeric ASV table. Additional USEARCH mapping is used for diagnostics and extraction and does not replace these counts. Unmapped USEARCH reads in `post_general_report.txt` are not a measure of native DADA2 read loss.

## Taxonomic confidence

Available ranks are retained independently; missing ranks are `NA`. Although the SINTAX command uses cutoff `0.8`, final tables parse its raw assignment field and can retain labels with confidence below `0.8`. Review the `*_conf` columns before interpreting species assignments. The species collapse table does not apply an additional confidence filter.

Species-level interpretation depends on the database, confidence and how well the 16S sequence distinguishes related species. Cluster metadata provides context about reference ambiguity.

## DADA2 processing

DADA2 receives cutadapt-trimmed and oriented reads and does not remove primers again. It uses PacBio error learning, independent sample inference (`pool=FALSE`) and consensus chimera removal. Filtering uses `truncLen=0`, `trimLeft=0`, `truncQ=2`, `minQ=3`, `maxN=0` and `rm.phix=FALSE`, together with the requested length and EE limits.

Primer matching follows this pipeline's cutadapt settings, so results need not exactly reproduce a QIIME2 CCS workflow. On native Windows, DADA2 FASTQ filtering runs in one process; the requested threads are passed to error learning and denoising. The supplied Conda YAML targets Linux/WSL.

## Output directories and diagnostics

Use a separate output directory for each run or method. Reusing a directory can leave files from previous runs alongside current outputs. Keep generated results outside the source directory.

Use `-V` to retain intermediates. Backend diagnostics are copied to `mcsmrt_files/` even when feature inference or downstream processing fails. Failures before backend dispatch may occur before that copy is made. If no reads or features remain, check CCS annotations, primers, thresholds and input format.

## Validation status

The supplied results demonstrate the output layout. The existing [validation record](VALIDATION.md) separately documents synthetic Python tests for dispatch, filtering, taxonomy joins, native DADA2 count preservation and extraction. It records unavailable external tools and skipped native checks in that validation environment. Example outputs do not establish that every backend has been independently validated in the current checkout.

# Citation

If you use this pipeline in published work, please cite the original MCSMRT publication:

Earl, J. P., Adappa, N. D., Krol, J., et al. (2018). Species-level bacterial community profiling of the healthy sinonasal microbiome using Pacific Biosciences sequencing of full-length 16S rRNA genes. *Microbiome*, **6**, 190. [doi: 10.1186/s40168-018-0569-2](https://doi.org/10.1186/s40168-018-0569-2).

Please also acknowledge the inference and classification tools used in your analysis and report their versions and database sources. If you would also like to acknowledge this modified implementation, a repository reference would be appreciated:

Genglin Guo. *PMCSMRT*. [GitHub repository](https://github.com/guogenglin/pmcsmrt).

# Contact

Maintainer of this Python version: **Genglin Guo**.

Questions, bug reports and suggestions are welcome through [GitHub Issues](https://github.com/guogenglin/pmcsmrt/issues) or by email at [gg599@drexel.edu](mailto:gg599@drexel.edu).
