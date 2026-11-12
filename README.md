# HGTlong

Read-resolved candidate shared genomic blocks from ONT alignments. **Version 0.1.3rc3: research release candidate.** Empirical evaluation is currently limited to the original research sample. Calls do not establish transfer direction, mechanism or age.

This README is the single unified user and release guide.

Author: Wei-Hua Chen (weihuachen@hust.edu.cn). Repository: https://github.com/whchenlab/HGTlong .

Copyright 2026 Weihua Chen. Licensed under [Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)](LICENSE). Commercial use is not granted under this license. External tools retain their own licenses.

## Contents

- [Demo and Linux real-data workflow](#demo-and-linux-real-data-workflow)
- [CPU cores and parallel execution](#cpu-cores-and-parallel-execution)
- [Calling semantics and implementation](#calling-semantics-and-implementation)
- [Development and verification](#development-and-verification)
- [Release scope](#release-scope)
- [Changes](#changes)

## Demo and Linux real-data workflow

This guide covers the supported HGTlong read-first workflow from MAG FASTA files and ONT FASTQ reads through read calling and block screening. It also explains the external annotation and ANI inputs needed before declaring a final study set. Version 0.1.3rc3 does not automate whole-MAG ANI or BLAST discovery; those handoffs are explicit below. Do not mistake read_supported calls for fully screened transfers.

## 1. Choose the appropriate workflow

**Demo: only Python >=3.9 and HGTlong are needed. No external bioinformatics tools or databases are required.** From the unpacked repository after installation:

```bash
python demo/run_demo.py --outdir demo-output
```

The demo uses curated synthetic alignments, not real ONT mapping. It checks the software, not biological accuracy or ANI. Open demo-output/read-first/report.html in a browser.

Real data require minimap2 for ONT/reference-block alignment and samtools for alignment processing. Heavy computation belongs in scheduler allocations on compute nodes, not login or bridge hosts.

## 2. Install on Linux

Use an existing Conda installation, or install Miniforge appropriate for the server architecture from https://github.com/conda-forge/miniforge . No root access is required for a user-owned environment. Follow the official installer instructions; activate Conda in the shell before proceeding.

From the unpacked HGTlong repository:

```bash
conda create -n hgtlong --override-channels   -c conda-forge -c bioconda --strict-channel-priority   python=3.11 pip minimap2=2.31 samtools=1.24
conda activate hgtlong
python -m pip install .
```

Alternatively, `CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml` installs these dependencies and the local package; run it from the repository root. Mamba can replace Conda. Keep channel order conda-forge then bioconda. On offline clusters, arrange an environment/container through the cluster administrators; the Python wheel alone does not provide external executables.

```bash
hgtlong --version
minimap2 --version
samtools --version
python demo/run_demo.py --outdir demo-output
```

Recommended minimap2 2.31 and samtools 1.24 match the recorded the original research sample benchmark. They are reproducibility pins, not established minimum versions. The package requires Python >=3.9; local package/demo testing used 3.9. The new Linux Conda specification has not yet been independently solved and tested as a release environment. Different tool versions can change mappings and calls. Pinning versions alone does not reproduce an experiment without matching inputs and parameters.

## 3. Prepare inputs

- One MAG per FASTA file in a single directory: .fa, .fna, .fasta, or their .gz forms. Use stable unique filenames without whitespace; the filename stem becomes the MAG ID.
- ONT reads from the sample corresponding to the MAG panel. Record sequencing chemistry, basecaller/model, filtering and any host-read removal when available.
- Stable, unique read names. Reused names across concatenated files can merge unrelated molecules. Do not rename individual alignments after mapping.
- Adequate scratch space for the BAM and query collation. Space/memory needs depend strongly on read volume and reference size; the example allocation below is not a sizing guarantee.

Use an absolute project directory and paths without whitespace for the simple shell examples:

```bash
export PROJECT=/path/to/project
mkdir -p "$PROJECT/input/MAGs" "$PROJECT/run" "$PROJECT/logs"
# Put MAG FASTAs in input/MAGs and ONT reads in input/reads.fastq.gz.
```

MAG completeness/contamination assessment is upstream of HGTlong. A limited panel can omit alternative placements; HGTlong cannot validate unsampled genomes.

## 4. Submit mapping and read-first calling

The supplied `workflows/read_first.slurm` is a portable Slurm example. Adjust memory, time, partition/account requirements and paths for your cluster. Activate the installed environment before submitting so its PATH is inherited:

```bash
conda activate hgtlong
export PROJECT=/path/to/project
mkdir -p "$PROJECT/logs"
sbatch --chdir="$PROJECT"   --output="$PROJECT/logs/hgtlong-%j.out"   --error="$PROJECT/logs/hgtlong-%j.err"   /path/to/HGTlong/workflows/read_first.slurm
```

For other schedulers, retain the shell commands but replace the scheduler directives and resource variables with the local equivalents. Do not execute the job script on a login/bridge host. In this project's mn01-connected cluster, submit via its scheduler to permitted compute nodes and explicitly exclude mn01.

The job builds a namespaced panel, maps all reads, retains secondary/supplementary alignments, collates alignments by query name, extracts evidence and calls events. It uses `set -euo pipefail` to stop on pipeline failures. Use a fresh run directory for reruns because build-panel refuses to overwrite a panel. Preserve panel/panel.tsv: original contig identifiers change during namespacing.

Mapping flags `-p 0 -N 100 -Y` retain alternatives needed for ambiguity assessment. Avoid MAPQ filtering before extraction. A pre-existing BAM is usable only if generated against exactly this panel with suitable alternative-alignments retained; otherwise remap. Ordinary coordinate-sorted SAM is not valid extraction input. Do not confuse collation with coordinate sorting.

## 5. Inspect read-supported candidates

Outputs in run/calls:

- events.tsv: all candidate statuses, coordinates and boundary support counts.
- event_evidence.tsv.gz: molecule memberships, identities and ambiguity information.
- report.html: standalone interactive browser; no web server required.
- run.json: parameters, hashes, version and status counts.

A read_supported status requires qualifying support at all four boundaries. partial/ambiguous calls are not final events. All event coordinates are zero-based, half-open. A/B is a label, not transfer direction. Preserve nonpassing calls for audit.

The report produced by call describes the read stage. Later table screens do not automatically regenerate it or change its displayed status. The final machine-readable screened/ANI tables are authoritative for those later stages.

## 6. Supply annotation masks and screen reference blocks

Prepare annotations against run/panel/panel.fa with an appropriate external annotation workflow. For the study-like rRNA exclusion, supply rRNA coordinates in BED format. Record tool/database versions and parameters. COG/mobile-element interpretation is optional and is not required by the core caller.

Required files:

- annotations/masks.bed: tab-separated contig, start, end, optional class. Coordinates zero-based, half-open, with no header (comments start with #). Contigs must match the namespaced panel.
- annotations/screened_contigs.txt: one namespaced contig per line, listing only contigs actually assessed successfully for the supplied mask classes.

Example format only (not valid substitute annotations):

```text
MAG_A__contig1<TAB>100<TAB>1600<TAB>rRNA
```

Replace <TAB> with real tabs. An empty BED is acceptable only when annotation genuinely completed without hits. Do not label every contig screened merely to bypass annotation. Do not indiscriminately mask mobile elements: that would exclude biologically interesting HGT loci. Any supplied mask overlap is a rejection at the screen stage.

Run the following within a compute allocation, from $PROJECT:

```bash
set -euo pipefail
hgtlong qc --events run/calls/events.tsv   --mask-bed annotations/masks.bed   --screened-contigs annotations/screened_contigs.txt   --out run/calls/events_qc.tsv
hgtlong sequence-qc --events run/calls/events_qc.tsv   --reference run/panel/panel.fa --minimap2 minimap2   --out run/calls/events_sequence_qc.tsv
hgtlong screen --events run/calls/events_sequence_qc.tsv   --min-block-identity 0.99 --min-block-coverage 0.99   --min-boundary-reads 3   --out run/calls/all_screened.tsv   --candidates run/calls/screened_candidates.tsv
```

Explicit coverage 0.99 is essential for the later the original research sample criteria: the historical CLI default is 0.95. Identity and coverage thresholds are fractions. Both blocks must meet coverage. Screening annotates all failures in screen_fail_reasons and never uses HDMI labels.

## 7. Whole-MAG ANI: required external handoff for study-like final calls

screened_candidates.tsv is **pre-ANI**, not a final study set. Extract unique magA/magB pairs and compute bidirectional whole-genome ANIb using the original full MAG FASTAs linked by panel_manifest.json. Do not compute ANI only over the candidate HGT blocks. Archive the exact ANI implementation, BLAST version, fragment/hit settings and input checksums.

The original benchmark used BLAST+ 2.13.0 for its separate ANI workflow, whereas native BLAST discovery used 2.17.0. Keep them in separate environments if reproducing that setup. Installing BLAST alone is not an ANIb implementation, and FastANI is not a drop-in reproduction of the study's ANIb calculation.

For each event retain only pairs with valid ANI **below 95% in both directions**. Missing/uncomputed ANI is unresolved and cannot pass. Preserve aligned fractions; low aligned fraction is an interpretation flag, not an automatic extra rejection in the original study. Join results by MAG IDs, retain both directional values, and record an explicit final retention reason for every candidate. This release supplies no automated ANI calculation or join command: obtaining those results is an external prerequisite, not something screen already performs.

## 8. HGTlong-blast scope

call-blocks accepts precomputed fixed block candidates and the same extracted evidence:

```bash
hgtlong call-blocks --blocks blocks.tsv   --evidence run/evidence.tsv.gz --outdir run/block-calls
```

Required TSV columns: event_id, magA, ctgA, a0, a1, magB, ctgB, b0, b1, orientation. See demo/blocks.tsv. IDs must be unique, endpoints ordered, orientation + or -, and the canonical A/B ordering must match extracted evidence. Reverse-strand endpoints remain increasing; strand is represented by orientation. Raw BLAST output is not this schema.

A full BLAST-based MAG search and conversion to this schema are not included in this release. Do not describe fixed-block validation alone as a complete HGTlong-blast rerun. Apply the same annotation, block and ANI checks afterwards.

## 9. Preserve provenance and troubleshoot

Record the environment and versions after installation:

```bash
conda list --explicit > run/conda-explicit.txt
python -m pip freeze > run/python-packages.txt
hgtlong --version > run/hgtlong.version.txt
minimap2 --version > run/minimap2.version.txt
samtools --version > run/samtools.version.txt
```

The explicit Conda export is platform-specific. Preserve job script/logs, input checksums, reference-panel manifest, evidence manifest, all filter tables and ANI/annotation provenance. The exported paths may need redaction before public sharing.

- Command not found: activate the same environment in the batch job; check command -v hgtlong/minimap2/samtools.
- Query-name collation error: collate the BAM as in the supplied job; retain all alignments of each molecule together.
- Missing NM tags: regenerate suitable alignments; absence must not be interpreted as perfect identity.
- Panel outputs already exist: use a fresh run directory instead of overwriting a reference used by existing BAMs.
- No calls: inspect extraction stats, input reference names, mapping quality, depth and all four boundary counts before changing thresholds.
- Annotation not_fully_assessed: confirm successful assessment and namespaced IDs, rather than forcing a pass.
- Out of memory/time/disk: inspect scheduler logs and request suitable resources; queue a revised job on a compute node.

Validation is currently limited to the original research sample plus software tests and the synthetic demo. A final retained locus remains a candidate transfer; its age, mechanism and direction require separate evidence.

## CPU cores and parallel execution

**HGTlong itself supports multiple CPU cores during evidence extraction.** It uses worker processes, not Python threads. Other stages have different parallelism; requesting CPUs from a scheduler does not automatically make every command use them.

| Stage | Multiple CPU cores? | Control and current behavior |
|---|---|---|
| Build MAG panel (`build-panel`) | No | Single-process FASTA parsing, namespacing and checksumming. |
| ONT mapping (external minimap2) | Yes | `-t N`; supplied Slurm example uses `-t 14` within a 16-CPU allocation. |
| BAM conversion/collation (external samtools) | Where supported | `-@ N` enables additional threads for supported operations, notably compression/I/O; it does not parallelize every operation. The example omits `-@`. |
| **Read-evidence extraction (`extract`)** | **Yes—multiprocessing** | `--workers N`; example uses 8. Default is up to 8 workers under Slurm, reserving two allocated CPUs, and 1 outside Slurm. Parsing/input and output remain potential bottlenecks. |
| Read-first event calling (`call`) | No | Currently single-process; no worker option. |
| Fixed-block validation (`call-blocks`) | No | Currently single-process; no worker option. |
| Annotation-mask overlap (`qc`) | No | Single-process. This reads supplied masks; it does not run annotation software. |
| Reference-block alignment (`sequence-qc`) | No in current implementation | Processes events sequentially and invokes minimap2 with `-t 1` for each event. No CLI parallelism option. |
| Final block screening (`screen`) | No | Single-process table filtering. Does not compute ANI. |
| HDMI concordance (`compare`) | No | Single-process comparison. |
| HTML generation | No | Single-process reporting during calling; browser rendering is separate. |
| External BLAST discovery, annotation and ANI | Tool/workflow-dependent | Not orchestrated by this release. Configure threads/processes in the chosen external workflow and respect scheduler allocation. |

For example, `hgtlong extract --sam alignments.sam --panel panel.tsv --workers 8 --out evidence.tsv.gz` uses eight worker processes plus coordinating/input-output work. More workers may increase memory use and need not improve performance if input or output is limiting. Do not set worker counts above the resources allocated to the job.

The provided 16-CPU Slurm job runs mapping first, then the collation/extraction pipeline, then event calling. Mapping and extraction are **separate phases**; their CPU counts are not added together. Within each shell pipeline, however, commands run concurrently, so leave capacity for samtools and input/output. Adding samtools `-@` threads requires revisiting that total budget.

The demo deliberately sets `--workers 1`: eleven synthetic reads are too small to benefit from multiprocessing. Real-data extraction can use multiple cores even though the later caller is single-process.

## Calling semantics and implementation

All event coordinates are zero-based, half-open; SAM positions are one-based. Event lengths are end minus start. A/B ordering is a canonical label, not donor/recipient direction.

Extraction requires an intact resident alignment covering at least 90% of query with identity at least 99%, MAPQ 20–254, and no gap above 500 bp. Competing placements including MAPQ-zero alternatives affect uniqueness. Partner fragments require at least 200 bp, at least 95% identity, no gap above 500 bp, and total query coverage at the partner MAG no greater than 75%. Ambiguous partner placements do not count as unambiguous support. NM tags are required to assess identity.

Read-first calling seeds resolved blocks from fully bracketed evidence and clusters compatible endpoints. Default boundary tolerance is 50 bp. Read support requires at least three distinct read names at each of four boundaries, resident reads in both genomic contexts, 100-bp anchors with at least 90% actual aligned bases per flank, and a block length ratio of at least 0.9. Distinct read names approximate distinct molecules. Fixed-block validation uses supplied boundaries with the same support requirements; it does not require a read-first full-block seed.

Outputs include events.tsv, event_evidence.tsv.gz, report.html and run.json. Status read_supported describes boundary evidence; partial and ambiguous calls are not final transfers. Support columns n_a0, n_a1, n_b0, n_b1 are distinct qualifying read counts. Identity and coverage values use fractions, not percentages. HTML reports retain up to 100 displayed molecules per event; the evidence table preserves all memberships.

### Screening implementation notes

The qc command needs an externally prepared BED mask and an explicit list of assessed contigs. Missing annotation must not be represented as a complete negative screen. sequence-qc compares called reference blocks using minimap2. screen requires read support, annotation_screen_pass, both_blocks_covered, and requested identity/coverage thresholds.

The historical screen default is 0.95 block coverage; to match the later the original research sample block screen explicitly use:

```
hgtlong screen --events annotated_sequence_qc.tsv --min-block-identity 0.99 --min-block-coverage 0.99 --min-boundary-reads 3 --out all_screened.tsv --candidates screened.tsv
```

This command does not calculate or enforce whole-MAG ANI. The study final sets additionally required ANI below 95% in both directions; missing ANI did not pass. Users must perform and record this separately. Synthetic demo genomes are unsuitable for benchmarking whole-genome ANI. Do not describe demo calls or screen output alone as reproducing the 130/110 study final sets.

This release contains no phylogenetic reconciliation, evolutionary dating, transfer-direction inference or proof of transfer. Genome incompleteness, repeats, contamination, mapping ambiguity and unsampled organisms remain relevant limitations.

## Development and verification

Run the tests from the repository root after installing HGTlong:

```bash
python -m unittest discover -s tests -v
```

An offline Python-only install can use the supplied wheel:

```bash
python -m pip install --no-index /path/to/hgtlong-0.1.3rc3-py3-none-any.whl
```


2026-11-12, macOS, Python 3.9.

- Built Python wheel and source distribution using setuptools.
- Installed wheel with --no-index into a fresh virtual environment.
- Ran the installed CLI from outside the source package.
- Both demo modes reproduced expected coordinates and 4/4/4/4 boundary support.
- 60 unittest tests passed against the installed package, including the demo integration test.
- Linux/Python 3.11/3.12 CI is configured but has not run yet.
- External mapping, BLAST discovery, annotation and ANI are outside this demonstration.

## Release scope

Version 0.1.3rc3 is a research prerelease.

Included: installable standard-library Python caller and reporting modules; read-first CLI; newly exposed fixed-block validation CLI; synthetic FASTA/FASTQ/SAM fixture; demo runner with expected-result assertions; regression tests; selected annotation utility used by a test; portable dependency specification; packaging and CI configuration.

Excluded: the original research sample reads/MAGs/results, third-party HDMI source, server configuration, credentials, large databases, cached environments and study-specific orchestration scripts. This is the reusable software/demo bundle, not the full manuscript reproduction archive. Mapping, BLAST discovery, ANI and annotation databases are not bundled.


Version 0.1.3rc3 differs from the internal 0.1.2 snapshot in packaging/docs, generic report labels, the call-blocks CLI wrapper, tooltip presentation, and optional post-call block-length filtering. The read-first discovery/filtering algorithms are unchanged. Existing historical screen defaults are preserved and documented explicitly.

## Changes

### 0.1.3rc3

- Prepare a portable review bundle from the internal 0.1.2 code.
- Add a fixed-block validation CLI around the existing HGTlong-blast validator.
- Add synthetic eleven-read demo, expected results and installation instructions.
- Replace historical sample-specific introduction with explicit release scope.
- Preserve existing biological thresholds and algorithms.

- Document Linux real-data workflow and scheduler template; align dependency pins with recorded benchmark versions (minimap2 2.31, samtools 1.24).

## hgtlong_demo and optional minimum block length

The public synthetic dataset is named **hgtlong_demo**. It is separate from the original research sample and is not an empirical validation dataset. It contains six reads spanning both boundaries, four reads spanning only one boundary (left and right in each resident context), and one null read. Each boundary has four supporting reads. One-boundary reads support only the endpoint they span; their truncated partner alignments do not redefine the block.

Use `--min-hgt-length BP` with `call`, `call-blocks`, or `screen` to require **both reference blocks** to have length at least BP. Default: no length cutoff. Exactly BP passes. This is distinct from extract's `--min-block` partner-alignment threshold (default 200 bp), which is unchanged.

```bash
hgtlong call --evidence run/evidence.tsv.gz --outdir run/calls --min-hgt-length 500
hgtlong call-blocks --blocks blocks.tsv --evidence run/evidence.tsv.gz --outdir run/block-calls --min-hgt-length 500
```

Calling applies this filter after evidence assignment, so excluded candidates do not release ambiguous reads to other candidates. Retained events and report contain only passing blocks. `length_excluded.tsv` and `length_filter.json` retain the rejected calls and cutoff/counts for audit. At the screen stage all rows remain in the audit table with `short_hgt_block` in screen_fail_reasons; only passing rows enter the candidates file.

For the Slurm template, optionally export `MIN_HGT_LENGTH=500` before submission; leave it unset for no cutoff. The demo runner accepts the same option; 1000 retains its block and 1001 excludes it.

Hovering or keyboard-focusing a read highlights its matching alignment in both A and B in magenta. Moving away, leaving focus, or pressing Escape clears the highlight. This links the same read ID, not separate molecules with similar sequences.
